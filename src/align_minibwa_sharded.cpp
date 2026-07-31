#include "align_minibwa_sharded.hpp"
#include "align_common.hpp"
#include "shard_debug.hpp"
#include "shard_progress.hpp"
#include "duckdb/common/file_system.hpp"

namespace duckdb {

// Build minibwa MiniBWAShardInfo from raw shard name/counts. Validates that
// both index files exist for each shard -- minibwa has no single-file
// magic-number probe like minimap2's mm_idx_is_idx(), so this checks file
// existence only (same advisory/TOCTOU caveat as align_minibwa.cpp's Bind).
static std::vector<MiniBWAShardInfo> BuildMiniBWAShardInfos(ClientContext &context, const std::string &table_name,
                                                            const std::string &shard_directory, FileSystem &fs) {
	auto raw_shards = ReadShardNameCounts(context, table_name);

	std::vector<MiniBWAShardInfo> shards;
	shards.reserve(raw_shards.size());

	for (const auto &raw : raw_shards) {
		MiniBWAShardInfo info;
		info.name = raw.name;
		info.read_count = raw.count;

		info.index_prefix = shard_directory;
		if (!info.index_prefix.empty() && info.index_prefix.back() != '/') {
			info.index_prefix += '/';
		}
		info.index_prefix += info.name;

		if (!fs.FileExists(info.index_prefix + ".l2b") || !fs.FileExists(info.index_prefix + ".mbw")) {
			throw BinderException("Shard index files do not exist for prefix: %s (expected %s.l2b and %s.mbw)",
			                      info.index_prefix, info.index_prefix, info.index_prefix);
		}

		shards.push_back(std::move(info));
	}

	return shards;
}

unique_ptr<FunctionData> AlignMiniBWAShardedTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                vector<LogicalType> &return_types,
                                                                vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 1) {
		throw BinderException("align_minibwa_sharded requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	auto shard_dir_param = input.named_parameters.find("shard_directory");
	if (shard_dir_param == input.named_parameters.end() || shard_dir_param->second.IsNull()) {
		throw BinderException("align_minibwa_sharded requires shard_directory parameter");
	}
	data->shard_directory = shard_dir_param->second.ToString();

	auto read_to_shard_param = input.named_parameters.find("read_to_shard");
	if (read_to_shard_param == input.named_parameters.end() || read_to_shard_param->second.IsNull()) {
		throw BinderException("align_minibwa_sharded requires read_to_shard parameter");
	}
	data->read_to_shard_table = read_to_shard_param->second.ToString();

	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.DirectoryExists(data->shard_directory)) {
		throw BinderException("Shard directory does not exist: %s", data->shard_directory);
	}

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);
	ValidateReadToShardSchema(context, data->read_to_shard_table, data->query_schema.id_type);

	data->subject_id_type = LogicalType::VARCHAR;
	data->types = GetAlignmentOutputTypes(data->query_schema.id_type, data->subject_id_type);

	ParseMiniBWAConfigParams(input.named_parameters, data->config);

	auto max_tps_param = input.named_parameters.find("max_threads_per_shard");
	if (max_tps_param != input.named_parameters.end() && !max_tps_param->second.IsNull()) {
		auto val = max_tps_param->second.GetValue<int32_t>();
		if (val < 1 || val > 64) {
			throw BinderException("max_threads_per_shard must be between 1 and 64 (got %d)", val);
		}
		data->max_threads_per_shard = static_cast<idx_t>(val);
	}

	auto debug_param = input.named_parameters.find("debug");
	if (debug_param != input.named_parameters.end() && !debug_param->second.IsNull()) {
		data->debug = debug_param->second.GetValue<bool>();
	}

	auto progress_param = input.named_parameters.find("progress");
	if (progress_param != input.named_parameters.end() && !progress_param->second.IsNull()) {
		data->progress = progress_param->second.GetValue<bool>();
	}

	auto include_shard_param = input.named_parameters.find("include_shard_name");
	if (include_shard_param != input.named_parameters.end() && !include_shard_param->second.IsNull()) {
		data->include_shard_name = include_shard_param->second.GetValue<bool>();
	}

	data->shards = BuildMiniBWAShardInfos(context, data->read_to_shard_table, data->shard_directory, fs);

	if (data->include_shard_name) {
		data->names.emplace_back("shard_name");
		data->types.emplace_back(LogicalType::VARCHAR);
	}

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMiniBWAShardedTableFunction::InitGlobal(ClientContext &context,
                                                                                  TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->shard_count = data.shards.size();
	gstate->max_threads_per_shard = data.max_threads_per_shard;

	idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	idx_t derived = (db_threads + data.max_threads_per_shard - 1) / data.max_threads_per_shard;
	gstate->max_active_shards = std::max<idx_t>(1, std::min(derived, gstate->shard_count));
	gstate->debug = data.debug;
	gstate->progress = data.progress;
	gstate->start_time = std::chrono::steady_clock::now();
	idx_t total = 0;
	for (const auto &shard : data.shards) {
		total += shard.read_count;
	}
	gstate->total_associations.store(total, std::memory_order_relaxed);
	SHARD_DBG_MEM(*gstate, "InitGlobal: shards=%zu db_threads=%zu max_tps=%zu max_active=%zu MaxThreads=%zu",
	              static_cast<size_t>(gstate->shard_count), static_cast<size_t>(db_threads),
	              static_cast<size_t>(gstate->max_threads_per_shard), static_cast<size_t>(gstate->max_active_shards),
	              static_cast<size_t>(gstate->MaxThreads()));
	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMiniBWAShardedTableFunction::InitLocal(ExecutionContext &context,
                                                                                TableFunctionInitInput &input,
                                                                                GlobalTableFunctionState *global_state) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	lstate->aligner = std::make_unique<miint::MiniBWAAligner>(data.config);
	return lstate;
}

static constexpr idx_t SHARD_READ_BATCH_SIZE = 2048;

std::shared_ptr<MiniBWAActiveShard> AlignMiniBWAShardedTableFunction::ClaimWork(ClientContext &context,
                                                                                GlobalState &gstate,
                                                                                const Data &bind_data,
                                                                                LocalState &lstate) {
	std::shared_ptr<MiniBWAActiveShard> active;
	idx_t shard_idx;

	SHARD_DBG(gstate, "ClaimWork: enter");

	{
		std::unique_lock<std::mutex> lock(gstate.lock);

		while (true) {
			// Phase 1: Try to join an existing active shard with capacity
			for (auto &shard : gstate.active_shards) {
				if (shard->ready.load(std::memory_order_acquire) && !shard->exhausted.load(std::memory_order_acquire) &&
				    shard->active_workers.load(std::memory_order_acquire) < gstate.max_threads_per_shard) {
					shard->active_workers.fetch_add(1, std::memory_order_acq_rel);
					auto &info = bind_data.shards[shard->shard_idx];
					SHARD_DBG(gstate, "ClaimWork: JOIN shard %zu '%s' (workers=%zu)",
					          static_cast<size_t>(shard->shard_idx), info.name.c_str(),
					          static_cast<size_t>(shard->active_workers.load(std::memory_order_relaxed)));
					return shard;
				}
			}

			// Phase 2: Try to claim a new shard if under the active shard limit
			bool has_unclaimed = gstate.next_shard_idx < bind_data.shards.size();
			idx_t max_active = gstate.max_active_shards;
			if (has_unclaimed && bind_data.shards[gstate.next_shard_idx].read_count <= SHARD_READ_BATCH_SIZE) {
				max_active = gstate.max_active_shards * gstate.max_threads_per_shard;
			}
			if (has_unclaimed && gstate.active_shards.size() < max_active) {
				shard_idx = gstate.next_shard_idx++;
				active = std::make_shared<MiniBWAActiveShard>();
				active->shard_idx = shard_idx;
				active->active_workers.store(1, std::memory_order_release);
				gstate.active_shards.push_back(active);
				SHARD_DBG(gstate, "ClaimWork: NEW shard %zu '%s' (active_shards=%zu)", static_cast<size_t>(shard_idx),
				          bind_data.shards[shard_idx].name.c_str(), static_cast<size_t>(gstate.active_shards.size()));
				break;
			}

			// Phase 3: Can't join or start - check if waiting is worthwhile
			bool any_not_exhausted = false;
			for (auto &shard : gstate.active_shards) {
				if (!shard->exhausted.load(std::memory_order_acquire)) {
					any_not_exhausted = true;
					break;
				}
			}

			if (!has_unclaimed && !any_not_exhausted) {
				SHARD_DBG(gstate, "ClaimWork: DONE (no more work)");
				return nullptr;
			}

			SHARD_DBG(gstate, "ClaimWork: WAIT (active=%zu, unclaimed=%s, any_alive=%s)",
			          static_cast<size_t>(gstate.active_shards.size()), has_unclaimed ? "yes" : "no",
			          any_not_exhausted ? "yes" : "no");
			gstate.cv.wait(lock);
		}
	}
	// Lock released

	// Phase 4: Load index OUTSIDE lock
	auto &shard_info = bind_data.shards[shard_idx];
	SHARD_DBG(gstate, "ClaimWork: LOADING index '%s'", shard_info.index_prefix.c_str());
	auto load_start = std::chrono::steady_clock::now();
	try {
		auto shared_idx = std::make_shared<miint::SharedMiniBWAIndex>(shard_info.index_prefix, bind_data.config);
		active->index = std::move(shared_idx);
	} catch (...) {
		SHARD_DBG(gstate, "ClaimWork: LOAD FAILED shard %zu", static_cast<size_t>(shard_idx));
		active->exhausted.store(true, std::memory_order_release);
		active->active_workers.fetch_sub(1, std::memory_order_acq_rel);
		{
			std::lock_guard<std::mutex> guard(gstate.lock);
			auto &shards = gstate.active_shards;
			shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
		}
		gstate.cv.notify_all();
		throw;
	}
	auto load_ms =
	    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - load_start).count();
	SHARD_DBG_MEM(gstate, "ClaimWork: LOADED index shard %zu '%s' in %ldms", static_cast<size_t>(shard_idx),
	              shard_info.name.c_str(), static_cast<long>(load_ms));

	// Phase 4b+4c: Materialize read IDs and pre-fetch sequences.
	idx_t seq_count;
	try {
		SHARD_DBG(gstate, "ClaimWork: MATERIALIZING IDs for shard %zu '%s'", static_cast<size_t>(shard_idx),
		          shard_info.name.c_str());
		auto ids_start = std::chrono::steady_clock::now();
		auto shard_read_ids =
		    ReadShardIds(context, bind_data.read_to_shard_table, shard_info.name, bind_data.query_schema.id_type);
		auto ids_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - ids_start).count();
		idx_t id_count = shard_read_ids.size();
		SHARD_DBG_MEM(gstate, "ClaimWork: MATERIALIZED %zu IDs for shard %zu in %ldms", static_cast<size_t>(id_count),
		              static_cast<size_t>(shard_idx), static_cast<long>(ids_ms));

		SHARD_DBG(gstate, "ClaimWork: PRE-FETCHING sequences for shard %zu '%s'", static_cast<size_t>(shard_idx),
		          shard_info.name.c_str());
		auto fetch_start = std::chrono::steady_clock::now();
		ReadBatchByIds(context, bind_data.query_table, bind_data.query_schema, shard_read_ids, 0, shard_read_ids.size(),
		               active->shard_sequences);
		auto fetch_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - fetch_start)
		        .count();
		seq_count = active->shard_sequences.size();
		active->batch_size = std::min(SHARD_READ_BATCH_SIZE, std::max<idx_t>(1, seq_count));
		SHARD_DBG_MEM(gstate, "ClaimWork: PRE-FETCHED %zu sequences for shard %zu in %ldms (batch_size=%zu)",
		              static_cast<size_t>(seq_count), static_cast<size_t>(shard_idx), static_cast<long>(fetch_ms),
		              static_cast<size_t>(active->batch_size));
	} catch (...) {
		SHARD_DBG(gstate, "ClaimWork: ID/SEQUENCE FETCH FAILED shard %zu", static_cast<size_t>(shard_idx));
		active->exhausted.store(true, std::memory_order_release);
		active->active_workers.fetch_sub(1, std::memory_order_acq_rel);
		{
			std::lock_guard<std::mutex> guard(gstate.lock);
			auto &shards = gstate.active_shards;
			shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
		}
		gstate.cv.notify_all();
		throw;
	}

	idx_t estimated = bind_data.shards[shard_idx].read_count;
	if (seq_count > estimated) {
		auto old_total = gstate.total_associations.fetch_add(seq_count - estimated, std::memory_order_relaxed);
		SHARD_DBG(gstate, "ClaimWork: adjusted total_associations %zu -> %zu (shard %zu: estimated=%zu, actual=%zu)",
		          static_cast<size_t>(old_total), static_cast<size_t>(old_total + seq_count - estimated),
		          static_cast<size_t>(shard_idx), static_cast<size_t>(estimated), static_cast<size_t>(seq_count));
	} else if (seq_count < estimated) {
		auto old_total = gstate.total_associations.fetch_sub(estimated - seq_count, std::memory_order_relaxed);
		SHARD_DBG(gstate, "ClaimWork: adjusted total_associations %zu -> %zu (shard %zu: estimated=%zu, actual=%zu)",
		          static_cast<size_t>(old_total), static_cast<size_t>(old_total - (estimated - seq_count)),
		          static_cast<size_t>(shard_idx), static_cast<size_t>(estimated), static_cast<size_t>(seq_count));
	}

	// Phase 5: Publish under lock to prevent lost wake-ups with CV
	active->total_reads = seq_count;
	active->start_time = std::chrono::steady_clock::now();
	{
		std::lock_guard<std::mutex> guard(gstate.lock);
		active->ready.store(true, std::memory_order_release);
	}
	gstate.cv.notify_all();
	if (gstate.progress) {
		shard_progress::Emit("minibwa",
		                     shard_progress::FormatShardStart(static_cast<uint64_t>(shard_idx) + 1,
		                                                      static_cast<uint64_t>(gstate.shard_count),
		                                                      shard_info.name, static_cast<int64_t>(seq_count)));
	}
	return active;
}

void AlignMiniBWAShardedTableFunction::ReleaseWork(GlobalState &gstate, LocalState &lstate) {
	auto active = lstate.current_active_shard;
	lstate.aligner->detach_shared_index();
	lstate.has_shard = false;

	auto prev_workers = active->active_workers.fetch_sub(1, std::memory_order_acq_rel);

	SHARD_DBG(gstate, "ReleaseWork: shard %zu (workers %zu->%zu, exhausted=%s)", static_cast<size_t>(active->shard_idx),
	          static_cast<size_t>(prev_workers), static_cast<size_t>(prev_workers - 1),
	          active->exhausted.load(std::memory_order_relaxed) ? "yes" : "no");

	{
		std::lock_guard<std::mutex> guard(gstate.lock);
		if (prev_workers == 1 && active->exhausted.load(std::memory_order_acquire)) {
			auto &shards = gstate.active_shards;
			shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
			SHARD_DBG_MEM(gstate, "ReleaseWork: REMOVED shard %zu (active_shards=%zu)",
			              static_cast<size_t>(active->shard_idx), static_cast<size_t>(shards.size()));
			if (gstate.progress) {
				const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(
				                             std::chrono::steady_clock::now() - active->start_time)
				                             .count();
				shard_progress::Emit("minibwa",
				                     shard_progress::FormatShardDone(
				                         static_cast<uint64_t>(active->shard_idx) + 1,
				                         static_cast<uint64_t>(gstate.shard_count), lstate.current_shard_name,
				                         active->total_reads,
				                         active->alignments_emitted.load(std::memory_order_relaxed), elapsed_s));
			}
		}
		lstate.current_active_shard = nullptr;
		gstate.cv.notify_all();
	}
}

void AlignMiniBWAShardedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	while (true) {
		idx_t available = local_state.result_buffer.size() - local_state.buffer_offset;

		if (available > 0) {
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, local_state.result_buffer, local_state.buffer_offset, output_count,
			                     bind_data.query_schema.id_type, bind_data.subject_id_type);
			if (bind_data.include_shard_name) {
				auto shard_col_idx = output.ColumnCount() - 1;
				auto &shard_vec = output.data[shard_col_idx];
				for (idx_t i = 0; i < output_count; i++) {
					FlatVector::GetData<string_t>(shard_vec)[i] =
					    StringVector::AddString(shard_vec, local_state.current_shard_name);
				}
			}
			local_state.buffer_offset += output_count;
			return;
		}

		if (!local_state.has_shard) {
			auto active = ClaimWork(context, global_state, bind_data, local_state);
			if (!active) {
				output.SetCardinality(0);
				return;
			}
			local_state.current_active_shard = active;
			local_state.has_shard = true;
			local_state.current_shard_name = bind_data.shards[active->shard_idx].name;
			local_state.aligner->attach_shared_index(active->index);
		}

		auto &active = local_state.current_active_shard;
		idx_t seq_count = active->shard_sequences.size();
		idx_t my_offset = active->next_batch_offset.fetch_add(active->batch_size, std::memory_order_acq_rel);

		if (my_offset >= seq_count) {
			SHARD_DBG(global_state, "Execute: shard %zu offset %zu >= seq_count %zu, exhausted",
			          static_cast<size_t>(active->shard_idx), static_cast<size_t>(my_offset),
			          static_cast<size_t>(seq_count));
			active->exhausted.store(true, std::memory_order_release);
			ReleaseWork(global_state, local_state);
			continue;
		}

		idx_t batch_count = std::min(active->batch_size, seq_count - my_offset);
		bool is_last_batch = (my_offset + batch_count >= seq_count);

		global_state.associations_processed.fetch_add(batch_count, std::memory_order_relaxed);

		auto query_batch = active->shard_sequences.SubRange(my_offset, batch_count);
		SHARD_DBG(global_state, "Execute: shard %zu SUB-RANGE offset=%zu count=%zu (is_last=%s)",
		          static_cast<size_t>(active->shard_idx), static_cast<size_t>(my_offset),
		          static_cast<size_t>(query_batch.size()), is_last_batch ? "yes" : "no");

		if (query_batch.empty()) {
			if (is_last_batch) {
				active->exhausted.store(true, std::memory_order_release);
			}
			ReleaseWork(global_state, local_state);
			continue;
		}

		if (is_last_batch) {
			active->exhausted.store(true, std::memory_order_release);
		}

		local_state.result_buffer.clear();
		local_state.result_buffer.shrink_to_fit();
		SHARD_DBG_MEM(global_state, "Execute: shard %zu buffer cleared+shrunk before align",
		              static_cast<size_t>(active->shard_idx));
		local_state.buffer_offset = 0;

		if (!query_batch.empty()) {
			auto align_start = std::chrono::steady_clock::now();
			local_state.aligner->align(query_batch, local_state.result_buffer);
			auto align_ms =
			    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
			        .count();
			FilterMappedOnly(local_state.result_buffer);
			if (global_state.progress) {
				active->alignments_emitted.fetch_add(local_state.result_buffer.size(), std::memory_order_relaxed);
			}
			SHARD_DBG_MEM(global_state, "Execute: shard %zu ALIGN %zu reads -> %zu results in %ldms",
			              static_cast<size_t>(active->shard_idx), static_cast<size_t>(query_batch.size()),
			              static_cast<size_t>(local_state.result_buffer.size()), static_cast<long>(align_ms));
		}
	}
}

double AlignMiniBWAShardedTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                                  const GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto total = gstate.total_associations.load(std::memory_order_relaxed);
	if (total == 0) {
		return 100.0;
	}
	auto processed = gstate.associations_processed.load(std::memory_order_relaxed);
	return std::min(100.0, 100.0 * static_cast<double>(processed) / static_cast<double>(total));
}

TableFunction AlignMiniBWAShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_minibwa_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["max_threads_per_shard"] = LogicalType::INTEGER;
	tf.named_parameters["debug"] = LogicalType::BOOLEAN;
	tf.named_parameters["progress"] = LogicalType::BOOLEAN;
	tf.named_parameters["include_shard_name"] = LogicalType::BOOLEAN;

	tf.table_scan_progress = Progress;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignMiniBWAShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
