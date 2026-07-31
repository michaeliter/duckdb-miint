#include "align_minibwa.hpp"
#include "align_common.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include <filesystem>

namespace duckdb {

void ParseMiniBWAConfigParams(const named_parameter_map_t &params, miint::MiniBWAConfig &config) {
	auto preset_param = params.find("preset");
	if (preset_param != params.end() && !preset_param->second.IsNull()) {
		config.preset = preset_param->second.ToString();
	}

	auto max_secondary_param = params.find("max_secondary");
	if (max_secondary_param != params.end() && !max_secondary_param->second.IsNull()) {
		config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
		if (config.max_secondary < 0) {
			throw InvalidInputException("max_secondary must be >= 0");
		}
	}
}

unique_ptr<FunctionData> AlignMiniBWATableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 1) {
		throw BinderException("align_minibwa requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	auto subject_param = input.named_parameters.find("subject_table");
	if (subject_param != input.named_parameters.end() && !subject_param->second.IsNull()) {
		data->subject_table = subject_param->second.ToString();
	}

	auto index_path_param = input.named_parameters.find("index_path");
	if (index_path_param != input.named_parameters.end() && !index_path_param->second.IsNull()) {
		data->index_path = index_path_param->second.ToString();
	}

	bool has_subject = !data->subject_table.empty();
	bool has_index = data->using_prebuilt_index();
	if (!has_subject && !has_index) {
		throw BinderException("align_minibwa requires either subject_table or index_path parameter");
	}
	if (has_subject && has_index) {
		throw BinderException(
		    "align_minibwa: Cannot specify both subject_table and index_path. "
		    "Use subject_table to build an index from sequences, or index_path to load a pre-built index.");
	}

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);

	ParseMiniBWAConfigParams(input.named_parameters, data->config);

	auto sa_bit_param = input.named_parameters.find("sa_bit");
	if (sa_bit_param != input.named_parameters.end() && !sa_bit_param->second.IsNull()) {
		data->config.sa_bit = sa_bit_param->second.GetValue<int32_t>();
		if (data->config.sa_bit <= 0) {
			throw InvalidInputException("sa_bit must be > 0");
		}
		if (has_index) {
			Printer::Print("WARNING: Parameter 'sa_bit' is ignored when using a pre-built index. "
			               "The SA sample rate is baked into the index.\n");
		}
	}

	if (has_index) {
		// Advisory existence check (TOCTOU, like align_minimap2's is_index_file
		// check) -- minibwa has no single-file magic-number probe like
		// mm_idx_is_idx; the index is a pair of files.
		auto &fs = FileSystem::GetFileSystem(context);
		if (!fs.FileExists(data->index_path + ".l2b") || !fs.FileExists(data->index_path + ".mbw")) {
			throw BinderException("minibwa index not found at prefix: %s (expected %s.l2b and %s.mbw)",
			                      data->index_path, data->index_path, data->index_path);
		}
	} else {
		// subject_table mode: load subjects now (matches align_minimap2's
		// pattern -- subjects are read once at bind time). The actual index
		// build (staging a temp FASTA, calling mb_idx_build) happens in
		// InitGlobal, same split as save_minibwa_index's Bind/InitGlobal.
		auto subject_schema = ValidateSequenceTableSchema(context, data->subject_table, /*allow_bigint=*/true);
		data->subjects = ReadSubjectTable(context, data->subject_table, subject_schema);
		if (data->subjects.empty()) {
			throw BinderException("Subject table '%s' is empty. Cannot align against an empty subject set.",
			                      data->subject_table);
		}
		// subject_id_type stays VARCHAR regardless of subject_schema.id_type --
		// see the class comment on Data for why minibwa can't preserve it.
	}

	data->types = GetAlignmentOutputTypes(data->query_schema.id_type, data->subject_id_type);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMiniBWATableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	if (data.using_prebuilt_index()) {
		try {
			gstate->shared_index = miint::MiniBWAAligner::BuildSharedIndex(data.index_path, data.config);
		} catch (const std::exception &e) {
			throw IOException("Failed to load minibwa index from '%s': %s", data.index_path, e.what());
		}
	} else {
		// Build a combined index from all subjects, staged through one temp
		// FASTA -- mb_idx_build() only accepts a file path (see the class
		// comment on Data). Mirrors align_minimap2's default (non-
		// per_subject_database) subject_table path: one index, built once,
		// all queries aligned against it.
		std::string temp_dir = MakeMiniBWATempDir("align_minibwa");
		std::string fasta_path = temp_dir + "/subjects.fa";
		std::string index_prefix = temp_dir + "/idx";
		try {
			WriteMiniBWASubjectsFasta(fasta_path, data.subjects, data.subject_table, "align_minibwa");
			int rc = mb_idx_build(fasta_path.c_str(), index_prefix.c_str(), data.config.sa_bit, /*n_thread=*/1,
			                      /*is_meth=*/0, /*seed=*/11);
			if (rc != 0) {
				throw IOException("align_minibwa: mb_idx_build failed (rc=%d) for staged FASTA: %s", rc, fasta_path);
			}
			gstate->shared_index = miint::MiniBWAAligner::BuildSharedIndex(index_prefix, data.config);
		} catch (...) {
			std::error_code ec;
			std::filesystem::remove_all(temp_dir, ec); // best-effort; can't propagate ec from an exception path
			throw;
		}
		// mb_idx_load() reads the whole index into heap memory (not mmap), so
		// once BuildSharedIndex has returned, the on-disk files are no longer
		// needed -- safe to remove immediately, unlike e.g. a daemon-backed
		// index that needs its files to persist for the session.
		std::error_code ec;
		std::filesystem::remove_all(temp_dir, ec);
	}

	gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	gstate->query_stream = std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema);

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMiniBWATableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();

	lstate->aligner = std::make_unique<miint::MiniBWAAligner>(data.config);
	lstate->aligner->attach_shared_index(gstate.shared_index);

	return lstate;
}

void AlignMiniBWATableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();

	while (true) {
		idx_t available = lstate.result_buffer.size() - lstate.buffer_offset;
		if (available > 0) {
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, lstate.result_buffer, lstate.buffer_offset, output_count,
			                     bind_data.query_schema.id_type, bind_data.subject_id_type);
			lstate.buffer_offset += output_count;
			return;
		}

		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		lstate.aligner->align(query_batch, lstate.result_buffer);
	}
}

TableFunction AlignMiniBWATableFunction::GetFunction() {
	auto tf = TableFunction("align_minibwa", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	tf.named_parameters["subject_table"] = LogicalType::VARCHAR;
	tf.named_parameters["index_path"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["sa_bit"] = LogicalType::INTEGER;

	// Alignment output order is non-deterministic (thread scheduling), so
	// NO_ORDER lets DuckDB parallelize CTAS pipelines instead of serializing
	// via preserve_insertion_order. Matches align_minimap2.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignMiniBWATableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
