#include "align_rammap.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb {

unique_ptr<FunctionData> AlignRammapTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                        vector<LogicalType> &return_types,
                                                        vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 1) {
		throw BinderException("align_rammap requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	auto subject_param = input.named_parameters.find("subject_table");
	if (subject_param != input.named_parameters.end() && !subject_param->second.IsNull()) {
		data->subject_table = subject_param->second.ToString();
	}
	if (data->subject_table.empty()) {
		throw BinderException("align_rammap requires subject_table parameter");
	}

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);

	ParseRammapConfigParams(input.named_parameters, data->config);

	auto subject_schema = ValidateSequenceTableSchema(context, data->subject_table, /*allow_bigint=*/true);
	data->subject_id_type = subject_schema.id_type;
	data->subjects = ReadSubjectTable(context, data->subject_table, subject_schema);

	data->types = GetAlignmentOutputTypes(data->query_schema.id_type, data->subject_id_type);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignRammapTableFunction::InitGlobal(ClientContext &context,
                                                                          TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	auto st = std::make_unique<StandardModeState>();
	try {
		st->shared_index = std::make_shared<miint::SharedRammapIndex>(data.subjects, data.config);
	} catch (const std::exception &e) {
		throw IOException("Failed to build rammap index: %s", e.what());
	}
	gstate->standard = std::move(st);
	gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());

	gstate->standard->query_stream =
	    std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema);

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignRammapTableFunction::InitLocal(ExecutionContext &context,
                                                                        TableFunctionInitInput &input,
                                                                        GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();
	lstate->aligner = std::make_unique<miint::RammapAligner>(gstate.standard->shared_index);
	gstate.init_local_count.fetch_add(1);
	return lstate;
}

void AlignRammapTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();
	auto &st = *gstate.standard;

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

		auto query_batch = st.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		lstate.aligner->align(query_batch, lstate.result_buffer);
	}
}

TableFunction AlignRammapTableFunction::GetFunction() {
	auto tf = TableFunction("align_rammap", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	tf.named_parameters["subject_table"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["w"] = LogicalType::INTEGER;
	tf.named_parameters["eqx"] = LogicalType::BOOLEAN;

	// Matches align_minimap2: alignment output order is non-deterministic
	// (thread scheduling), so NO_ORDER lets DuckDB parallelize CTAS
	// pipelines instead of serializing via preserve_insertion_order.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignRammapTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
