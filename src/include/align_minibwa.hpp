#pragma once
#include "MiniBWAAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <memory>
#include <vector>

namespace duckdb {

// index_path only in this phase (pre-built index, loaded via mb_idx_load).
// subject_table (build-on-the-fly) is deferred to pair with
// save_minibwa_index, since both need the same "stage to temp FASTA, call
// mb_idx_build()" logic -- see the integration plan's Phase 3.
class AlignMiniBWATableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string index_path;
		miint::MiniBWAConfig config;
		SequenceTableSchema query_schema;

		// Subject names in minibwa's index files are opaque bytes -- always VARCHAR.
		LogicalType subject_id_type = LogicalType::VARCHAR;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data() : names(GetAlignmentOutputNames()), types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		idx_t num_threads = 1;
		std::shared_ptr<miint::SharedMiniBWAIndex> shared_index;
		std::unique_ptr<QuerySequenceStream> query_stream;

		idx_t MaxThreads() const override {
			return num_threads;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::MiniBWAAligner> aligner;
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
