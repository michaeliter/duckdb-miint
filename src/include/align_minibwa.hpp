#pragma once
#include "MiniBWAAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "save_minibwa_index.hpp" // MakeMiniBWATempDir / WriteMiniBWASubjectsFasta, shared for subject_table mode
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

// Shared with align_minibwa_sharded.cpp. Unlike ParseMinimap2ConfigParams
// (align_common.hpp), there's no warn_prebuilt_index case -- MiniBWAConfig
// has no analogue of minimap2's k/w baked-into-the-index parameters.
void ParseMiniBWAConfigParams(const named_parameter_map_t &params, miint::MiniBWAConfig &config);

// Exactly one of subject_table (build on-the-fly) or index_path (prebuilt)
// must be given, mirroring align_minimap2's contract. Unlike minimap2's
// mm_idx_str(), minibwa's index builder only accepts a FASTA file path, so
// subject_table mode still round-trips through a temp FASTA (same staging
// helpers save_minibwa_index.cpp uses) -- which means, unlike minimap2,
// subject_table mode's output `reference` column is VARCHAR always, never
// mirroring the subject table's original id type: the on-disk index's
// contig names are the FASTA-header string form regardless of how the
// index was built, so there is no path back to a BIGINT/UUID reference
// column the way minimap2's in-memory mm_idx_str() index preserves.
//
// per_subject_database (rebuild the index per subject, align every query
// against each) is not implemented -- for minibwa this would mean a
// disk round trip per subject (temp FASTA + suffix-array build) rather
// than minimap2's cheap in-memory mm_idx_str() call per subject.
class AlignMiniBWATableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string subject_table; // OPTIONAL (either this or index_path required)
		std::string index_path;    // OPTIONAL: path prefix to a prebuilt index
		miint::MiniBWAConfig config;
		SequenceTableSchema query_schema;
		std::vector<miint::AlignmentSubject> subjects; // Pre-loaded at bind time (empty if using index_path)

		bool using_prebuilt_index() const {
			return !index_path.empty();
		}

		// Subject names in minibwa's index files are opaque bytes -- always
		// VARCHAR, for both index_path and subject_table (see the class
		// comment above for why subject_table can't preserve the original type).
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
