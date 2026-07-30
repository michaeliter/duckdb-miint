#pragma once
#include "MiniBWAAligner.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

// save_minibwa_index(subject_table, output_path, [sa_bit], [threads]) builds
// a minibwa index from a single-end subject table and persists it at
// `output_path` (a prefix: writes <output_path>.l2b and <output_path>.mbw).
// The minibwa analogue of save_minimap2_index -- except mb_idx_build() (the
// ext/minibwa fork patch backing this) only accepts a FASTA file path, not
// in-memory sequences, so subjects are staged to a temp FASTA first (same
// mkdtemp/RAII-cleanup shape as align_bowtie2's temp index dir).
//
// Known limitation: minibwa's index round-trips subjects through an on-disk
// FASTA (kseq-parsed), unlike minimap2's mm_idx_str() which builds directly
// from in-memory C strings. A read_id containing whitespace will be
// truncated at the first space when re-read back as a contig name (kseq
// convention: everything after the first space becomes the comment, not the
// name) -- consistently, so align_minibwa's `reference` output stays
// internally consistent with the index, but will differ from the original
// subject_table read_id in that edge case.
class SaveMiniBWAIndexTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string subject_table;
		std::string output_path;
		miint::MiniBWAConfig config;
		std::vector<miint::AlignmentSubject> subjects;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names({"success", "index_path", "num_subjects"}),
		      types({LogicalType::BOOLEAN, LogicalType::VARCHAR, LogicalType::BIGINT}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		bool done = false;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	struct LocalState : public LocalTableFunctionState {};

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
