#pragma once
/*
 * align_rammap: v1 scope is subject_table + shared-index multi-thread mode
 * only (no index_path / per_subject_database — those exist for
 * align_minimap2 but rammap_ffi.rs doesn't yet expose an index-file loader
 * or single-subject rebuild path). See docs/internals/embedded-tools.md /
 * the align_rammap parity plan for the follow-up scope.
 */
#include "RammapAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <chrono>
#include <memory>
#include <vector>

namespace duckdb {

// Parse rammap config parameters from named_parameters (preset, max_secondary, k, w, eqx).
// Mirrors ParseMinimap2ConfigParams (align_common.hpp) but writes into a
// miint::RammapConfig, whose numeric fields default to 0 ("use preset
// default") rather than a hardcoded value — see RammapAligner.hpp.
inline void ParseRammapConfigParams(const named_parameter_map_t &params, miint::RammapConfig &config) {
	auto preset_param = params.find("preset");
	if (preset_param != params.end() && !preset_param->second.IsNull()) {
		config.preset = preset_param->second.ToString();
	}
	auto max_secondary_param = params.find("max_secondary");
	if (max_secondary_param != params.end() && !max_secondary_param->second.IsNull()) {
		config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
	}
	auto k_param = params.find("k");
	if (k_param != params.end() && !k_param->second.IsNull()) {
		config.k = k_param->second.GetValue<int32_t>();
	}
	auto w_param = params.find("w");
	if (w_param != params.end() && !w_param->second.IsNull()) {
		config.w = w_param->second.GetValue<int32_t>();
	}
	auto eqx_param = params.find("eqx");
	if (eqx_param != params.end() && !eqx_param->second.IsNull()) {
		config.eqx = eqx_param->second.GetValue<bool>();
	}
}

class AlignRammapTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string subject_table;
		miint::RammapConfig config;
		SequenceTableSchema query_schema;
		std::vector<miint::AlignmentSubject> subjects; // Pre-loaded at bind time
		LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data() : names(GetAlignmentOutputNames()), types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	struct StandardModeState {
		std::shared_ptr<miint::SharedRammapIndex> shared_index;
		std::unique_ptr<QuerySequenceStream> query_stream;
	};

	struct GlobalState : public GlobalTableFunctionState {
		idx_t num_threads = 1;
		std::atomic<idx_t> init_local_count {0};
		std::unique_ptr<StandardModeState> standard;

		idx_t MaxThreads() const override {
			return num_threads;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::RammapAligner> aligner;
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
