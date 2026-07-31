#pragma once
#include "MiniBWAAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "align_minibwa.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

// Reference sharding with pre-routed reads, mirroring
// AlignMinimap2ShardedTableFunction's architecture exactly (shared
// in-process index, cooperative multi-thread-per-shard scheduler) rather
// than bowtie2-sharded's one-shard-per-worker-daemon model -- minibwa, like
// minimap2, is a statically linked library whose mb_map*() calls take a
// `const mb_idx_t *` and are safe to call concurrently against one shared
// index (see MiniBWAAligner.hpp). Shard index files are
// <shard_directory>/<shard_name>.l2b + .mbw (minibwa's two-file index,
// vs minimap2's single .mmi).

struct MiniBWAShardInfo {
	std::string name;
	std::string index_prefix; // e.g. "/path/shards/shard_001" (no extension)
	idx_t read_count;
};

// A shard currently being processed by one or more threads. Mirrors
// ActiveShard (align_minimap2_sharded.hpp) field-for-field.
struct MiniBWAActiveShard {
	idx_t shard_idx;
	idx_t batch_size;
	miint::SequenceRecordBatch shard_sequences;
	std::shared_ptr<miint::SharedMiniBWAIndex> index;
	std::atomic<idx_t> next_batch_offset {0};
	std::atomic<idx_t> active_workers {0};
	std::atomic<bool> exhausted {false};
	std::atomic<bool> ready {false};
	std::atomic<idx_t> alignments_emitted {0};
	idx_t total_reads = 0;
	std::chrono::steady_clock::time_point start_time;
};

class AlignMiniBWAShardedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string shard_directory;
		std::string read_to_shard_table;
		SequenceTableSchema query_schema;
		miint::MiniBWAConfig config;
		std::vector<MiniBWAShardInfo> shards; // Sorted by read_count DESC (largest first)
		idx_t max_threads_per_shard = 4;
		bool debug = false;
		bool progress = false;
		bool include_shard_name = false;

		// Subject-side id type: sharded mode always loads prebuilt indexes
		// whose subject names are opaque bytes, same contract as
		// align_minibwa(index_path:=...). VARCHAR once Bind runs.
		LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::condition_variable cv;
		idx_t next_shard_idx = 0;
		idx_t shard_count = 0;
		idx_t max_threads_per_shard = 4;
		idx_t max_active_shards = 1;
		bool debug = false;
		bool progress = false;
		std::chrono::steady_clock::time_point start_time;
		std::vector<std::shared_ptr<MiniBWAActiveShard>> active_shards;
		std::atomic<idx_t> total_associations {0};
		std::atomic<idx_t> associations_processed {0};

		idx_t MaxThreads() const override {
			return max_active_shards * max_threads_per_shard;
		}

		GlobalState() = default;
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::MiniBWAAligner> aligner;
		std::shared_ptr<MiniBWAActiveShard> current_active_shard;
		bool has_shard = false;
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
		std::string current_shard_name;

		LocalState() = default;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	static double Progress(ClientContext &context, const FunctionData *bind_data,
	                       const GlobalTableFunctionState *global_state);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);

private:
	static std::shared_ptr<MiniBWAActiveShard> ClaimWork(ClientContext &context, GlobalState &gstate,
	                                                     const Data &bind_data, LocalState &lstate);
	static void ReleaseWork(GlobalState &gstate, LocalState &lstate);
};

} // namespace duckdb
