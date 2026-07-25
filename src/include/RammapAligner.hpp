#pragma once

#include "Minimap2Aligner.hpp" // reuse miint::AlignmentSubject
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "rammap_ffi.h"
#include <htslib-1.22.1/htslib/sam.h>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Configuration for rammap alignment. Unlike Minimap2Config, the numeric
// overrides default to 0 ("use the preset's own default") rather than a
// hardcoded value — matches the CLI's "only override if the flag was
// passed" semantics (see rammap_ffi.rs), which is what byte-parity with the
// standalone `rammap` CLI depends on.
struct RammapConfig {
	std::string preset = "sr";
	int32_t k = 0;             // 0 = preset default
	int32_t w = 0;             // 0 = preset default
	bool eqx = true;           // =/X CIGAR ops
	int32_t max_secondary = 0; // 0 = preset default; maps to rammap's filtering.best_n
};

struct RammapIndexDeleter {
	void operator()(RammapIndexHandle *h) const {
		rammap_index_free(h);
	}
};
struct RammapCtxDeleter {
	void operator()(RammapCtxHandle *h) const {
		rammap_ctx_free(h);
	}
};
using RammapIndexPtr = std::unique_ptr<RammapIndexHandle, RammapIndexDeleter>;
using RammapCtxPtr = std::unique_ptr<RammapCtxHandle, RammapCtxDeleter>;

// Shared, immutable rammap index for multi-thread-per-shard alignment.
// Multiple RammapAligner instances may reference the same SharedRammapIndex
// concurrently (each aligner owns its own RammapCtxHandle). Mirrors
// SharedMinimap2Index.
class SharedRammapIndex {
public:
	SharedRammapIndex(const std::vector<AlignmentSubject> &subjects, const RammapConfig &config);
	~SharedRammapIndex() = default;

	SharedRammapIndex(const SharedRammapIndex &) = delete;
	SharedRammapIndex &operator=(const SharedRammapIndex &) = delete;

	RammapIndexHandle *raw() const {
		return index_.get();
	}
	// Synthetic @SQ-only header used to parse rammap's SAM text output back
	// into bam1_t records (rammap_align_query returns data lines only, no
	// header — see rammap_ffi.rs).
	sam_hdr_t *header() const {
		return hdr_.get();
	}

private:
	struct SamHdrDeleter {
		void operator()(sam_hdr_t *h) const {
			if (h) {
				sam_hdr_destroy(h);
			}
		}
	};
	RammapIndexPtr index_;
	std::unique_ptr<sam_hdr_t, SamHdrDeleter> hdr_;
};

// Main aligner class. NOT thread-safe: each thread must have its own
// RammapAligner instance. Multiple instances may share a SharedRammapIndex
// concurrently. Mirrors Minimap2Aligner, single-end only for v1 (parity is
// proven read-by-read before paired-end is added).
class RammapAligner {
public:
	explicit RammapAligner(std::shared_ptr<SharedRammapIndex> shared_index);
	~RammapAligner() = default;

	RammapAligner(const RammapAligner &) = delete;
	RammapAligner &operator=(const RammapAligner &) = delete;
	RammapAligner(RammapAligner &&) noexcept = default;
	RammapAligner &operator=(RammapAligner &&) noexcept = default;

	// Align a single-end read against the shared index, appending result(s)
	// to `output`. Parses rammap's own SAM text via htslib's sam_parse1 (the
	// same parser MIINT uses for on-disk SAM files) so tag/position/CIGAR
	// handling is identical to the read_alignments() path — no field is
	// reconstructed in C++.
	void align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output);

	// Align a batch of queries. v1 is single-end only: a paired-end query
	// (queries.is_paired with a non-empty sequence2 for some row) throws —
	// paired-end rammap alignment is a follow-up (see rammap_ffi.rs notes on
	// align_and_format_pair not yet being exposed).
	void align(const SequenceRecordBatch &queries, SAMRecordBatch &output);

private:
	std::shared_ptr<SharedRammapIndex> shared_index_;
	RammapCtxPtr ctx_;
};

} // namespace miint
