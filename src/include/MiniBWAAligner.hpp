#pragma once

#include "Minimap2Aligner.hpp" // reuse miint::AlignmentSubject, miint::AlignmentStats
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <minibwa/minibwa.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Configuration for minibwa alignment. Presets mirror the CLI's own:
// "sr"/"adap" (short/adaptive read, MB_F_PE on by default upstream -- cleared
// explicitly here since this phase is single-end only) and "lr" (long read).
// max_secondary maps to mb_opt_t::best_n; minibwa's own sr/adap default (50)
// is intentionally kept rather than matching minimap2's lower default (5) --
// see options.c's mb_opt_preset().
struct MiniBWAConfig {
	std::string preset = "sr";
	int max_secondary = 50;
	bool is_meth = false; // BS-seq mode; plumbed for a future phase, not wired into align_minibwa yet
	int sa_bit = 4;       // SA sample-rate exponent for index building (CLI default)
};

// Custom deleters
struct MiniBWAIndexDeleter {
	void operator()(mb_idx_t *idx) const;
};
struct MiniBWATbufDeleter {
	void operator()(mb_tbuf_t *tbuf) const;
};

using MiniBWAIndexPtr = std::unique_ptr<mb_idx_t, MiniBWAIndexDeleter>;
using MiniBWATbufPtr = std::unique_ptr<mb_tbuf_t, MiniBWATbufDeleter>;

// Shared, immutable minibwa index for multi-thread-per-shard alignment.
// mb_map()/mb_map_batch() take a `const mb_idx_t *`, so the index is safe to
// read concurrently from many threads; each thread needs its own mb_tbuf_t
// (see MiniBWAAligner below). Mirrors SharedMinimap2Index.
class SharedMiniBWAIndex {
public:
	SharedMiniBWAIndex(const std::string &index_prefix, const MiniBWAConfig &config);
	~SharedMiniBWAIndex();

	SharedMiniBWAIndex(const SharedMiniBWAIndex &) = delete;
	SharedMiniBWAIndex &operator=(const SharedMiniBWAIndex &) = delete;

	const mb_idx_t *index() const {
		return index_.get();
	}
	const mb_opt_t &opt() const {
		return opt_;
	}
	const std::vector<std::string> &subject_names() const {
		return subject_names_;
	}

private:
	MiniBWAIndexPtr index_;
	mb_opt_t opt_;
	std::vector<std::string> subject_names_;
};

// Main aligner class.
// NOT thread-safe: each thread must have its own MiniBWAAligner instance
// (specifically its own mb_tbuf_t arena). Multiple instances may share a
// SharedMiniBWAIndex concurrently.
class MiniBWAAligner {
public:
	explicit MiniBWAAligner(const MiniBWAConfig &config);
	~MiniBWAAligner() = default;

	MiniBWAAligner(const MiniBWAAligner &) = delete;
	MiniBWAAligner &operator=(const MiniBWAAligner &) = delete;
	MiniBWAAligner(MiniBWAAligner &&) noexcept = default;
	MiniBWAAligner &operator=(MiniBWAAligner &&) noexcept = default;

	// Load a prebuilt index (<prefix>.l2b / <prefix>.mbw, written by
	// mb_idx_build() -- see save_minibwa_index.cpp).
	void load_index(const std::string &index_prefix);

	// Attach a shared index (clears any owned index)
	void attach_shared_index(std::shared_ptr<SharedMiniBWAIndex> shared_idx);
	// Detach the shared index (does not destroy it; other aligners may still reference it)
	void detach_shared_index();

	static std::shared_ptr<SharedMiniBWAIndex> BuildSharedIndex(const std::string &index_prefix,
	                                                            const MiniBWAConfig &config);

	// Align single-end queries against the current index, appending to output.
	// Paired-end support (mb_map_batch_pe) is a later phase; queries.sequences2
	// is ignored here even if present.
	void align(const SequenceRecordBatch &queries, SAMRecordBatch &output);

	static void InitOptions(const MiniBWAConfig &config, mb_opt_t &opt);
	static std::vector<std::string> EnumerateContigNames(const mb_idx_t *idx);

private:
	MiniBWAConfig config_;
	mb_opt_t opt_;
	MiniBWAIndexPtr index_;              // owned index (mutually exclusive with shared_index_)
	std::vector<std::string> subject_names_;
	MiniBWATbufPtr tbuf_;
	std::shared_ptr<SharedMiniBWAIndex> shared_index_;

	const mb_idx_t *active_index() const;
	const mb_opt_t &active_opt() const;
	const std::vector<std::string> &active_subject_names() const;

	void align_batch(const std::vector<std::string> &read_ids, const std::vector<std::string> &sequences,
	                 SAMRecordBatch &output);

	void hit_to_sam(const mb_hit_t &hit, const std::string &read_id, int32_t query_len, SAMRecordBatch &batch) const;
	std::string cigar_string(const mb_hit_t &hit, int32_t query_len, uint16_t sam_flags,
	                         AlignmentStats *stats_out) const;
	uint16_t calculate_flags(const mb_hit_t &hit, bool is_unmapped) const;
	const std::string &get_reference_name(int64_t tid) const;
};

} // namespace miint
