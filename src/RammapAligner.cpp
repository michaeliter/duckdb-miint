#include <RammapAligner.hpp>
#include <cstring>
#include <stdexcept>

namespace miint {

namespace {

// RAII wrapper for the char* rammap_align_query returns (Rust-allocated,
// released via rammap_string_free — not plain free()).
struct RammapStringPtr {
	char *ptr = nullptr;
	explicit RammapStringPtr(char *p) : ptr(p) {
	}
	~RammapStringPtr() {
		rammap_string_free(ptr);
	}
	RammapStringPtr(const RammapStringPtr &) = delete;
	RammapStringPtr &operator=(const RammapStringPtr &) = delete;
};

struct Bam1Deleter {
	void operator()(bam1_t *b) const {
		bam_destroy1(b);
	}
};

} // namespace

SharedRammapIndex::SharedRammapIndex(const std::vector<AlignmentSubject> &subjects, const RammapConfig &config) {
	if (subjects.empty()) {
		throw std::runtime_error("SharedRammapIndex: subjects cannot be empty");
	}

	std::vector<const char *> names;
	std::vector<const uint8_t *> seqs;
	std::vector<size_t> lens;
	names.reserve(subjects.size());
	seqs.reserve(subjects.size());
	lens.reserve(subjects.size());
	for (const auto &s : subjects) {
		names.push_back(s.read_id.c_str());
		seqs.push_back(reinterpret_cast<const uint8_t *>(s.sequence.data()));
		lens.push_back(s.sequence.size());
	}

	RammapIndexHandle *raw =
	    rammap_index_build(names.data(), seqs.data(), lens.data(), subjects.size(), config.preset.c_str(), config.k,
	                       config.w, config.eqx ? 1 : 0, config.max_secondary);
	if (!raw) {
		throw std::runtime_error("rammap_index_build failed (unrecognized preset '" + config.preset + "'?)");
	}
	index_.reset(raw);

	// Synthetic @SQ-only header so sam_parse1 can resolve RNAME to a tid
	// when parsing rammap's SAM text output back into bam1_t records.
	std::string header_text;
	for (const auto &s : subjects) {
		header_text += "@SQ\tSN:" + s.read_id + "\tLN:" + std::to_string(s.length()) + "\n";
	}
	sam_hdr_t *hdr = sam_hdr_parse(header_text.length(), header_text.c_str());
	if (!hdr) {
		throw std::runtime_error("SharedRammapIndex: failed to parse synthetic SAM header");
	}
	hdr_.reset(hdr);
}

RammapAligner::RammapAligner(std::shared_ptr<SharedRammapIndex> shared_index)
    : shared_index_(std::move(shared_index)), ctx_(rammap_ctx_new()) {
	if (!shared_index_) {
		throw std::runtime_error("RammapAligner: shared_index cannot be null");
	}
	if (!ctx_) {
		throw std::runtime_error("RammapAligner: rammap_ctx_new failed");
	}
}

void RammapAligner::align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output) {
	char *sam_text = rammap_align_query(shared_index_->raw(), ctx_.get(), read_id.c_str(),
	                                    reinterpret_cast<const uint8_t *>(sequence.data()), sequence.size());
	if (!sam_text) {
		throw std::runtime_error("rammap_align_query failed for read '" + read_id + "'");
	}
	RammapStringPtr owned_text(sam_text);

	// sam_text is 0+ newline-terminated SAM lines (no header). Parse each
	// with htslib's own text-SAM parser — the same one MIINT uses for
	// on-disk SAM files (see SAMReader.cpp) — so tag/position/CIGAR
	// handling is identical to that path; nothing is reconstructed here.
	const char *line_start = sam_text;
	while (*line_start != '\0') {
		const char *line_end = std::strchr(line_start, '\n');
		size_t line_len = line_end ? static_cast<size_t>(line_end - line_start) : std::strlen(line_start);
		if (line_len > 0) {
			kstring_t ks = KS_INITIALIZE;
			kputsn(line_start, line_len, &ks);

			std::unique_ptr<bam1_t, Bam1Deleter> aln(bam_init1());
			if (sam_parse1(&ks, shared_index_->header(), aln.get()) < 0) {
				ks_free(&ks);
				throw std::runtime_error("failed to parse rammap SAM output for read '" + read_id + "'");
			}
			sam_utils::parse_record_to_batch(aln.get(), shared_index_->header(), output, /*include_seq_qual=*/false);
			ks_free(&ks);

			// Bridge two conventions that don't otherwise meet: parse_record_to_batch
			// (SAMRecord.hpp) tracks absent int tags via the batch's tag_*_valid
			// vectors (0 = "no value pushed"), but align_common.hpp's
			// OutputSAMRecordBatch — reused here for consistency with
			// align_minimap2's output schema — ignores those _valid vectors
			// entirely and instead treats -1 as the null sentinel (see
			// SetAlignResultInt64Nullable). Without this, an absent tag (e.g.
			// rammap never emits XS/YS/XN/XM/XO/XG) would render as 0 instead
			// of NULL. Normalize just the row parse_record_to_batch appended.
			size_t last = output.size() - 1;
			if (!output.tag_as_valid[last]) {
				output.tag_as_values[last] = -1;
			}
			if (!output.tag_xs_valid[last]) {
				output.tag_xs_values[last] = -1;
			}
			if (!output.tag_ys_valid[last]) {
				output.tag_ys_values[last] = -1;
			}
			if (!output.tag_xn_valid[last]) {
				output.tag_xn_values[last] = -1;
			}
			if (!output.tag_xm_valid[last]) {
				output.tag_xm_values[last] = -1;
			}
			if (!output.tag_xo_valid[last]) {
				output.tag_xo_values[last] = -1;
			}
			if (!output.tag_xg_valid[last]) {
				output.tag_xg_values[last] = -1;
			}
			if (!output.tag_nm_valid[last]) {
				output.tag_nm_values[last] = -1;
			}
		}
		if (!line_end) {
			break;
		}
		line_start = line_end + 1;
	}
}

void RammapAligner::align(const SequenceRecordBatch &queries, SAMRecordBatch &output) {
	if (queries.empty()) {
		return;
	}
	output.reserve(output.size() + queries.size());
	for (size_t i = 0; i < queries.size(); i++) {
		bool actually_paired = queries.is_paired && i < queries.sequences2.size() && !queries.sequences2[i].empty();
		if (actually_paired) {
			throw std::runtime_error("align_rammap: paired-end queries are not yet supported (read_id='" +
			                         queries.read_ids[i] + "')");
		}
		if (queries.sequences1[i].empty()) {
			continue; // matches Minimap2Aligner::align_single's empty-sequence skip
		}
		align_single(queries.read_ids[i], queries.sequences1[i], output);
	}
}

} // namespace miint
