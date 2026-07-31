#include "MiniBWAAligner.hpp"
#include <algorithm>
#include <cstdlib>
#include <stdexcept>

// mb_hit_t / mb_extra_t use the same BAM CIGAR op encoding, op string
// ("MIDNSHP=XB"), and O(1) NM formula (blen - mlen + n_ambi) as minimap2's
// mm_reg1_t -- both descend from the same klib/ksw2 lineage. The CIGAR walk
// below mirrors Minimap2Aligner::cigar_string exactly for this reason.
//
// minibwa's arena allocator (mb_tbuf_init(0), used here to match the CLI's
// own default) aborts the whole process on genuine OOM (kalloc.c's
// morecore() -> panic() -> abort()). Switching to mb_tbuf_init(1) (raw
// malloc/free, no arena) does not meaningfully avoid this: minibwa's higher
// call sites do not null-check allocation results either way, so an actual
// OOM crashes one way or the other. DuckDB's own memory_limit accounting is
// what keeps normal queries from reaching real system OOM in the first
// place; this is a known, inherent limitation of the vendored library, not
// unique to minibwa among the embedded C aligners here.

namespace miint {

void MiniBWAIndexDeleter::operator()(mb_idx_t *idx) const {
	if (idx) {
		mb_idx_destroy(idx);
	}
}

void MiniBWATbufDeleter::operator()(mb_tbuf_t *tbuf) const {
	if (tbuf) {
		mb_tbuf_destroy(tbuf);
	}
}

void MiniBWAAligner::InitOptions(const MiniBWAConfig &config, mb_opt_t &opt) {
	mb_opt_init(&opt);
	if (mb_opt_preset(&opt, config.preset.c_str()) != 0) {
		throw std::runtime_error("Unknown minibwa preset: " + config.preset);
	}
	// This phase is single-end only; mb_opt_preset("sr"/"adap") sets MB_F_PE
	// by default (see ext/minibwa/options.c). Paired-end support is a later
	// phase (mb_map_batch_pe).
	opt.flag &= ~MB_F_PE;
	// best_n (internal retention) intentionally left at the preset's own
	// value -- see the MiniBWAConfig comment on why this must not be
	// max_secondary.
	opt.out_n = config.max_secondary;
	opt.out_s = config.secondary_score_ratio;
}

std::vector<std::string> MiniBWAAligner::EnumerateContigNames(const mb_idx_t *idx) {
	// minibwa's public API has no mb_idx_n_ctg(); mb_idx_ctg_name() returns
	// NULL past the last contig (see ext/minibwa/map-algo.c), so enumerate
	// until that sentinel instead of needing a further fork patch.
	std::vector<std::string> names;
	for (int32_t tid = 0;; tid++) {
		const char *name = mb_idx_ctg_name(idx, tid);
		if (!name) {
			break;
		}
		names.emplace_back(name);
	}
	return names;
}

// SharedMiniBWAIndex
SharedMiniBWAIndex::SharedMiniBWAIndex(const std::string &index_prefix, const MiniBWAConfig &config) : opt_() {
	MiniBWAAligner::InitOptions(config, opt_);
	mb_idx_t *raw_idx = mb_idx_load(index_prefix.c_str(), config.is_meth ? 1 : 0);
	if (!raw_idx) {
		throw std::runtime_error("Cannot load minibwa index: " + index_prefix);
	}
	index_.reset(raw_idx);
	subject_names_ = MiniBWAAligner::EnumerateContigNames(index_.get());
}

SharedMiniBWAIndex::~SharedMiniBWAIndex() = default;

// MiniBWAAligner
MiniBWAAligner::MiniBWAAligner(const MiniBWAConfig &config) : config_(config), opt_(), tbuf_(mb_tbuf_init(0)) {
	InitOptions(config_, opt_);
}

void MiniBWAAligner::attach_shared_index(std::shared_ptr<SharedMiniBWAIndex> shared_idx) {
	index_.reset();
	subject_names_.clear();
	shared_index_ = std::move(shared_idx);
}

void MiniBWAAligner::detach_shared_index() {
	shared_index_.reset();
}

std::shared_ptr<SharedMiniBWAIndex> MiniBWAAligner::BuildSharedIndex(const std::string &index_prefix,
                                                                     const MiniBWAConfig &config) {
	return std::make_shared<SharedMiniBWAIndex>(index_prefix, config);
}

void MiniBWAAligner::load_index(const std::string &index_prefix) {
	mb_idx_t *raw_idx = mb_idx_load(index_prefix.c_str(), config_.is_meth ? 1 : 0);
	if (!raw_idx) {
		throw std::runtime_error("Cannot load minibwa index: " + index_prefix);
	}
	index_.reset(raw_idx);
	subject_names_ = EnumerateContigNames(index_.get());
	shared_index_.reset();
}

const mb_idx_t *MiniBWAAligner::active_index() const {
	if (shared_index_) {
		return shared_index_->index();
	}
	return index_.get();
}

const mb_opt_t &MiniBWAAligner::active_opt() const {
	if (shared_index_) {
		return shared_index_->opt();
	}
	return opt_;
}

const std::vector<std::string> &MiniBWAAligner::active_subject_names() const {
	if (shared_index_) {
		return shared_index_->subject_names();
	}
	return subject_names_;
}

void MiniBWAAligner::align(const SequenceRecordBatch &queries, SAMRecordBatch &output) {
	if (queries.empty()) {
		return;
	}
	if (!active_index()) {
		throw std::runtime_error("No index loaded. Call load_index() or attach_shared_index() first.");
	}

	output.reserve(output.size() + queries.size());

	if (queries.is_paired) {
		// The query table's schema having a sequence2 column (e.g.
		// read_fastx()'s general schema always includes one) doesn't mean
		// this batch is actually paired data -- single-end input read via
		// read_fastx() with no mate file leaves sequence2 NULL/empty for
		// every row. Fall back to single-end for an all-empty batch, same
		// as Minimap2Aligner's per-row graceful degradation; a genuine mix
		// of paired and unpaired rows in one batch still throws in
		// align_paired_batch, since minibwa's flat-interleaved-pairs
		// contract has no per-row opt-out the way minimap2's does.
		bool any_paired = false;
		for (const auto &s2 : queries.sequences2) {
			if (!s2.empty()) {
				any_paired = true;
				break;
			}
		}
		if (any_paired) {
			align_paired_batch(queries, output);
		} else {
			align_batch(queries.read_ids, queries.sequences1, output);
		}
	} else {
		align_batch(queries.read_ids, queries.sequences1, output);
	}
}

void MiniBWAAligner::align_batch(const std::vector<std::string> &read_ids, const std::vector<std::string> &sequences,
                                 SAMRecordBatch &output) {
	std::vector<int32_t> qlens;
	std::vector<const char *> seq_ptrs;
	std::vector<const char *> name_ptrs;
	qlens.reserve(sequences.size());
	seq_ptrs.reserve(sequences.size());
	name_ptrs.reserve(sequences.size());

	for (size_t i = 0; i < sequences.size(); i++) {
		qlens.push_back(static_cast<int32_t>(sequences[i].length()));
		seq_ptrs.push_back(sequences[i].c_str());
		name_ptrs.push_back(read_ids[i].c_str());
	}

	std::vector<int32_t> n_hit(sequences.size(), 0);
	mb_hit_t **hits = mb_map_batch(&active_opt(), active_index(), static_cast<int32_t>(sequences.size()), qlens.data(),
	                               seq_ptrs.data(), n_hit.data(), tbuf_.get(), name_ptrs.data());
	if (!hits) {
		// mb_map_batch returns NULL for n_seq <= 0 or a meth-index mismatch;
		// neither applies here (checked non-empty above, is_meth left off).
		return;
	}

	const mb_opt_t &opt = active_opt();
	for (size_t i = 0; i < sequences.size(); i++) {
		int32_t n_sec = 0; // count of secondaries emitted so far for this read
		for (int32_t j = 0; j < n_hit[i]; j++) {
			const mb_hit_t &hit = hits[i][j];
			bool is_primary = (hit.parent == hit.id);
			// Exact gate from map-main.c: primary always passes; a secondary
			// only if there's still room in the output budget (opt.out_n).
			if (!(is_primary || n_sec < opt.out_n)) {
				continue;
			}
			if (!is_primary) {
				const mb_hit_t &parent = hits[i][hit.parent];
				if (parent.p && hit.p) {
					if (hit.p->dp_max < opt.out_s * parent.p->dp_max) {
						continue;
					}
				} else if (hit.score < opt.out_s * parent.score) {
					continue;
				}
			}
			hit_to_sam(hit, read_ids[i], qlens[i], output, -1, false, false, -1, 0, 0);
			if (!is_primary) {
				n_sec++;
			}
		}
		for (int32_t j = 0; j < n_hit[i]; j++) {
			free(hits[i][j].p);
		}
		free(hits[i]);
	}
	free(hits);
}

void MiniBWAAligner::align_paired_batch(const SequenceRecordBatch &queries, SAMRecordBatch &output) {
	size_t n_seq = queries.size();
	for (size_t i = 0; i < n_seq; i++) {
		if (queries.sequences2[i].empty()) {
			throw std::runtime_error("MiniBWAAligner: paired-end batch has an empty sequence2 at row " +
			                         std::to_string(i) + " (read_id=" + queries.read_ids[i] +
			                         "); minibwa's mb_map_batch_pe requires every row to be a genuine pair");
		}
	}

	std::vector<int32_t> qlens(2 * n_seq);
	std::vector<const char *> seq_ptrs(2 * n_seq);
	std::vector<const char *> name_ptrs(2 * n_seq);
	for (size_t i = 0; i < n_seq; i++) {
		qlens[2 * i] = static_cast<int32_t>(queries.sequences1[i].length());
		qlens[2 * i + 1] = static_cast<int32_t>(queries.sequences2[i].length());
		seq_ptrs[2 * i] = queries.sequences1[i].c_str();
		seq_ptrs[2 * i + 1] = queries.sequences2[i].c_str();
		name_ptrs[2 * i] = queries.read_ids[i].c_str();
		name_ptrs[2 * i + 1] = queries.read_ids[i].c_str();
	}

	mb_opt_t pe_opt = active_opt();
	pe_opt.flag |= MB_F_PE;

	std::vector<int32_t> n_hit(2 * n_seq, 0);
	mb_hit_t **hits = mb_map_batch_pe(&pe_opt, active_index(), static_cast<int32_t>(2 * n_seq), qlens.data(),
	                                  seq_ptrs.data(), n_hit.data(), tbuf_.get(), name_ptrs.data());
	if (!hits) {
		return;
	}

	for (size_t p = 0; p < n_seq; p++) {
		size_t r0 = 2 * p, r1 = 2 * p + 1;

		const mb_hit_t *primary0 = nullptr;
		const mb_hit_t *primary1 = nullptr;
		for (int32_t j = 0; j < n_hit[r0]; j++) {
			if (hits[r0][j].parent == hits[r0][j].id) {
				primary0 = &hits[r0][j];
				break;
			}
		}
		for (int32_t j = 0; j < n_hit[r1]; j++) {
			if (hits[r1][j].parent == hits[r1][j].id) {
				primary1 = &hits[r1][j];
				break;
			}
		}

		bool mate_mapped0 = (primary1 != nullptr && primary1->tid >= 0);
		bool mate_mapped1 = (primary0 != nullptr && primary0->tid >= 0);
		bool mate_rev0 = (primary1 != nullptr) && primary1->rev;
		bool mate_rev1 = (primary0 != nullptr) && primary0->rev;
		int64_t mate_tid0 = mate_mapped0 ? primary1->tid : -1;
		int64_t mate_tid1 = mate_mapped1 ? primary0->tid : -1;
		int64_t mate_pos0 = mate_mapped0 ? (primary1->ts + 1) : 0;
		int64_t mate_pos1 = mate_mapped1 ? (primary0->ts + 1) : 0;

		int32_t tlen = 0;
		if (primary0 && primary1 && primary0->tid >= 0 && primary0->tid == primary1->tid) {
			int64_t leftmost = std::min(primary0->ts, primary1->ts);
			int64_t rightmost = std::max(primary0->te, primary1->te);
			tlen = static_cast<int32_t>(rightmost - leftmost);
			if (primary0->ts > primary1->ts) {
				tlen = -tlen;
			}
		}

		// Same output gate + score-ratio filter as align_batch's single-end
		// path (see the comment there), applied independently per mate --
		// matches map-main.c, which also resets n_sec per segment, not per
		// fragment.
		int32_t n_sec0 = 0;
		for (int32_t j = 0; j < n_hit[r0]; j++) {
			const mb_hit_t &hit = hits[r0][j];
			bool is_primary = (hit.parent == hit.id);
			if (!(is_primary || n_sec0 < pe_opt.out_n)) {
				continue;
			}
			if (!is_primary) {
				const mb_hit_t &parent = hits[r0][hit.parent];
				if (parent.p && hit.p) {
					if (hit.p->dp_max < pe_opt.out_s * parent.p->dp_max) {
						continue;
					}
				} else if (hit.score < pe_opt.out_s * parent.score) {
					continue;
				}
			}
			hit_to_sam(hit, queries.read_ids[p], qlens[r0], output, 0, mate_mapped0, mate_rev0, mate_tid0, mate_pos0,
			          tlen);
			if (!is_primary) {
				n_sec0++;
			}
		}
		int32_t n_sec1 = 0;
		for (int32_t j = 0; j < n_hit[r1]; j++) {
			const mb_hit_t &hit = hits[r1][j];
			bool is_primary = (hit.parent == hit.id);
			if (!(is_primary || n_sec1 < pe_opt.out_n)) {
				continue;
			}
			if (!is_primary) {
				const mb_hit_t &parent = hits[r1][hit.parent];
				if (parent.p && hit.p) {
					if (hit.p->dp_max < pe_opt.out_s * parent.p->dp_max) {
						continue;
					}
				} else if (hit.score < pe_opt.out_s * parent.score) {
					continue;
				}
			}
			hit_to_sam(hit, queries.read_ids[p], qlens[r1], output, 1, mate_mapped1, mate_rev1, mate_tid1, mate_pos1,
			          -tlen);
			if (!is_primary) {
				n_sec1++;
			}
		}
	}

	for (size_t i = 0; i < 2 * n_seq; i++) {
		for (int32_t j = 0; j < n_hit[i]; j++) {
			free(hits[i][j].p);
		}
		free(hits[i]);
	}
	free(hits);
}

void MiniBWAAligner::hit_to_sam(const mb_hit_t &hit, const std::string &read_id, int32_t query_len,
                                SAMRecordBatch &batch, int segment_idx, bool mate_mapped, bool mate_rev,
                                int64_t mate_tid, int64_t mate_pos, int32_t tlen) const {
	bool is_paired = (segment_idx >= 0);
	bool is_unmapped = (hit.tid < 0);

	uint16_t flags = calculate_flags(hit, is_unmapped, segment_idx, mate_mapped, mate_rev);
	batch.read_ids.push_back(read_id);
	batch.flags.push_back(flags);

	AlignmentStats stats;
	if (is_unmapped) {
		batch.references.push_back("*");
		batch.positions.push_back(0);
		batch.stop_positions.push_back(0);
		batch.mapqs.push_back(0);
		batch.cigars.push_back("*");
	} else {
		batch.references.push_back(get_reference_name(hit.tid));
		batch.positions.push_back(hit.ts + 1); // 1-based
		batch.stop_positions.push_back(hit.te + 1);
		batch.mapqs.push_back(static_cast<uint8_t>(hit.mapq));
		batch.cigars.push_back(cigar_string(hit, query_len, flags, &stats));
	}

	if (is_paired && mate_mapped && mate_tid >= 0) {
		const std::string &mate_ref = get_reference_name(mate_tid);
		if (!is_unmapped && mate_ref == batch.references.back()) {
			batch.mate_references.push_back("=");
		} else {
			batch.mate_references.push_back(mate_ref);
		}
		batch.mate_positions.push_back(mate_pos);
	} else {
		batch.mate_references.push_back("*");
		batch.mate_positions.push_back(0);
	}
	batch.template_lengths.push_back(is_paired ? tlen : 0);

	batch.tag_as_values.push_back(is_unmapped ? -1 : (hit.p ? hit.p->dp_score : -1));
	batch.tag_xs_values.push_back(hit.subsc > 0 ? hit.subsc : -1);
	batch.tag_ys_values.push_back(-1); // opposite-mate score not tracked by mb_hit_t
	batch.tag_xn_values.push_back(-1); // not available from minibwa
	batch.tag_xm_values.push_back(is_unmapped ? -1 : stats.mismatches);
	batch.tag_xo_values.push_back(is_unmapped ? -1 : stats.gap_opens);
	batch.tag_xg_values.push_back(is_unmapped ? -1 : stats.gap_extends);

	// NM: O(1) formula, matches Minimap2Aligner (same field names/semantics
	// on mb_hit_t/mb_extra_t as mm_reg1_t/mm_extra_t).
	int64_t nm = hit.blen - hit.mlen + (hit.p ? hit.p->n_ambi : 0);
	batch.tag_nm_values.push_back(is_unmapped ? -1 : nm);

	std::string yt;
	if (!is_paired) {
		yt = "UU";
	} else if (mate_mapped && !is_unmapped && hit.proper_pair) {
		yt = "CP"; // concordant pair (mb_pair already decided this -- see hit.proper_pair)
	} else if (mate_mapped && !is_unmapped) {
		yt = "DP"; // discordant pair
	} else {
		yt = "UP"; // one mate unmapped
	}
	batch.tag_yt_values.push_back(yt);

	// MD tag requires mb_write_MD(), a private function (mbpriv.h) needing
	// direct l2b_t access this library does not expose publicly. Left NULL,
	// same treatment bowtie2 gives tag_sa as a documented parity stub.
	batch.tag_md_values.push_back("");
	batch.tag_sa_values.push_back("");
}

std::string MiniBWAAligner::cigar_string(const mb_hit_t &hit, int32_t query_len, uint16_t sam_flags,
                                         AlignmentStats *stats_out) const {
	if (!hit.p || hit.p->n_cigar == 0) {
		return "*";
	}

	if (hit.qs < 0 || hit.qe < 0 || hit.qs > query_len || hit.qe > query_len || hit.qs > hit.qe) {
		return "*";
	}

	// mb_map's CIGAR covers only the aligned region -- soft/hard clips must
	// be synthesized from qs/qe/qlen (documented in api-test/ex-one.c).
	int clip_front = hit.rev ? (query_len - hit.qe) : hit.qs;
	int clip_back = hit.rev ? hit.qs : (query_len - hit.qe);

	char clip_char = (sam_flags & 0x800) ? 'H' : 'S';

	std::string result;
	result.reserve((hit.p->n_cigar + 2) * 11);

	auto append_int = [&result](uint32_t val) {
		char buf[10];
		int len = 0;
		do {
			buf[len++] = '0' + (val % 10);
			val /= 10;
		} while (val);
		for (int i = len - 1; i >= 0; --i) {
			result.push_back(buf[i]);
		}
	};

	if (clip_front > 0) {
		append_int(clip_front);
		result.push_back(clip_char);
	}

	for (uint32_t i = 0; i < hit.p->n_cigar; i++) {
		uint32_t op = hit.p->cigar[i] & 0xf;
		uint32_t len = hit.p->cigar[i] >> 4;
		append_int(len);
		result.push_back(MB_CIGAR_STR[op]);
		if (stats_out) {
			switch (op) {
			case MB_CIGAR_X_MISMATCH:
				stats_out->mismatches += len;
				break;
			case MB_CIGAR_INS:
			case MB_CIGAR_DEL:
				stats_out->gap_opens += 1;
				stats_out->gap_extends += (len > 1) ? (len - 1) : 0;
				break;
			default:
				break;
			}
		}
	}

	if (clip_back > 0) {
		append_int(clip_back);
		result.push_back(clip_char);
	}

	return result;
}

uint16_t MiniBWAAligner::calculate_flags(const mb_hit_t &hit, bool is_unmapped, int segment_idx, bool mate_mapped,
                                         bool mate_rev) const {
	uint16_t flags = 0;
	bool is_paired = (segment_idx >= 0);

	if (is_paired) {
		flags |= 0x1;
		flags |= (segment_idx == 0) ? 0x40 : 0x80;
		if (mate_mapped && !is_unmapped && hit.proper_pair) {
			flags |= 0x2;
		}
		if (!mate_mapped) {
			flags |= 0x8;
		}
		if (mate_rev) {
			flags |= 0x20;
		}
	}

	if (is_unmapped) {
		flags |= 0x4;
	} else if (hit.rev) {
		flags |= 0x10;
	}

	// Secondary/supplementary: parent != id -> secondary; parent == id but
	// !sam_pri -> supplementary. Same convention as minimap2 (shared lineage
	// -- mb_set_parent/mb_set_sam_pri in map-algo.c mirror minimap2's).
	if (hit.parent != hit.id) {
		flags |= 0x100;
	} else if (!hit.sam_pri) {
		flags |= 0x800;
	}

	return flags;
}

const std::string &MiniBWAAligner::get_reference_name(int64_t tid) const {
	auto &names = active_subject_names();
	if (tid < 0 || static_cast<size_t>(tid) >= names.size()) {
		static const std::string unknown = "*";
		return unknown;
	}
	return names[static_cast<size_t>(tid)];
}

} // namespace miint
