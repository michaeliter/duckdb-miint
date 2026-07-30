#include "MiniBWAAligner.hpp"
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
	opt.best_n = config.max_secondary;
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
	align_batch(queries.read_ids, queries.sequences1, output);
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

	for (size_t i = 0; i < sequences.size(); i++) {
		int secondary_count = 0;
		for (int32_t j = 0; j < n_hit[i]; j++) {
			const mb_hit_t &hit = hits[i][j];
			bool is_primary = (hit.parent == hit.id);
			if (!is_primary) {
				if (secondary_count >= config_.max_secondary) {
					continue;
				}
				secondary_count++;
			}
			hit_to_sam(hit, read_ids[i], qlens[i], output);
		}
		for (int32_t j = 0; j < n_hit[i]; j++) {
			free(hits[i][j].p);
		}
		free(hits[i]);
	}
	free(hits);
}

void MiniBWAAligner::hit_to_sam(const mb_hit_t &hit, const std::string &read_id, int32_t query_len,
                                SAMRecordBatch &batch) const {
	bool is_unmapped = (hit.tid < 0);

	uint16_t flags = calculate_flags(hit, is_unmapped);
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

	// Single-end only in this phase: no mate.
	batch.mate_references.push_back("*");
	batch.mate_positions.push_back(0);
	batch.template_lengths.push_back(0);

	batch.tag_as_values.push_back(is_unmapped ? -1 : (hit.p ? hit.p->dp_score : -1));
	batch.tag_xs_values.push_back(hit.subsc > 0 ? hit.subsc : -1);
	batch.tag_ys_values.push_back(-1); // not applicable, single-end
	batch.tag_xn_values.push_back(-1); // not available from minibwa
	batch.tag_xm_values.push_back(is_unmapped ? -1 : stats.mismatches);
	batch.tag_xo_values.push_back(is_unmapped ? -1 : stats.gap_opens);
	batch.tag_xg_values.push_back(is_unmapped ? -1 : stats.gap_extends);

	// NM: O(1) formula, matches Minimap2Aligner (same field names/semantics
	// on mb_hit_t/mb_extra_t as mm_reg1_t/mm_extra_t).
	int64_t nm = hit.blen - hit.mlen + (hit.p ? hit.p->n_ambi : 0);
	batch.tag_nm_values.push_back(is_unmapped ? -1 : nm);

	batch.tag_yt_values.push_back("UU"); // single-end only in this phase

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

uint16_t MiniBWAAligner::calculate_flags(const mb_hit_t &hit, bool is_unmapped) const {
	uint16_t flags = 0;

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
