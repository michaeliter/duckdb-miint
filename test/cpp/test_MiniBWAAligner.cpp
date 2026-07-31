#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "MiniBWAAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace miint;

// mb_idx_build() only accepts a FASTA file path (unlike minimap2's in-memory
// mm_idx_str), so tests write a temp FASTA, build an index prefix from it,
// then clean up both -- same convention test_Minimap2Aligner.cpp uses for
// its .mmi files, adapted for minibwa's two-file (.l2b/.mbw) index.
static void build_test_index(const std::string &fasta_path, const std::vector<std::pair<std::string, std::string>> &subjects,
                             const std::string &index_prefix) {
	std::ofstream out(fasta_path);
	for (const auto &s : subjects) {
		out << '>' << s.first << '\n' << s.second << '\n';
	}
	out.close();
	int rc = mb_idx_build(fasta_path.c_str(), index_prefix.c_str(), /*sa_bit=*/4, /*n_thread=*/1, /*is_meth=*/0,
	                      /*seed=*/11);
	REQUIRE(rc == 0);
}

static void remove_index_files(const std::string &prefix) {
	std::remove((prefix + ".l2b").c_str());
	std::remove((prefix + ".mbw").c_str());
}

static SequenceRecordBatch make_query_batch(const std::string &read_id, const std::string &sequence) {
	SequenceRecordBatch batch(false);
	batch.read_ids.push_back(read_id);
	batch.comments.push_back("");
	batch.sequences1.push_back(sequence);
	batch.quals1.push_back(QualScore(""));
	return batch;
}

static SequenceRecordBatch make_paired_query_batch(const std::string &read_id, const std::string &seq1,
                                                    const std::string &seq2) {
	SequenceRecordBatch batch(true);
	batch.read_ids.push_back(read_id);
	batch.comments.push_back("");
	batch.sequences1.push_back(seq1);
	batch.sequences2.push_back(seq2);
	batch.quals1.push_back(QualScore(""));
	batch.quals2.push_back(QualScore(""));
	return batch;
}

TEST_CASE("MiniBWAConfig defaults", "[MiniBWAAligner]") {
	MiniBWAConfig config;
	REQUIRE(config.preset == "sr");
	// Matches minibwa CLI's own --outn default (emit zero secondaries), NOT
	// minimap2Config's default of 5 -- see the comment on MiniBWAConfig and
	// the out_n/out_s fix in MiniBWAAligner.cpp's history.
	REQUIRE(config.max_secondary == 0);
	REQUIRE(config.secondary_score_ratio == Catch::Approx(0.8f));
	REQUIRE(config.sa_bit == 4);
	REQUIRE(config.is_meth == false);
}

TEST_CASE("MiniBWAAligner construction does not throw", "[MiniBWAAligner]") {
	MiniBWAConfig config;
	REQUIRE_NOTHROW(MiniBWAAligner(config));
}

TEST_CASE("InitOptions never lets config override best_n", "[MiniBWAAligner]") {
	// The whole point of the out_n/out_s fix: max_secondary must map to
	// out_n only. best_n (internal retention) stays at the preset's value
	// (50 for sr/adap) regardless of what max_secondary is set to.
	MiniBWAConfig config;
	config.max_secondary = 0;
	mb_opt_t opt;
	MiniBWAAligner::InitOptions(config, opt);
	REQUIRE(opt.best_n == 50);
	REQUIRE(opt.out_n == 0);

	config.max_secondary = 7;
	MiniBWAAligner::InitOptions(config, opt);
	REQUIRE(opt.best_n == 50); // unchanged
	REQUIRE(opt.out_n == 7);
}

TEST_CASE("InitOptions rejects unknown preset", "[MiniBWAAligner]") {
	MiniBWAConfig config;
	config.preset = "not_a_real_preset";
	mb_opt_t opt;
	REQUIRE_THROWS(MiniBWAAligner::InitOptions(config, opt));
}

TEST_CASE("InitOptions clears MB_F_PE by default", "[MiniBWAAligner]") {
	// mb_opt_preset("sr"/"adap") sets MB_F_PE upstream; InitOptions must
	// clear it since MiniBWAAligner::align() decides pairing itself from
	// the query batch, not from a sticky opt flag.
	MiniBWAConfig config;
	mb_opt_t opt;
	MiniBWAAligner::InitOptions(config, opt);
	REQUIRE((opt.flag & MB_F_PE) == 0);
}

TEST_CASE("MiniBWAAligner single-end exact match alignment", "[MiniBWAAligner]") {
	std::string ref = "TTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATGCATCGATCGTAGC";
	build_test_index("data/shards/test_minibwa_se.fa", {{"reference", ref}}, "data/shards/test_minibwa_se");

	MiniBWAConfig config;
	config.preset = "sr";
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_se");

	auto query = make_query_batch("query1", ref.substr(20, 50));
	SAMRecordBatch output;
	aligner.align(query, output);

	REQUIRE(output.size() == 1);
	REQUIRE(output.references[0] == "reference");
	REQUIRE(output.positions[0] == 21); // 1-based
	REQUIRE(output.cigars[0] == "50M");
	REQUIRE(output.mapqs[0] == 60);
	REQUIRE(output.tag_nm_values[0] == 0);
	REQUIRE((output.flags[0] & 0x10) == 0); // forward strand
	REQUIRE(output.tag_md_values[0] == ""); // documented gap: no MD tag support

	std::remove("data/shards/test_minibwa_se.fa");
}

TEST_CASE("MiniBWAAligner CIGAR includes soft clips for partial query alignment", "[MiniBWAAligner]") {
	std::string ref = "TTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATGCATCGATCGTAGC";
	build_test_index("data/shards/test_minibwa_clip.fa", {{"reference", ref}}, "data/shards/test_minibwa_clip");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_clip");

	// 5bp junk prefix + 50bp exact match + 5bp junk suffix
	std::string clipped = "AAAAA" + ref.substr(20, 50) + "GGGGG";
	auto query = make_query_batch("clipped_read", clipped);
	SAMRecordBatch output;
	aligner.align(query, output);

	REQUIRE(output.size() == 1);
	REQUIRE(output.cigars[0].find('S') != std::string::npos);
	// clip + match lengths must sum to the query length
	int total = 0, num = 0;
	for (char c : output.cigars[0]) {
		if (c >= '0' && c <= '9') {
			num = num * 10 + (c - '0');
		} else {
			total += num;
			num = 0;
		}
	}
	REQUIRE(total == static_cast<int>(clipped.size()));

	std::remove("data/shards/test_minibwa_clip.fa");
}

TEST_CASE("MiniBWAAligner unmapped read produces no output row", "[MiniBWAAligner]") {
	// Matches Minimap2Aligner::align_single's existing behavior: a read with
	// zero candidate seeds emits nothing at all, not an unmapped placeholder
	// row (see align_minibwa.cpp's commit history for the full explanation).
	std::string ref = "TTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATGCATCGATCGTAGC";
	build_test_index("data/shards/test_minibwa_unmap.fa", {{"reference", ref}}, "data/shards/test_minibwa_unmap");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_unmap");

	auto query = make_query_batch("nomatch", "GGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCC");
	SAMRecordBatch output;
	aligner.align(query, output);

	REQUIRE(output.size() == 0);

	std::remove("data/shards/test_minibwa_unmap.fa");
}

TEST_CASE("MiniBWAAligner paired-end alignment reports proper pair", "[MiniBWAAligner]") {
	std::string ref = "CTAGCATGGCTAGCATGGCTAGCATGGCTATTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGAT"
	                  "CGTAGCTAGCTAGCATGCATCGATCGTAGCTAGCATCGATGCTAGCATCGATGCATCGATGCTAGCGGGCTATCGATCGATCGATGCATGCATGCAT"
	                  "GCATCGATGCTAGCTAGCTAGCATCGATGCTAGCATGCATGCATGCATCGATGCTAGCATGCATGCATGCATGC";
	build_test_index("data/shards/test_minibwa_pe.fa", {{"reference", ref}}, "data/shards/test_minibwa_pe");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_pe");

	std::string mate1 = ref.substr(30, 50);
	std::string mate2_fwd = ref.substr(200, 50);
	std::string mate2(mate2_fwd.rbegin(), mate2_fwd.rend());
	for (auto &c : mate2) {
		switch (c) {
		case 'A':
			c = 'T';
			break;
		case 'C':
			c = 'G';
			break;
		case 'G':
			c = 'C';
			break;
		case 'T':
			c = 'A';
			break;
		default:
			break;
		}
	}

	auto query = make_paired_query_batch("pair1", mate1, mate2);
	SAMRecordBatch output;
	aligner.align(query, output);

	REQUIRE(output.size() == 2);
	// mate1 (first in pair): flag 0x40 set
	size_t r1 = (output.flags[0] & 0x40) ? 0 : 1;
	size_t r2 = 1 - r1;
	REQUIRE(output.positions[r1] == 31);
	REQUIRE(output.positions[r2] == 201);
	REQUIRE((output.flags[r1] & 0x2) != 0); // proper pair
	REQUIRE((output.flags[r2] & 0x2) != 0);
	REQUIRE(output.template_lengths[r1] == 220);
	REQUIRE(output.template_lengths[r2] == -220);
	REQUIRE(output.tag_yt_values[r1] == "CP");
	REQUIRE(output.tag_yt_values[r2] == "CP");

	std::remove("data/shards/test_minibwa_pe.fa");
}

TEST_CASE("MiniBWAAligner falls back to single-end for schema-paired/data-unpaired batch", "[MiniBWAAligner]") {
	// Regression test for the read_fastx()-general-schema bug: is_paired
	// true (schema has a sequence2 column) but every row's sequence2 is
	// empty must NOT throw -- it must align as single-end.
	std::string ref = "TTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATGCATCGATCGTAGC";
	build_test_index("data/shards/test_minibwa_fallback.fa", {{"reference", ref}}, "data/shards/test_minibwa_fallback");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_fallback");

	SequenceRecordBatch query(true); // is_paired=true, but sequences2 stays empty
	query.read_ids.push_back("query1");
	query.comments.push_back("");
	query.sequences1.push_back(ref.substr(20, 50));
	query.sequences2.push_back(""); // empty, matching read_fastx()'s no-mate-file shape
	query.quals1.push_back(QualScore(""));
	query.quals2.push_back(QualScore(""));

	SAMRecordBatch output;
	REQUIRE_NOTHROW(aligner.align(query, output));
	REQUIRE(output.size() == 1);
	REQUIRE(output.positions[0] == 21);

	std::remove("data/shards/test_minibwa_fallback.fa");
}

TEST_CASE("MiniBWAAligner paired batch with a genuine mix throws", "[MiniBWAAligner]") {
	std::string ref = "TTCGATCGGATCCTAGCATGCTAGCTAGCATCGATGCTAGCTGACTGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATGCATCGATCGTAGC";
	build_test_index("data/shards/test_minibwa_mix.fa", {{"reference", ref}}, "data/shards/test_minibwa_mix");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_mix");

	SequenceRecordBatch query(true);
	query.read_ids = {"r1", "r2"};
	query.comments = {"", ""};
	query.sequences1 = {ref.substr(0, 50), ref.substr(20, 50)};
	query.sequences2 = {"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC", ""}; // mixed: one paired, one not
	query.quals1 = {QualScore(""), QualScore("")};
	query.quals2 = {QualScore(""), QualScore("")};

	SAMRecordBatch output;
	REQUIRE_THROWS(aligner.align(query, output));

	std::remove("data/shards/test_minibwa_mix.fa");
}

TEST_CASE("EnumerateContigNames handles multiple contigs", "[MiniBWAAligner]") {
	build_test_index("data/shards/test_minibwa_multi.fa",
	                 {{"contig_a", "ACGTACGTACGTACGTACGTACGTACGTACGT"}, {"contig_b", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"}},
	                 "data/shards/test_minibwa_multi");

	MiniBWAConfig config;
	MiniBWAAligner aligner(config);
	aligner.load_index("data/shards/test_minibwa_multi");

	// EnumerateContigNames needs a real mb_idx_t*, which is only reachable
	// through the aligner's private state; exercise it indirectly via
	// alignment against each contig to confirm both names round-trip.
	auto q_a = make_query_batch("qa", "ACGTACGTACGTACGTACGTACGTACGTACGT");
	auto q_b = make_query_batch("qb", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT");
	SAMRecordBatch out_a, out_b;
	aligner.align(q_a, out_a);
	aligner.align(q_b, out_b);
	REQUIRE(out_a.size() == 1);
	REQUIRE(out_a.references[0] == "contig_a");
	REQUIRE(out_b.size() == 1);
	REQUIRE(out_b.references[0] == "contig_b");

	std::remove("data/shards/test_minibwa_multi.fa");
}

TEST_CASE("Two aligners sharing SharedMiniBWAIndex produce identical results", "[MiniBWAAligner]") {
	std::string ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                  "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	build_test_index("data/shards/test_minibwa_shared.fa", {{"reference", ref}}, "data/shards/test_minibwa_shared");

	MiniBWAConfig config;
	auto shared_idx = MiniBWAAligner::BuildSharedIndex("data/shards/test_minibwa_shared", config);

	MiniBWAAligner aligner1(config);
	aligner1.attach_shared_index(shared_idx);
	MiniBWAAligner aligner2(config);
	aligner2.attach_shared_index(shared_idx);

	auto query = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
	SAMRecordBatch batch1, batch2;
	aligner1.align(query, batch1);
	aligner2.align(query, batch2);

	REQUIRE(batch1.size() == batch2.size());
	REQUIRE(batch1.size() == 1);
	REQUIRE(batch1.positions[0] == batch2.positions[0]);
	REQUIRE(batch1.cigars[0] == batch2.cigars[0]);

	std::remove("data/shards/test_minibwa_shared.fa");
}

TEST_CASE("Concurrent alignment on shared minibwa index from two threads", "[MiniBWAAligner]") {
	std::string ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                  "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	build_test_index("data/shards/test_minibwa_concurrent.fa", {{"reference", ref}}, "data/shards/test_minibwa_concurrent");

	MiniBWAConfig config;
	auto shared_idx = MiniBWAAligner::BuildSharedIndex("data/shards/test_minibwa_concurrent", config);

	auto query = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch1, batch2;
	bool ok1 = false, ok2 = false;

	std::thread t1([&]() {
		MiniBWAAligner aligner(config);
		aligner.attach_shared_index(shared_idx);
		aligner.align(query, batch1);
		ok1 = (batch1.size() == 1 && batch1.references[0] == "reference");
	});
	std::thread t2([&]() {
		MiniBWAAligner aligner(config);
		aligner.attach_shared_index(shared_idx);
		aligner.align(query, batch2);
		ok2 = (batch2.size() == 1 && batch2.references[0] == "reference");
	});
	t1.join();
	t2.join();

	REQUIRE(ok1);
	REQUIRE(ok2);
	REQUIRE(batch1.positions[0] == batch2.positions[0]);
	REQUIRE(batch1.cigars[0] == batch2.cigars[0]);

	std::remove("data/shards/test_minibwa_concurrent.fa");
}

// Clean up temporary index files created by tests above.
TEST_CASE("Cleanup temp minibwa index files", "[MiniBWAAligner]") {
	std::vector<std::string> prefixes = {
	    "data/shards/test_minibwa_se",       "data/shards/test_minibwa_clip",     "data/shards/test_minibwa_unmap",
	    "data/shards/test_minibwa_pe",       "data/shards/test_minibwa_fallback", "data/shards/test_minibwa_mix",
	    "data/shards/test_minibwa_multi",    "data/shards/test_minibwa_shared",   "data/shards/test_minibwa_concurrent"};
	for (const auto &prefix : prefixes) {
		std::remove((prefix + ".l2b").c_str());
		std::remove((prefix + ".mbw").c_str());
	}
}
