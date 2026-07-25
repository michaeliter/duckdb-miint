#include <catch2/catch_test_macros.hpp>
#include "RammapAligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace miint;

// Helper to create a single unpaired query batch
static SequenceRecordBatch make_query_batch(const std::string &read_id, const std::string &sequence) {
	SequenceRecordBatch batch(false); // unpaired
	batch.read_ids.push_back(read_id);
	batch.comments.push_back("");
	batch.sequences1.push_back(sequence);
	batch.quals1.push_back(QualScore(""));
	return batch;
}

static std::vector<AlignmentSubject> two_subjects() {
	return {
	    {"ref1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
	             "TTAAGGCC"},
	    {"ref2", "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT"
	             "AATTAA"},
	};
}

TEST_CASE("RammapConfig defaults", "[RammapAligner]") {
	RammapConfig config;
	REQUIRE(config.preset == "sr");
	REQUIRE(config.k == 0);
	REQUIRE(config.w == 0);
	REQUIRE(config.eqx == true);
	REQUIRE(config.max_secondary == 0); // 0 = use preset default
}

TEST_CASE("SharedRammapIndex builds from subjects", "[RammapAligner]") {
	RammapConfig config;
	REQUIRE_NOTHROW(std::make_shared<SharedRammapIndex>(two_subjects(), config));
}

TEST_CASE("SharedRammapIndex rejects empty subjects", "[RammapAligner]") {
	RammapConfig config;
	std::vector<AlignmentSubject> empty;
	REQUIRE_THROWS_AS(SharedRammapIndex(empty, config), std::runtime_error);
}

TEST_CASE("SharedRammapIndex rejects unrecognized preset", "[RammapAligner]") {
	RammapConfig config;
	config.preset = "not-a-real-preset";
	REQUIRE_THROWS_AS(SharedRammapIndex(two_subjects(), config), std::runtime_error);
}

TEST_CASE("RammapAligner single-end alignment - exact match", "[RammapAligner]") {
	RammapConfig config;
	auto shared_index = std::make_shared<SharedRammapIndex>(two_subjects(), config);
	RammapAligner aligner(shared_index);

	SAMRecordBatch output;
	aligner.align_single("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", output);

	REQUIRE(output.size() >= 1);
	REQUIRE(output.read_ids[0] == "query1");
	REQUIRE(output.references[0] == "ref1");
	REQUIRE(output.positions[0] == 1); // 1-based, matches SAM POS convention
	REQUIRE(output.cigars[0] != "*");
	REQUIRE(output.tag_nm_valid[0] == true);
	REQUIRE(output.tag_nm_values[0] == 0); // exact match
}

TEST_CASE("RammapAligner unmapped read", "[RammapAligner]") {
	RammapConfig config;
	auto shared_index = std::make_shared<SharedRammapIndex>(two_subjects(), config);
	RammapAligner aligner(shared_index);

	SAMRecordBatch output;
	// Random noise unlikely to align against either subject.
	aligner.align_single("query_noise", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", output);

	REQUIRE(output.size() == 1);
	REQUIRE(output.references[0] == "*");
	// A SAM CIGAR of "*" parses to zero cigar ops; sam_utils::cigar_to_string
	// (the shared htslib-backed helper reused here — see RammapAligner.cpp)
	// re-serializes zero ops as "", not "*". This matches read_alignments()'s
	// existing behavior for on-disk unmapped SAM records, not a rammap quirk.
	REQUIRE(output.cigars[0].empty());
}

TEST_CASE("RammapAligner::align batch, single-end only", "[RammapAligner]") {
	RammapConfig config;
	auto shared_index = std::make_shared<SharedRammapIndex>(two_subjects(), config);
	RammapAligner aligner(shared_index);

	auto batch = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
	SAMRecordBatch output;
	REQUIRE_NOTHROW(aligner.align(batch, output));
	REQUIRE(output.size() >= 1);
}

TEST_CASE("RammapAligner::align rejects paired-end queries", "[RammapAligner]") {
	RammapConfig config;
	auto shared_index = std::make_shared<SharedRammapIndex>(two_subjects(), config);
	RammapAligner aligner(shared_index);

	SequenceRecordBatch batch(true); // paired
	batch.read_ids.push_back("query1");
	batch.comments.push_back("");
	batch.sequences1.push_back("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
	batch.sequences2.push_back("TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA");
	batch.quals1.push_back(QualScore(""));
	batch.quals2.push_back(QualScore(""));

	SAMRecordBatch output;
	REQUIRE_THROWS_AS(aligner.align(batch, output), std::runtime_error);
}
