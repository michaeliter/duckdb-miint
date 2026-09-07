#pragma once

/*
 * KreppIndexBuilder - C++ wrapper for krepp index construction.
 *
 * The mirror of KreppPlacer: that one loads an index krepp's CLI wrote, this
 * one writes an index KreppPlacer can load, from the same linked v0.9.1. Same
 * PIMPL discipline for the same reason - krepp's headers reach exactly one
 * translation unit, because their klib basenames collide with minimap2's and
 * their seq_nt4_table is a different table of a different length.
 *
 * INPUT IS FILES, NOT MEMORY, AND THAT IS KREPP'S SHAPE
 * ====================================================
 * IndexMultiple::read_input_file takes a path and nothing else. Handed a FASTA
 * it switches to "per sequence" mode, where every record is its own reference -
 * and obtain_build_tree then REFUSES a guide tree in that mode outright
 * (index.cpp: "A guide tree (-t) is incompatible with a per sequence
 * indexing."). So an index with a backbone can only be built the other way:
 * a TSV mapping each reference name to a file of its own.
 *
 * That is why the caller materialises one FASTA per reference. It is not a
 * convenience; it is the only shape that keeps the backbone.
 *
 * ERRORS
 * ======
 * krepp reports fatal errors with error_exit, which is std::exit - a malformed
 * tree would take the DuckDB process down with no SQL error at all.
 * InstallKreppErrorHandler redirects those into KreppFatalError. It is safe
 * here because this build compiles krepp without -fopenmp (CMakeLists.txt sets
 * _WOPENMP=0 and adds no flag), so krepp's `#pragma omp` directives are ignored
 * and an exception never unwinds out of a structured block.
 */

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>

namespace miint {

// What krepp's error_exit becomes once InstallKreppErrorHandler has run.
// `code` is the exit status krepp would have used; it is EXIT_FAILURE for
// every call site in the library sources, and carried anyway so a caller can
// tell them apart if that ever stops being true.
class KreppFatalError : public std::runtime_error {
public:
	KreppFatalError(const std::string &message, int code) : std::runtime_error(message), code_(code) {
	}
	int code() const noexcept {
		return code_;
	}

private:
	int code_;
};

// Install the handler that turns krepp's fatal errors into KreppFatalError.
// Idempotent, and safe to call from several threads.
//
// Called once at extension load rather than lazily at first use, which is what
// krepp's contract asks for: error_exit reads the handler without a lock, so
// installing it while another thread might already be inside krepp is a data
// race. Loading the extension is the one moment when no krepp call can be in
// flight.
void InstallKreppErrorHandler();

// Everything krepp's IndexConfig needs, minus the CLI. Defaults match krepp's
// own; `w` and `h` are optional so "not given" (derive from k) stays distinct
// from "given as 0" (invalid, and krepp says so).
struct KreppIndexOptions {
	std::string index_dir;      // directory to create and write the index into
	std::string input_map_path; // TSV: reference name <TAB> path to its FASTA
	std::string newick_path;    // backbone tree; empty means krepp generates one
	uint8_t k = 29;
	std::optional<uint8_t> w; // unset: k + 6
	std::optional<uint8_t> h; // unset: k - 16
	uint32_t m = 4;
	uint32_t r = 1;
	bool frac = true;
	uint32_t sdust_t = 0;
	uint32_t sdust_w = 0;
};

// Build the index. Runs krepp's own index pipeline in the order krepp's main
// runs it, and writes the same files to `index_dir`.
//
// Throws KreppFatalError for anything krepp itself rejects - an unreadable
// input, a malformed tree, an out-of-range parameter. Single-threaded: krepp's
// num_threads global is only ever assigned by its CLI layer, and this build has
// no OpenMP regardless.
void BuildKreppIndex(const KreppIndexOptions &options);

} // namespace miint
