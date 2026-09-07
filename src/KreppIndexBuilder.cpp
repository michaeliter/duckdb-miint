#include "KreppIndexBuilder.hpp"

#include <mutex>

// krepp's own headers. krepp.hpp is deliberately not included: it pulls in
// CLI11 and the command-line layer, neither of which applies here. Everything
// the builder needs moved into index.hpp/cpp upstream (bo1929/krepp#13) for
// exactly this reason.
#include "common.hpp"
#include "index.hpp"

namespace miint {

void InstallKreppErrorHandler() {
	static std::once_flag installed;
	std::call_once(installed, [] {
		set_error_handler([](const std::string &message, int code) { throw KreppFatalError(message, code); });
	});
}

void BuildKreppIndex(const KreppIndexOptions &options) {
	IndexConfig config;
	config.input = options.input_map_path;
	config.index_dir = options.index_dir;
	config.nwk_path = options.newick_path;
	config.k = options.k;
	config.w = options.w;
	config.h = options.h;
	config.m = options.m;
	config.r = options.r;
	config.frac = options.frac;
	config.sdust_t = options.sdust_t;
	config.sdust_w = options.sdust_w;

	// The order krepp's main uses for its `index` subcommand, unchanged. It is
	// not interchangeable: set_nrows/set_lshf seed the hash function the whole
	// build hangs off, read_input_file decides per-sequence vs per-file mode,
	// and obtain_build_tree consults that decision.
	IndexMultiple index(config);
	index.set_nrows();
	index.set_lshf();
	index.read_input_file();
	index.obtain_build_tree();
	index.build_index();
	index.save_index();
}

} // namespace miint
