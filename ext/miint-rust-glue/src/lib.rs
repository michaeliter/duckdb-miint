// Umbrella staticlib for rype + sylph. The crate body is intentionally empty:
// cargo treats both upstream crates as rlib dependencies, rustc emits one
// staticlib that links them through a single Rust std, and the `#[no_mangle]
// extern "C"` symbols in rype/sylph survive into the archive's export table.
// `extern crate ... as _` forces the linkage without exposing names.

extern crate rype as _;

#[cfg(feature = "with-sylph")]
extern crate sylph as _;

// Unlike rype/sylph (which each ship their own arrow-ffi extern "C" surface
// upstream), rammap-core has no FFI of its own — this crate owns it. See
// rammap_ffi.rs for why (parity with the standalone rammap CLI).
#[cfg(feature = "with-rammap")]
mod rammap_ffi;
