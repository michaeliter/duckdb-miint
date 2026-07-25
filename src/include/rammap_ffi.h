#pragma once
/*
 * Hand-written C declarations for the rammap FFI surface implemented in
 * ext/miint-rust-glue/src/rammap_ffi.rs (built into libmiint_rust_glue.a,
 * linked via the `rype` CMake target — see the comment near that target's
 * declaration in CMakeLists.txt for why there is no separate `rammap`
 * target). Mirrors the pattern used for ext/sylph/sylph.h.
 *
 * See rammap_ffi.rs for the parity rationale: these functions route through
 * rammap-core's own `pipeline::align_and_format_query` (the same function
 * the standalone `rammap` CLI calls), returning SAM text verbatim rather
 * than structured fields for C++ to re-serialize.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct RammapIndexHandle RammapIndexHandle;
typedef struct RammapCtxHandle RammapCtxHandle;

/* eqx: nonzero = use =/X CIGAR ops. k_override/w_override/best_n_override:
 * <= 0 means "use the preset's default". Returns NULL on error (bad preset
 * name, empty input, or invalid UTF-8 in a subject name). */
RammapIndexHandle *rammap_index_build(const char *const *names, const uint8_t *const *seqs, const size_t *seq_lens,
                                      size_t n, const char *preset, int32_t k_override, int32_t w_override,
                                      int eqx, int32_t best_n_override);

/* Releases an index built by rammap_index_build. NULL is a no-op. Must not
 * be called while any rammap_align_query call against it is in flight. */
void rammap_index_free(RammapIndexHandle *handle);

/* Allocates per-thread scratch state. NOT thread-safe: exactly one context
 * per worker thread, never shared. */
RammapCtxHandle *rammap_ctx_new(void);

/* Releases a context allocated by rammap_ctx_new. NULL is a no-op. */
void rammap_ctx_free(RammapCtxHandle *handle);

/* Aligns one single-end read and returns its formatted SAM record(s)
 * (newline-terminated, 0+ lines) exactly as rammap-core's pipeline
 * produces them. Returns NULL on invalid arguments. The returned pointer
 * must be released with rammap_string_free. `index` may be shared
 * (read-only) across threads; `ctx` must not be used concurrently. */
char *rammap_align_query(RammapIndexHandle *index, RammapCtxHandle *ctx, const char *qname, const uint8_t *qseq,
                         size_t qseq_len);

/* Releases a string returned by rammap_align_query. NULL is a no-op. */
void rammap_string_free(char *s);

#ifdef __cplusplus
}
#endif
