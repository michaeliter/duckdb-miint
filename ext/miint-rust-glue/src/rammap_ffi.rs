// FFI surface for rammap-core, consumed by MIINT's RammapAligner
// (src/RammapAligner.cpp) via src/include/rammap_ffi.h.
//
// PARITY STRATEGY: this module deliberately does NOT expose rammap's
// structured `Mapping`/`AlnResult` fields for C++ to re-serialize into SAM.
// That reconstruction-in-C++ pattern is exactly what causes align_minimap2's
// output to drift from standalone minimap2. Instead, `rammap_align_query`
// calls `rammap::align::pipeline::align_and_format_query` — the same
// function the standalone `rammap` CLI calls to produce its SAM output (see
// rammap/src/main.rs) — and returns the resulting SAM text verbatim. C++
// parses that text with MIINT's existing SAM parsing, so as long as the
// MapOptions/OutputConfig constructed here match the CLI's for a given
// preset, output is byte-identical by construction.
//
// `rammap_index_build`'s option derivation mirrors rammap/src/main.rs::run()
// exactly for the subset of CLI behavior MIINT exposes today (preset, k, w,
// eqx, best_n / "max_secondary"): apply_preset_str -> optional k/w override
// -> unconditional chn_pen_gap/chn_pen_skip recompute -> Index::build ->
// cal_mid_occ(2e-4, ...). Splice junction files (--junc-bed/--spsc) and
// paired-end alignment are not yet exposed — follow-up work.

use rammap::align::extend::AlignmentContext;
use rammap::align::index::Index;
use rammap::align::map::{MapContext, MapOptions};
use rammap::align::pipeline::{align_and_format_query, OutputConfig, ReadInfo};
use rammap::api::apply_preset_str;
use std::ffi::{c_char, CStr, CString};
use std::os::raw::c_int;
use std::ptr;
use std::slice;

/// Owns a built rammap `Index` plus the `MapOptions`/`OutputConfig` derived
/// from the preset — everything shared (read-only) across every query
/// aligned against this subject set. One instance is shared across worker
/// threads; mirrors `miint::SharedMinimap2Index`.
pub struct RammapIndexHandle {
    index: Index,
    opt: MapOptions,
    out_cfg: OutputConfig,
}

/// Per-thread scratch state (chaining/DP scratch buffers). NOT thread-safe:
/// exactly one `RammapCtxHandle` per worker thread, never shared. Mirrors
/// `miint::Minimap2Aligner`'s `mm_tbuf_t`.
pub struct RammapCtxHandle {
    align_ctx: AlignmentContext,
    map_ctx: MapContext,
}

/// Build an index from in-memory subject sequences and derive `MapOptions`
/// from `preset`.
///
/// `k_override`/`w_override`/`best_n_override` <= 0 mean "use the preset's
/// default" (matches the CLI's `Option<T>` "only override if the flag was
/// passed" semantics). `eqx` toggles `=`/`X` CIGAR ops. Returns null on
/// error (unrecognized preset, empty input, or invalid UTF-8 in a name).
///
/// # Safety
/// `names`, `seqs`, and `seq_lens` must be valid arrays of length `n`.
/// `names[i]` must be a NUL-terminated C string readable for the duration of
/// this call. `seqs[i]` must be readable for `seq_lens[i]` bytes. `preset`
/// must be a valid NUL-terminated C string. The caller owns the returned
/// pointer and must release it exactly once via `rammap_index_free`.
#[no_mangle]
pub unsafe extern "C" fn rammap_index_build(
    names: *const *const c_char,
    seqs: *const *const u8,
    seq_lens: *const usize,
    n: usize,
    preset: *const c_char,
    k_override: c_int,
    w_override: c_int,
    eqx: bool,
    best_n_override: c_int,
) -> *mut RammapIndexHandle {
    if names.is_null() || seqs.is_null() || seq_lens.is_null() || n == 0 || preset.is_null() {
        return ptr::null_mut();
    }
    let preset_str = match CStr::from_ptr(preset).to_str() {
        Ok(s) => s,
        Err(_) => return ptr::null_mut(),
    };

    // --- Mirrors rammap/src/main.rs::run(): preset, then CLI overrides ---
    let mut k = 15usize;
    let mut w = 10usize;
    let mut is_hpc = false;
    let mut opt = MapOptions::default();
    if apply_preset_str(&mut opt, &mut k, &mut w, &mut is_hpc, preset_str).is_err() {
        return ptr::null_mut();
    }
    if k_override > 0 {
        k = k_override as usize;
    }
    if w_override > 0 {
        w = w_override as usize;
    }
    if best_n_override > 0 {
        opt.filtering.best_n = best_n_override;
    }
    // Recompute chaining penalties after overrides — unconditional in the
    // CLI ("Recompute chaining penalties after all CLI overrides").
    opt.chaining.chn_pen_gap = (opt.chaining.chain_gap_scale as f64 * 0.01 * (k as f64)) as f32;
    opt.chaining.chn_pen_skip = opt.filtering.chain_skip_scale * opt.scoring.match_score as f32 * 0.01;

    let names_slice = slice::from_raw_parts(names, n);
    let seqs_slice = slice::from_raw_parts(seqs, n);
    let lens_slice = slice::from_raw_parts(seq_lens, n);

    let mut records: Vec<(String, Vec<u8>)> = Vec::with_capacity(n);
    for i in 0..n {
        let name = match CStr::from_ptr(names_slice[i]).to_str() {
            Ok(s) => s.to_string(),
            Err(_) => return ptr::null_mut(),
        };
        if seqs_slice[i].is_null() {
            return ptr::null_mut();
        }
        let seq = slice::from_raw_parts(seqs_slice[i], lens_slice[i]).to_vec();
        records.push((name, seq));
    }

    let mut index = Index::build(records, w, k, is_hpc, usize::MAX);
    index.index = 0;
    // Some presets (e.g. "sr") set a fixed mid_occ directly and skip
    // recalculation entirely — mirror the CLI's `if mid_occ == 0` guard
    // exactly (main.rs); unconditionally overwriting it here would silently
    // diverge from CLI output for those presets.
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ = index.cal_mid_occ(2e-4, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }

    // output_sam=true, do_md=true: MIINT's tag_md column is populated
    // directly from rammap's own MD:Z tag (no C++-side reconstruction).
    let out_cfg = OutputConfig {
        do_cigar: true,
        do_cs: false,
        cs_long: false,
        do_md: true,
        do_ds: false,
        eqx,
        output_sam: true,
        rg_id: None,
        split_mode: false,
    };

    Box::into_raw(Box::new(RammapIndexHandle { index, opt, out_cfg }))
}

/// # Safety
/// `handle` must be a pointer previously returned by `rammap_index_build`
/// (or null, which is a no-op), not already freed, and not in concurrent
/// use by any `rammap_align_query` call.
#[no_mangle]
pub unsafe extern "C" fn rammap_index_free(handle: *mut RammapIndexHandle) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

/// Allocate per-thread scratch state. The caller owns the returned pointer
/// and must release it exactly once via `rammap_ctx_free`.
#[no_mangle]
pub extern "C" fn rammap_ctx_new() -> *mut RammapCtxHandle {
    Box::into_raw(Box::new(RammapCtxHandle {
        align_ctx: AlignmentContext::new(),
        map_ctx: MapContext::new(),
    }))
}

/// # Safety
/// `handle` must be a pointer previously returned by `rammap_ctx_new` (or
/// null) and not already freed.
#[no_mangle]
pub unsafe extern "C" fn rammap_ctx_free(handle: *mut RammapCtxHandle) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

/// Align one single-end read and return its formatted SAM record(s)
/// (newline-terminated, 0+ lines — an unmapped read still produces one
/// unmapped SAM line, matching CLI behavior) exactly as
/// `rammap::align::pipeline::align_and_format_query` produces them.
/// Returns null on invalid arguments or invalid UTF-8 in `qname`.
///
/// # Safety
/// `index` must be a valid pointer from `rammap_index_build`, readable
/// (shared) concurrently from multiple threads. `ctx` must be a valid
/// pointer from `rammap_ctx_new` used by exactly one thread at a time.
/// `qname` must be a valid NUL-terminated C string; `qseq` must be readable
/// for `qseq_len` bytes. The returned pointer must be released exactly once
/// via `rammap_string_free`.
#[no_mangle]
pub unsafe extern "C" fn rammap_align_query(
    index: *mut RammapIndexHandle,
    ctx: *mut RammapCtxHandle,
    qname: *const c_char,
    qseq: *const u8,
    qseq_len: usize,
) -> *mut c_char {
    if index.is_null() || ctx.is_null() || qname.is_null() || qseq.is_null() {
        return ptr::null_mut();
    }
    let idx = &mut *index;
    let c = &mut *ctx;

    let qname_str = match CStr::from_ptr(qname).to_str() {
        Ok(s) => s,
        Err(_) => return ptr::null_mut(),
    };
    let qseq_slice = slice::from_raw_parts(qseq, qseq_len);

    let read = ReadInfo {
        qname: qname_str,
        qseq: qseq_slice,
        qual: None,
        comment: None,
        n_seg: 1,
        seg_idx: 0,
    };

    let (sam_text, _stats) = align_and_format_query(
        &idx.opt,
        &idx.index,
        &read,
        &mut c.align_ctx,
        &mut c.map_ctx,
        None, // junc_db: splice junction BED/SPSC not yet exposed (follow-up)
        None, // jump_db: splice jump extension not yet exposed (follow-up)
        &idx.out_cfg,
    );

    match CString::new(sam_text) {
        Ok(s) => s.into_raw(),
        // SAM text should never contain interior NULs; fail loud rather
        // than silently truncate if it somehow does.
        Err(_) => ptr::null_mut(),
    }
}

/// # Safety
/// `s` must be a pointer previously returned by `rammap_align_query` (or
/// null) and not already freed.
#[no_mangle]
pub unsafe extern "C" fn rammap_string_free(s: *mut c_char) {
    if !s.is_null() {
        drop(CString::from_raw(s));
    }
}
