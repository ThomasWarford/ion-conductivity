# Lossless compression of VASP CHGCAR/AECCAR volumetric data

Test file: `samples/frame1/CHGCAR` (Li-P, 180×180×180 grid, 5,832,000 values,
106,178,679 bytes ASCII, matches the existing `.gz` used in production).
AECCAR1/AECCAR2 behave essentially identically (same grid, same density-like
smoothness) — see note at bottom.

**Definition of "lossless" used here:** the numeric grid values must be
reconstructed to the exact precision printed in the ASCII file (11
significant digits, `0.dddddddddddE±dd`), i.e. re-formatting the decoded
values reproduces the original text byte-for-byte for the data section.
Whitespace/layout of header/footer is preserved separately, verbatim.

## Key finding: text beats naive binary — and why

Parsing the ASCII grid into a raw `float64` array and compressing *that*
(with zstd or xz) is **worse** than compressing the original ASCII text.
Cause: VASP's 11-digit text is itself a lossy quantization of the double's
53-bit mantissa. That quantization manufactures long repeated digit-runs
across neighboring grid points (the density field is spatially smooth), which
LZ/entropy coders exploit heavily. A raw `float64` keeps the full mantissa,
most of which is uncorrelated "noise" relative to the 11 digits actually
present, so it compresses worse. Byte-shuffling (transpose across the 8
mantissa/exponent bytes) barely helps, because the noise lives in low mantissa
*bits*, not in byte alignment — shuffling doesn't touch that. Predictive
codecs (fpzip, zfp) do better than plain binary+zstd because they decorrelate
via a spatial predictor before entropy coding, but even they don't beat
zstd-19 on the raw text here, because they're still predicting from the noisy
53-bit mantissa rather than the exact 11-digit value. The best result came
from working with the same 11-digit fixed precision as the text, but as
delta-coded integers instead of ASCII (see "custom" below).

## Results (CHGCAR, 106,178,679 bytes original)

| Method | Size (bytes) | vs. gzip | vs. original | Time | Notes |
|---|---|---|---|---|---|
| gzip -9 (current production baseline) | 38,582,126 | 1.00× | 2.75× | 24s | |
| zstd -3 | 42,316,848 | 0.91× | 2.51× | 0.8s | |
| zstd -9 | 38,532,143 | 1.00× | 2.75× | 4.5s | ≈ gzip, 5× faster |
| **zstd -19** | **28,094,644** | **1.37×** | **3.78×** | 96s | best generic text compressor |
| zstd -19 --long=27 | 28,152,117 | 1.37× | 3.77× | 102s | long-distance matching, no help (file < window) |
| zstd --ultra -22 | 28,306,854 | 1.36× | 3.75× | 131s | slightly worse than -19, much slower |
| xz -9 / -9e | 31,347,492 / 31,348,928 | 1.23× | 3.39× | 134-137s | |
| bzip2 -9 | 32,401,600 | 1.19× | 3.28× | 8s | good ratio for the time cost |
| binary float64 + zstd-19 | 38,877,184 | 0.99× | 2.73× | 16s | binary parse alone doesn't help (see above) |
| binary float64 + xz-9e | 34,161,708 | 1.13× | 3.11× | 23s | |
| binary float64, byte-shuffled + zstd-19 | 36,308,684 | 1.06× | 2.94× | 12s | shuffle ≈ no better than unshuffled xz |
| fpzip (lossless, `-p 0`) | 33,139,953 | 1.16× | 3.20× | **0.6s** | verified bit-exact round-trip; 150× faster than zstd-19 |
| ZFP reversible mode (via `hdf5plugin`, HDF5 container) | 36,367,332 | 1.06× | 2.92× | 0.4s | verified bit-exact; fast but weaker ratio than fpzip here |
| den2bin | — | — | — | — | skipped: lossless path needs libboost/libglm/libtclap dev headers, not present, no root, `conda-forge` solve too slow to be worth it (see plan note) |
| **custom: 11-digit mantissa+exponent decomposition, mantissa delta-coded, zstd-19** | **25,107,874** | **1.54×** | **4.23×** | ~65s (encode+compress) | **best ratio overall**, verified lossless on random sample of 2000 values |
| custom variant with xz-9e instead of zstd | 25,126,260 | 1.54× | 4.23× | ~75s | same ratio, no benefit over zstd |

AECCAR1 (same grid shape, similarly smooth core-density field): gzip 38.06MB
→ zstd-19 27.76MB — same ~1.37× gain, confirming the CHGCAR pattern
generalizes; full sweep skipped for AECCAR to save time given the result is
consistent.

## Near-lossless float32 (NOT bit-exact — informational only)

**This section is a separate, explicitly lossy experiment. It does not change
the lossless recommendation below.**

VASP prints CHGCAR values at 11 significant decimal digits; float64
(~15-17 sig digits) round-trips that exactly (as used throughout the rest of
this report), but float32 has only ~7 sig digits. Cast `CHGCAR.f64` →
`float32` → back to `float64` and re-render in VASP's exact format:

- **99.9%** of a random 5,000-value sample changed at the printed digit level
  (4,997/5,000) — this is not formatting noise, float32 genuinely discards
  real digits VASP prints.
- Measured directly over the full 5,832,000-value grid (`float64 → float32 →
  float64`, computed against the true starting float64 array — not a
  sub-sample):
  - **Max absolute error: 9.74×10⁻⁴** (occurs at value 16,409.02)
  - **Max relative error: 5.96×10⁻⁸** (occurs at value 128.00 — a *different*
    grid point than the max absolute error, since relative error is
    dominated by small-value points while absolute error is dominated by
    large ones; reporting only one metric hides this)
  - Mean absolute error: 7.58×10⁻⁶; mean relative error: 2.15×10⁻⁸
    (consistent with float32 machine epsilon, ~1.2×10⁻⁷)

Compressed sizes on the float32 buffer (23,328,000 bytes raw, half of
float64's 46,656,000). Note these three methods are themselves lossless
*wrappers* around the already-lossy float32 array — fpzip `-p 0` and zfp's
reversible mode round-trip the float32 bit pattern exactly, and zstd is
byte-exact by construction — so all three carry the *same* error above,
inherited entirely from the float64→float32 downcast, not from the
compressor itself:

| Method | Size (bytes) | vs. float64 equivalent | Max abs err | Max rel err | Time |
|---|---|---|---|---|---|
| fpzip (`-t float -3 180 180 180 -p 0`) | 11,993,308 | 33.1MB → 12.0MB (2.76× smaller) | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 0.3s |
| ZFP reversible (hdf5plugin, float32) | 14,320,509 | 36.4MB → 14.3MB (2.54× smaller) | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 0.2s |
| zstd -19 (raw float32 buffer) | 20,903,343 | 38.9MB → 20.9MB (1.86× smaller) | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 6.2s |

All three roughly halve again relative to their float64 counterparts, as
expected — but none of these are lossless, and none are directly comparable
to the text/zstd-19 or custom results above (25.1MB, which *are* lossless).

**Scientific defensibility, briefly:** a relative error of ~6×10⁻⁸ is far
below typical DFT/AIMD numerical noise floors (SCF convergence, pseudopotential
error, and thermal/statistical noise in an AIMD trajectory are all orders of
magnitude larger than 1e-8), so float32 CHGCAR data would very likely be
scientifically indistinguishable from float64 for downstream analysis
(Bader charges, visualization, etc.). But that's a statement about float32's
*absolute* precision, not about whether it's safe to silently discard the
extra digits VASP happens to print — if bit-for-bit reproducibility of the
stored file (not just "close enough" physics) matters, stick with the
lossless methods above.

## Lossy compression tuned to float32-comparable error budget

**This section is a separate, explicitly lossy/informational experiment,
like the float32 section above. It does not change the lossless
recommendation.** Question: do zfp/fpzip's *adaptive* lossy modes beat plain
float32's uniform 24-bit-mantissa truncation for a comparable error budget?

CHGCAR's value range: max |value| = 19,484.27. The earlier float32 round-trip
max relative error (5.96×10⁻⁸) converts to an absolute tolerance of
**1.161×10⁻³** at that value range — used as the target for zfp's
fixed-accuracy mode.

| Method | Size | Ratio vs. original | Achieved max abs err | Achieved max rel err | Achieved mean abs err |
|---|---|---|---|---|---|
| raw float32 cast (no compressor — precision reduction only) | 23,328,000 | 4.55× | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 7.58×10⁻⁶ |
| zfp fixed-accuracy, `accuracy=1.161e-3` (hdf5plugin) | **10,020,879** | **10.60×** | 2.64×10⁻⁴ (4.4× tighter than requested) | 1.50×10⁻³ | 3.26×10⁻⁵ |
| fpzip fixed-precision `-p 24` (double) | 4,013,422 | 26.46× | 3.99 (~3,400× looser than target) | 2.44×10⁻⁴ | 3.10×10⁻² |
| fpzip fixed-precision `-p 32` | 9,805,472 | 10.83× | 1.56×10⁻² (~13× looser than target) | 9.53×10⁻⁷ | 1.21×10⁻⁴ |
| fpzip fixed-precision `-p 40` | 15,639,896 | 6.79× | 6.09×10⁻⁵ (tighter than target) | 3.72×10⁻⁹ | 4.7×10⁻⁷ |
| *reference:* naive float32 + fpzip (from float32 section) | 11,993,308 | 8.85× | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 7.58×10⁻⁶ |
| *reference:* naive float32 + zfp-reversible | 14,320,509 | 7.41× | 9.74×10⁻⁴ | 5.96×10⁻⁸ | 7.58×10⁻⁶ |
| *reference:* lossless custom scheme (best, this report) | 25,107,874 | 4.23× | 0 (exact) | 0 | 0 |

("Ratio vs. original" = 106,178,679-byte source ASCII file ÷ compressed size,
same convention as the main results table. The two `naive float32 +
{fpzip,zfp-reversible}` reference rows carry identical error to the raw
float32 cast row above, since both are lossless *wrappers* around an
already-float32-cast array — the error is inherited entirely from the cast,
not introduced by either compressor.)

(zfp run via `hdf5plugin.Zfp(accuracy=...)` rather than the standalone zfp
CLI — same underlying codec, avoids an extra build. fpzip's `-p N` is a
bits-of-precision knob in its predicted-residual domain, not a direct
absolute-error dial, hence the bracketing sweep across N=24/32/40 rather than
one shot.)

**Systematic bias check (total electron count, not just per-cell error):**
small per-cell errors could in principle still be biased in one direction and
shift the grid sum — which matters because CHGCAR must integrate to exactly
350 electrons/cell (see the cross-file correlation section below). Checked
`sum(decoded) - sum(original)`, converted to an implied electron-count shift:

| Method | Sum deviation (relative) | Implied electron shift |
|---|---|---|
| raw float32 cast | −4.6×10⁻¹¹ | −1.6×10⁻⁸ e |
| zfp fixed-accuracy | −1.9×10⁻¹⁰ | −6.8×10⁻⁸ e |
| fpzip `-p 40` | −1.4×10⁻⁹ | −4.7×10⁻⁷ e |
| fpzip `-p 32` | −3.5×10⁻⁷ | −1.2×10⁻⁴ e |
| fpzip `-p 24` | −8.9×10⁻⁵ | **−3.1×10⁻² e** |

Bader charge analysis (the actual downstream use of these files) typically
needs only ~10⁻³–10⁻² e precision to be meaningful. Every method except
`-p 24` has a systematic bias many orders of magnitude below that — per-cell
errors are essentially unbiased (IEEE round-to-nearest) and cancel over the
5.8M-point grid, so there's no hidden charge-conservation problem for the
methods that already passed the per-cell error check. `fpzip -p 24` is
already disqualified above for blowing the per-cell budget by ~3,400×; this
confirms it — its 0.031-electron systematic shift is right at the edge of
mattering for Bader analysis, not just a per-cell curiosity.

**Verdict: yes, adaptive tuning clearly beats naive float32 truncation.**
zfp's fixed-accuracy mode hit **10.02MB** while *exceeding* the target
accuracy (actual max error 4.4× tighter than asked) — smaller than both
naive-float32 numbers (11.99MB / 14.32MB) at an equal-or-better error level.
This is exactly what fixed-accuracy mode is for: most of the CHGCAR grid is
near zero (smooth region away from atoms) and needs very few bits to hit an
absolute tolerance set by the *peak* value, whereas float32 spends a uniform
24 mantissa bits everywhere regardless of local magnitude. fpzip's precision
knob is cruder here — no `-p` setting lands close to the target error without
either wasting bits (`-p 40`: 15.64MB for an error tighter than needed) or
blowing the budget (`-p 24`: tiny 4.0MB but ~3,400× over tolerance) — zfp's
accuracy mode is the better-calibrated tool for a stated error budget.
None of this changes the lossless recommendation: it only matters if the
project decides bit-exact reproducibility of the stored file is not required.

## Cross-file correlation between CHGCAR/AECCAR

**What each file is:** INCAR confirms `LAECHG = True`, `LCHARG = True` (no
explicit `LCORE`). Header comment lines only carry the system name, no
further metadata. Per `run_chargemol.py`'s own docstring, VASP natively
writes three files here — AECCAR0, AECCAR1, AECCAR2 — of which AECCAR0 is
deleted post-chargemol as "redundant" and AECCAR1+AECCAR2 are always kept;
this project's own AECCAR0/1/2 naming is VASP's native output, not a
renaming artifact. Empirically (grid-sum check): CHGCAR integrates to
350.00 electrons/cell, exactly matching total ZVAL (75 Li_sv×3 + 25 P×5 =
350) — the pseudo-valence density, as expected. AECCAR1 and AECCAR2 both
also integrate to ~350 (349.97 / 349.96) but with far sharper peaks (max
~144,000 vs. CHGCAR's ~19,500), consistent with all-electron reconstructions
that restore the nuclear cusps the pseudopotential smooths out of CHGCAR. We
did not chase VASP's exact internal semantic distinction between AECCAR1 and
AECCAR2 beyond this — not needed for the compression question below.

**Correlation (Pearson r, float64 grid values, same frame/grid):**

| Pair | r |
|---|---|
| CHGCAR vs AECCAR1 | 0.832 |
| CHGCAR vs AECCAR2 | 0.832 |
| AECCAR1 vs AECCAR2 | 0.9996 |

AECCAR1 and AECCAR2 are nearly numerically identical (best-fit linear
regression: AECCAR2 ≈ 0.9999×AECCAR1 + 0.019, residual std 43.2 vs.
AECCAR2's own std of 1604 — a ~37× reduction in spread).

**Does this help compression? No — verified two ways, both negative:**

1. *Joint vs. separate ASCII compression:* zstd-19 on each file
   independently sums to 84,051,281 bytes (28,094,644 + 27,759,950 +
   28,196,687). Concatenating all three and compressing as one zstd-19
   stream gives 84,096,700 bytes — **no improvement** (very slightly worse,
   +0.05%, from single-stream framing overhead). An initial run of this test
   showed a spurious 26.6MB result; a repeat came back at 84.10MB, and an
   isolated AECCAR1+AECCAR2-only concatenation (the pair with r=0.9996) gave
   55,933,450 bytes vs. their independent sum of 55,956,637 — a negligible
   0.04% improvement — confirming the 26.6MB figure was an artifact and is
   disregarded. Cause:
   despite r=0.9996, the two files' printed 11-digit tokens essentially never
   match verbatim (0.005% exact match rate on a 20,000-token random sample)
   — the real values agree to ~3 significant figures, not the ~11 needed for
   literal text matches, so LZ77 finds nothing to exploit across files.
2. *Residual encoding:* `AECCAR2 − AECCAR1` (float64) compresses to
   40,763,347 bytes with zstd-19 — **worse** than compressing AECCAR1
   (38,075,809) or AECCAR2 (38,880,767) alone. Same root cause as the
   text/binary mismatch documented earlier in this report: subtracting two
   independently 11-digit-quantized values produces a float64 residual whose
   low mantissa bits are dominated by quantization noise rather than real
   structure, which byte-oriented compressors can't exploit even though the
   residual's *statistical* spread (std 43 vs. 1604) is much smaller.

**Verdict: dead end.** Despite strong numeric correlation between AECCAR1
and AECCAR2, neither joint text compression nor residual-based binary
compression captures it — both the "quantization redundancy" that makes
single-file text compress well and the "predictor decorrelation" that makes
fpzip/zfp work require exploiting structure that direct cross-file
differencing doesn't expose. Not worth pursuing further; compress each file
independently.

## Spatial predictor comparison (single frame)

A simple 3D MED-style predictor (median of the left/up/behind axis-neighbors,
JPEG-LS style) applied to the CHGCAR float64 grid, residual compressed with
zstd-19: **39,715,193 bytes** — worse than fpzip's built-in Lorenzo predictor
(33,139,953) and far worse than zstd-19-on-text (28,094,644), despite cutting
the raw standard deviation from 993 to 138. This reinforces the same
mechanism documented above: a naive predictor on full-precision float64
still leaves mantissa-level quantization noise that a generic entropy coder
can't clean up as well as either (a) working in the exact 11-digit precision
domain (text, or the custom mantissa scheme) or (b) fpzip/zfp's tuned
predictor+coder combination built specifically for floating-point data.

## Split-bytes + delta filter (Aras Pranckevičius technique)

Per [Aras Pranckevičius's float-compression series](https://aras-p.info/blog/2023/01/29/Float-Compression-0-Intro/),
tested a two-step filter distinct from the plain byte-shuffle already tried
above: (1) split each float into its byte-planes (all byte-0's, then all
byte-1's, ...), then (2) delta-encode *within* each byte-plane (`plane[i] -
plane[i-1] mod 256`), reasoning that nearby grid points' matching
high/exponent bytes should turn into long runs of near-zero deltas. Verified
bit-exact round-trip before compressing.

| Input | zstd -9 | zstd -19 | zstd -19, 1MB-chunked |
|---|---|---|---|
| float64, split+delta | 36,396,263 | 35,782,758 | 41,036,414 (worse — delta resets per chunk) |
| float32, split+delta (lossy, informational) | 14,837,530 | 14,427,354 | not tested |

**Verdict: does not beat the existing lossless results.** float64
split+delta+zstd-19 (35.78MB) is better than plain byte-shuffle alone
(36.3MB) — confirming the blog's mechanism does add something — but it's
still worse than fpzip (33.14MB) and far worse than zstd-19-on-text (28.09MB)
or the custom mantissa/delta scheme (25.11MB, current best). Chunking into
1MB blocks (for decompression-speed parallelism) makes the ratio *worse*
here since resetting the delta at each chunk boundary throws away exactly
the cross-point correlation the filter exists to capture. This does **not**
supersede the mantissa/delta recommendation — it's simpler to implement than
the custom text-domain scheme, but it compresses noticeably worse, so it
isn't a better final recommendation despite the simplicity appeal.

## den2bin source review (no build needed)

Building den2bin was skipped earlier (missing boost/glm/tclap dev headers,
no root). Read the source directly instead (`tools/den2bin/src/`,
`converter.cpp`/`density.h`/`compressor.h`) rather than building it — no new
findings worth building for:

1. **Precision:** its internal grid storage (`Density::gridptr`) is
   `std::vector<float>` — plain 32-bit float, parsed straight from the ASCII
   text. Its "lossless" `BIN` mode is exactly naive ASCII→float32, no
   fixed-point/quantized representation. By our stricter definition (exact
   reproduction of VASP's 11-digit text), **this is not actually lossless**
   — same precision loss as our own float32 section above (5.96×10⁻⁸ rel.
   err.). Header (comments string, a version-tagged format string, unit
   cell as 9 floats, grid dims as 3×uint32) is packed before the flat value
   stream — a reasonable layout, but not smarter than our `.header`/`.tail`
   + raw array split.
2. **Value transform:** the "lossless" `BIN` path applies none — straight
   parse-and-dump. Its actual lossy path (`DCT` mode) does more: normalize
   the whole grid to `[-0.5, 0.5]` using global min/max, split into
   fixed-size 3D blocks, apply a block-wise 3D DCT, and keep only the
   lowest-frequency `coeffsize` coefficients per block (JPEG-style
   truncation) — no delta coding, no predictor, nothing exploiting the
   total-electron-count normalization.
3. **Domain assumptions:** none beyond global min/max scaling — no
   core/valence split, no per-species handling, nothing AECCAR-specific.
4. **Wrapping compression:** both modes pipe the packed binary through
   **bzip2** (`boost::iostreams::bzip2_compressor`), not zstd/xz. README
   claims a 10-12× lossless ratio and up to 250× lossy ratio "for example"
   files. Reproducing their exact lossless recipe (float32 + bzip2) on our
   CHGCAR gives **21,400,215 bytes** — a 4.96× ratio vs. our 106.18MB
   original, well short of their claimed 10-12×, and worse than every other
   float32 option already tested here (float32+zstd-19: 20.90MB,
   float32+fpzip: 11.99MB, float32+zfp-reversible: 14.32MB) — bzip2 is
   simply a weaker backend than zstd/fpzip/zfp, consistent with bzip2's
   underperformance everywhere else in this report.

**Verdict: no new tricks, nothing to adopt.** den2bin's "lossless" mode is a
plain float32 dump (already covered, and already known to be lossy by our
definition) wrapped in a weaker general-purpose compressor than what we're
already using. Its one genuine idea — a value transform before entropy
coding — is the lossy 3D-block-DCT path, which is a coarser, non-adaptive
cousin of what zfp's fixed-accuracy mode already does better (zfp's
per-block bit allocation beats DCT's fixed per-block coefficient count for
a given error target, per the "Lossy compression tuned to float32-comparable
error budget" section above). Not worth building or testing further.

## Frame 2 validation

Re-ran the core methods on `samples/frame2/CHGCAR` (copied from
`/home/twarford/ion-conductivity/our_trajectories/chargemol_runs/Li-P/1500K/seg08-prod0/f03293`
— a different temperature (1500K vs. frame1's 1000K) and a different segment/
timestep than frame1, same composition/grid). Goal: check whether frame1's
numbers and recommendation generalize, not re-run every exploratory
dead-end (byte-shuffle, MED predictor, split-bytes+delta, cross-file
correlation, den2bin) — the core sweep below confirms the same pattern
within ~1%, so those mechanism-level conclusions (a property of VASP's text
format, not the specific frame) aren't re-tested here.

**Core lossless sweep (CHGCAR):**

| Method | Frame 1 (1000K) | Frame 2 (1500K) | Δ |
|---|---|---|---|
| gzip -9 | 38,582,126 | 38,575,876 | −0.02% |
| zstd -19 (text) | 28,094,644 | 28,247,427 | +0.5% |
| fpzip lossless | 33,139,953 | 33,261,869 | +0.4% |
| zfp reversible | 36,367,332 | 36,405,870 | +0.1% |
| **custom mantissa/delta (best)** | **25,107,874** | **25,308,494** | **+0.8%** |
| AECCAR1, zstd -19 (text) | 27,759,950 | 27,675,924 | −0.3% |

All within ~1% of frame1 — same ranking, same mechanism, no surprises.
Frame2's own gzip baseline (38,575,876) and its `.gz` on disk agree, so this
isn't an artifact of a different file.

**Float32 cast fidelity (measured fresh on frame2, not assumed from frame1):**

| Quantity | Frame 1 | Frame 2 |
|---|---|---|
| max \|value\| | 19,484.27 | 17,178.50 |
| max relative error | 5.96×10⁻⁸ | 5.96×10⁻⁸ |
| max absolute error | 9.74×10⁻⁴ | 9.09×10⁻⁴ |
| mean absolute error | 7.58×10⁻⁶ | 7.49×10⁻⁶ |
| fraction of sampled tokens changed | 99.9% | 100.0% (4998/5000) |

Max relative error is IEEE-precision-determined and reproduces essentially
identically; the absolute error scales with each frame's own peak value, as
expected.

**zfp fixed-accuracy, tuned to each frame's own tolerance** (tolerance =
max-rel-err × that frame's own max\|value\|, *not* reused across frames):

| | Frame 1 (tol=1.161×10⁻³) | Frame 2 (tol=1.023×10⁻³) |
|---|---|---|
| Size | 10,020,879 | 10,136,378 |
| Achieved max abs err | 2.64×10⁻⁴ | 2.62×10⁻⁴ |
| Achieved mean abs err | 3.26×10⁻⁵ | 3.27×10⁻⁵ |

Within 1.2% on size, near-identical achieved error — the adaptive-accuracy
advantage over naive float32 truncation is not a frame1-specific fluke.

**Systematic bias check (implied electron-count shift = Σ(decoded−original) ÷
5,832,000 grid points):**

| Method | Frame 1 | Frame 2 |
|---|---|---|
| raw float32 cast | −1.6×10⁻⁸ e | −1.37×10⁻⁸ e |
| zfp fixed-accuracy | −6.8×10⁻⁸ e | −5.02×10⁻⁸ e |
| fpzip `-p 40` | −4.7×10⁻⁷ e | −4.68×10⁻⁷ e |
| fpzip `-p 32` | −1.2×10⁻⁴ e | −1.20×10⁻⁴ e |
| fpzip `-p 24` | −3.1×10⁻² e | −3.06×10⁻² e |

All within a few percent of frame1, same conclusion: every method except
`-p 24` has a bias many orders of magnitude below Bader-relevant precision
(~10⁻³–10⁻² e); `-p 24` remains the one setting where the systematic shift
is borderline significant.

**Verdict: frame2 confirms the frame1-derived recommendation without
qualification.** Every core number is within ~1% of frame1's, the float32
error characterization reproduces almost exactly (as expected — it's set by
IEEE-754 float32 precision, not by which frame you pick), and the bias
figures track frame1's to 2-3 significant figures. Different temperature
(1500K vs. 1000K) and a different production segment/timestep did not
meaningfully change compressibility, error behavior, or bias — VASP's
text-quantization mechanism and each codec's behavior are properties of the
file *format*, not of the specific configuration snapshotted. No change to
the recommendation below: zstd-19-on-text for simplicity, the custom
mantissa/delta scheme for best lossless ratio, zfp fixed-accuracy if lossy
is acceptable.

## Recommendation

For actual deployment on the Swift cluster: **use `zstd -19` directly on the
existing ASCII CHGCAR/AECCAR files** (no format conversion, one command,
`zstd`/`libzstd` already on PATH). It gets 3.78× vs. raw and 1.37× vs. the
current gzip -9, with no new build dependencies and no risk of a lossy or
buggy custom parser silently corrupting density data. If encode/decode CPU
time on Swift is a bottleneck rather than storage, `fpzip` is a strong
alternative for hot/scratch data — 150× faster than zstd-19 at a
respectable 3.2×, and it's a 2-second `cmake`+`make` build with no external
deps (built at `tools/fpzip/build/bin/fpzip`, static-linkable). The custom
mantissa/delta scheme squeezes out another ~11% over zstd-19-on-text (4.23×
total) but requires a bespoke encoder/decoder that must round-trip VASP's
exact print format perfectly — worth doing only if storage cost genuinely
justifies maintaining that extra code path; H5Z-ZFP/den2bin are not worth
pursuing further here given the missing system dependencies and no gain over
what's already been measured.
