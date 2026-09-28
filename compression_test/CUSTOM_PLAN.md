# Custom approach — plan

Observed so far: parsing ASCII -> raw float64 binary and compressing that with
zstd/xz loses to compressing the ASCII text directly (zstd-19 text: 28.09MB vs
xz-9e binary: 34.16MB, from a 106.18MB CHGCAR). Byte-shuffling the float64
array barely helps (36.3MB). Reason: VASP prints exactly 11 significant
decimal digits (`0.dddddddddddE+dd`). That is a lossy *quantization* of the
true double, and it manufactures long repeated digit-runs across nearby grid
points (smooth density field) that LZ/entropy coders exploit well. A raw
float64 retains the full 53-bit mantissa, most of which is uncorrelated
"noise" relative to the 11-digit precision actually present, so it compresses
worse than the quantized text, and per-byte shuffling doesn't fix that because
the noise is in the low mantissa bits, not byte alignment.

**Plan:** decompose each value into its exact fixed-point representation
implied by the text (11-digit integer mantissa + 2-digit exponent + sign),
store mantissas and exponents as separate integer arrays (removes ASCII
formatting overhead: spaces, `.`, `E`, sign chars — ~45% of the text is
non-digit decoration), delta-encode the mantissa array along the fastest grid
axis (x, since CHGCAR is smooth and stored x-fastest) to shrink the magnitude
of the common case, then zstd -19 each stream. This keeps the exact precision
already present (still lossless vs. the text) while dropping formatting
overhead and giving zstd smaller/more repetitive symbols than either raw text
or raw float64 mantissas.
