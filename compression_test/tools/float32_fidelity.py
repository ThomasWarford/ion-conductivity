#!/usr/bin/env python3
"""Check how much precision float32 casting throws away vs VASP's printed text."""
import numpy as np, math

def vasp_fmt(v):
    # VASP prints 0.dddddddddddE+dd (11-digit mantissa, sign, 2-digit exponent)
    s = f"{v:.11E}"  # e.g. -4.28929118500E+01 or 4.28929118500E+01
    sign = '-' if s.startswith('-') else ''
    s = s.lstrip('-')
    mant, exp = s.split('E')
    d0, rest = mant.split('.')
    digits = (d0 + rest)  # combine leading digit + fractional digits -> 12 digits total (d0 + 11)
    # renormalize to VASP's 0.xxxxxxxxxxx form: shift exponent by +1, digits[0] is the original integer part
    e = int(exp) + 1
    return f"{sign}0.{digits[:11]}E{e:+03d}"

lines = open('samples/frame1/CHGCAR').readlines()
dims_idx = None
for i, l in enumerate(lines):
    p = l.split()
    if len(p) == 3 and all(x.isdigit() for x in p):
        dims_idx = i
        nx, ny, nz = map(int, p)
        break
n = nx * ny * nz
nlines = math.ceil(n / 5)
data_start = dims_idx + 1
data_end = data_start + nlines
orig_toks = ''.join(lines[data_start:data_end]).split()
assert len(orig_toks) == n

f64 = np.fromfile('samples/frame1/CHGCAR.f64', dtype=np.float64)
f32 = f64.astype(np.float32)
back = f32.astype(np.float64)

rel_err = np.abs((back - f64) / np.where(f64 != 0, f64, 1))
print(f"max relative error: {rel_err.max():.3e}")
print(f"mean relative error: {rel_err.mean():.3e}")

rng = np.random.default_rng(0)
idxs = rng.choice(n, size=5000, replace=False)
changed = 0
for i in idxs:
    orig = orig_toks[i].replace('-.', '-0.')
    if orig.startswith('.'):
        orig = '0' + orig
    recon = vasp_fmt(back[i])
    if recon != orig:
        changed += 1
print(f"fraction of sampled tokens changed by f64->f32->f64 round trip: {changed}/{len(idxs)} = {changed/len(idxs)*100:.1f}%")
