#!/usr/bin/env python3
"""Custom lossless encoder for VASP CHGCAR-style grids.
Decompose text values (0.dddddddddddE+dd, 11-digit mantissa) into
sign+mantissa (int64) and exponent (int8) streams, delta-encode mantissa
along fastest axis, write both streams for external compression (zstd/xz).
"""
import sys, re
import numpy as np

def main():
    infile, nx, ny, nz, outprefix = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    n = nx * ny * nz
    with open(infile) as f:
        text = f.read()
    # tokens like -0.12345678901E+02 or 0.12345678901E-05
    toks = text.split()
    assert len(toks) == n, f"{len(toks)} vs {n}"

    mant = np.empty(n, dtype=np.int64)
    exp = np.empty(n, dtype=np.int8)
    pat = re.compile(r'^(-?)0?\.(\d{11})E([+-]\d{2})$')
    for i, t in enumerate(toks):
        m = pat.match(t)
        sign, digits, e = m.groups()
        v = int(digits)
        if sign == '-':
            v = -v
        mant[i] = v
        exp[i] = int(e)

    mant_delta = mant.copy()
    mant_delta[1:] -= mant[:-1]  # simple delta along flattened (x-fastest) order

    mant_delta.tofile(outprefix + '.mant_delta.i64')
    exp.tofile(outprefix + '.exp.i8')
    print(f"mant_delta bytes={mant_delta.nbytes} exp bytes={exp.nbytes} total={mant_delta.nbytes+exp.nbytes}")

if __name__ == '__main__':
    main()
