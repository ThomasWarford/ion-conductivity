#!/usr/bin/env python3
"""Parse a VASP CHGCAR/AECCAR main density grid (ASCII) into raw binary float arrays.
Usage: chgcar_to_bin.py <infile> <outprefix>
Writes <outprefix>.f64 and <outprefix>.f32, plus prints stats.
Also writes <outprefix>.header (everything before the grid, verbatim) and
<outprefix>.tail (augmentation occupancies + anything after, verbatim) so the
file can be losslessly reconstructed later if desired.
"""
import sys
import numpy as np

def main():
    infile, outprefix = sys.argv[1], sys.argv[2]
    with open(infile, 'r') as f:
        lines = f.readlines()

    # find grid dims line: first line after "Direct"/"Cartesian" that is exactly 3 ints
    dims_idx = None
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            dims_idx = i
            nx, ny, nz = map(int, parts)
            break
    if dims_idx is None:
        raise SystemExit("could not find grid dims line")

    n = nx * ny * nz
    import math
    nlines = math.ceil(n / 5)
    data_start = dims_idx + 1
    data_end = data_start + nlines  # exclusive

    header_text = ''.join(lines[:data_start])
    data_text = ''.join(lines[data_start:data_end])
    tail_text = ''.join(lines[data_end:])

    arr = np.fromstring(data_text, dtype=np.float64, sep=' ')
    assert arr.size == n, f"parsed {arr.size} values, expected {n}"

    arr.tofile(outprefix + '.f64')
    arr32 = arr.astype(np.float32)
    arr32.tofile(outprefix + '.f32')

    with open(outprefix + '.header', 'w') as f:
        f.write(header_text)
    with open(outprefix + '.tail', 'w') as f:
        f.write(tail_text)

    # round-trip check: does float32 preserve the printed text exactly?
    # VASP prints with 11 sig-fig mantissa e.g. 0.42892911850E+02
    sample = arr[:5]
    print(f"nx,ny,nz = {nx},{ny},{nz}  n={n}")
    print(f"data_text bytes = {len(data_text)}")
    print(f"f64 bytes = {arr.nbytes}, f32 bytes = {arr32.nbytes}")
    # check max relative error introduced by float32
    rel_err = np.max(np.abs((arr32.astype(np.float64) - arr) / np.where(arr != 0, arr, 1)))
    print(f"max relative error f64->f32: {rel_err:.3e}")

if __name__ == '__main__':
    main()
