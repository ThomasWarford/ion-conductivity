#!/usr/bin/env python3
"""Byte-plane split + within-plane delta filter (Aras Pranckevicius technique).
encode: arr (float32/float64) -> byte-split (all byte0s, then byte1s, ...) ->
delta each plane mod 256.
decode: inverse. Verifies round trip bit-exact.
"""
import sys
import numpy as np

def encode(arr_bytes_2d):
    # arr_bytes_2d: shape (N, nbytes) uint8, one row per value
    planes = arr_bytes_2d.T.copy()  # (nbytes, N)
    deltas = np.empty_like(planes)
    deltas[:, 0] = planes[:, 0]
    deltas[:, 1:] = planes[:, 1:] - planes[:, :-1]  # wraps mod 256 automatically (uint8)
    return deltas

def decode(deltas):
    planes = np.cumsum(deltas.astype(np.uint16), axis=1).astype(np.uint8)
    # cumsum on uint16 then cast back to uint8 reproduces mod-256 wraparound cumsum
    return planes.T.copy()  # (N, nbytes)

def main():
    infile, dtype_str, outfile = sys.argv[1], sys.argv[2], sys.argv[3]
    dtype = np.float64 if dtype_str == 'double' else np.float32
    nbytes = 8 if dtype_str == 'double' else 4
    arr = np.fromfile(infile, dtype=dtype)
    raw = arr.view(np.uint8).reshape(-1, nbytes)
    enc = encode(raw)
    enc.tofile(outfile)

    # round-trip verify
    dec = decode(enc)
    assert np.array_equal(dec, raw), "round trip MISMATCH"
    arr2 = dec.reshape(-1).view(dtype)
    assert np.array_equal(arr, arr2), "value round trip MISMATCH"
    print(f"OK: {infile} -> {outfile}, bytes={enc.nbytes}, round-trip verified bit-exact")

if __name__ == '__main__':
    main()
