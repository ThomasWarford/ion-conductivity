#!/usr/bin/env python3
"""Simple 3D MED-style spatial predictor (median of left/up/behind neighbors)
on a CHGCAR float64 grid; writes the float64 residual for external zstd test."""
import sys
import numpy as np

def main():
    infile, nx, ny, nz, outfile = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    arr = np.fromfile(infile, dtype=np.float64).reshape(nz, ny, nx)

    left = np.pad(arr, ((0,0),(0,0),(1,0)), mode='edge')[:, :, :-1]
    up = np.pad(arr, ((0,0),(1,0),(0,0)), mode='edge')[:, :-1, :]
    behind = np.pad(arr, ((1,0),(0,0),(0,0)), mode='edge')[:-1, :, :]

    stacked = np.stack([left, up, behind], axis=0)
    pred = np.median(stacked, axis=0)

    resid = arr - pred
    resid.tofile(outfile)
    print(f"residual std={resid.std():.4f} (orig std={arr.std():.4f}), max abs={np.abs(resid).max():.4f}")

if __name__ == '__main__':
    main()
