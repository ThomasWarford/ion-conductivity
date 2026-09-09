"""Read MP amorphous_diffusivity (NCSD) trajectory files and assemble the
continuous equilibrium-volume production trajectory.

One .jsonl.gz file == one composition. Each JSON line is one MD segment
(serialized pymatgen Trajectory). Segments, ordered by their Mongo ObjectId
timestamp, are: 3 EOS-scan volumes, then a variable number of contiguous runs
at the equilibrium volume (iterative 4 ps equilibration + 10 ps production
chunks). The equilibrium-volume segments join frame-to-frame with no gap.

Timestep is 2 fs, every step written (see project memory).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np

try:
    import orjson

    def _loads(s):
        return orjson.loads(s)
except ImportError:  # pragma: no cover
    import json

    def _loads(s):
        return json.loads(s)

FS_PER_FRAME = 2.0
PS_PER_FRAME = FS_PER_FRAME / 1000.0


@dataclass
class Traj:
    system: str            # chemical_system tag, e.g. "Ag-Li"
    temperature: float     # K
    contribs_id: str
    symbols: np.ndarray    # (natoms,) element strings
    frac: np.ndarray       # (nframes, natoms, 3) fractional coords, wrapped
    lattice: np.ndarray    # (3, 3) Angstrom, constant
    energy: np.ndarray     # (nframes,) e_0_energy per frame, eV (NaN if absent)
    n_equil_frames: int = 0  # leading frames from <4 ps equilibration segments
    ps_per_frame: float = PS_PER_FRAME

    @property
    def nframes(self) -> int:
        return self.frac.shape[0]

    @property
    def natoms(self) -> int:
        return self.frac.shape[1]

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.lattice)))

    @property
    def total_ps(self) -> float:
        return (self.nframes - 1) * self.ps_per_frame

    def cart_unwrapped(self) -> np.ndarray:
        """(nframes, natoms, 3) cartesian displacement trajectory, PBC-unwrapped,
        referenced so frame 0 sits at the wrapped position."""
        df = np.diff(self.frac, axis=0)
        df -= np.round(df)
        disp = np.zeros_like(self.frac)
        np.cumsum(df, axis=0, out=disp[1:])
        disp += self.frac[0]
        return disp @ self.lattice


def _seg_records(path):
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line.strip():
                yield _loads(line)


def _decode_segment(rec):
    tr = rec["trajectory"]
    oid = rec["metadata"]["_id"]
    ts = int(oid[:8], 16)
    lat = np.asarray(tr["lattice"], dtype=float)
    lat = lat[0] if lat.ndim == 3 else lat
    frac = np.asarray(tr["coords"], dtype=float)
    if tr.get("coords_are_displacement"):
        raise ValueError("unexpected displacement coords")
    syms = np.array([s["element"] if isinstance(s, dict) else s
                     for s in tr["species"]])
    fp = tr.get("frame_properties") or []
    en = np.array([f.get("e_0_energy", f.get("total", np.nan)) for f in fp], dtype=float)
    if len(en) < frac.shape[0]:
        en = np.append(en, np.full(frac.shape[0] - len(en), np.nan))
    return dict(ts=ts, lat=lat, frac=frac, syms=syms, energy=en[:frac.shape[0]],
                temperature=float(rec["metadata"]["temperature"]),
                contribs_id=rec["metadata"]["contribs_id"])


def load_production(path, vol_rtol=2e-3, join_atol=1e-3) -> Traj:
    """Assemble the contiguous equilibrium-volume trajectory from a file.

    Returns every equilibrium-volume frame (iterative equilibration + production)
    as one continuous trajectory; strip a burn-in downstream.
    """
    segs = sorted((_decode_segment(r) for r in _seg_records(path)), key=lambda s: s["ts"])
    if not segs:
        raise ValueError(f"no segments in {path}")

    eq_vol = float(abs(np.linalg.det(segs[-1]["lat"])))
    eq = [s for s in segs
          if abs(abs(np.linalg.det(s["lat"])) - eq_vol) <= vol_rtol * eq_vol]

    # keep only the maximal contiguous (frame-joining) tail ending at the last seg
    chain = [eq[-1]]
    for s in reversed(eq[:-1]):
        nxt = chain[0]
        d = s["frac"][-1] - nxt["frac"][0]
        d -= np.round(d)
        if np.sqrt(np.mean((d @ s["lat"]) ** 2)) <= join_atol:
            chain.insert(0, s)
        else:
            break

    parts = [chain[0]["frac"]] + [s["frac"][1:] for s in chain[1:]]
    frac = np.concatenate(parts, axis=0)
    frac -= np.floor(frac)
    en_parts = [chain[0]["energy"]] + [s["energy"][1:] for s in chain[1:]]
    energy = np.concatenate(en_parts, axis=0)

    # frames belonging to leading short (<= ~4 ps == 2000-frame) runs: treat as
    # iterative equilibration, the rest as production (matches the MP pipeline).
    n_equil = 0
    for s in chain:
        if s["frac"].shape[0] <= 2500:
            n_equil += s["frac"].shape[0] - (1 if n_equil else 0)
        else:
            break

    import os
    system = os.path.basename(os.path.dirname(path)).split("chemical_system=")[-1]
    return Traj(system=system, temperature=segs[0]["temperature"],
                contribs_id=segs[0]["contribs_id"], symbols=chain[0]["syms"],
                frac=frac, lattice=chain[0]["lat"], energy=energy,
                n_equil_frames=n_equil)


def segment_summary(path) -> list[dict]:
    out = []
    for s in sorted((_decode_segment(r) for r in _seg_records(path)), key=lambda s: s["ts"]):
        out.append(dict(ts=s["ts"], nframes=s["frac"].shape[0],
                        volume=float(abs(np.linalg.det(s["lat"]))),
                        natoms=len(s["syms"])))
    return out
