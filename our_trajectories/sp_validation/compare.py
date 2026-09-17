"""Compare our single-point VASP energies/forces against
  (a) the NCSD reference frame in each run dir's ref.json, and
  (b) our own fully-converged 'gold' tier.

    conda run -n dft python compare.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Outcar, Vasprun

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
TIERS = ("exact", "tight", "gold")


def read_ours(d: Path) -> dict | None:
    if (d / "vasprun.xml").exists():
        try:
            vr = Vasprun(d / "vasprun.xml", parse_potcar_file=False,
                         exception_on_bad_xml=False)
            step = vr.ionic_steps[-1]
            return dict(e0=vr.final_energy, e_fr=step["e_fr_energy"],
                        e_wo=step["e_wo_entrp"], forces=np.array(step["forces"]),
                        nscf=len(step["electronic_steps"]),
                        conv=vr.converged_electronic)
        except Exception:
            pass
    if (d / "OUTCAR").exists():
        oc = Outcar(d / "OUTCAR")
        f = np.array(oc.read_table_pattern(
            r"TOTAL-FORCE \(eV/Angst\)\n\s*-+\n",
            r"\s+".join([r"([+-]?\d+\.\d+)"] * 6),
            r"\s*-+\n\s+total drift", last_one_only=True), dtype=float)
        return dict(e0=oc.final_energy, e_fr=np.nan, e_wo=np.nan,
                    forces=f[:, 3:] if f.size else np.array([]),
                    nscf=-1, conv=None)
    return None


def fstats(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    d = a - b
    return (float(np.sqrt(np.mean(d ** 2))), float(np.mean(np.abs(d))),
            float(np.abs(d).max()))


def main() -> None:
    cases = (RUNS / "cases.txt").read_text().split()
    rows = []
    for case in cases:
        ref = json.loads((RUNS / case / "exact" / "ref.json").read_text())
        n = ref["natoms"]
        ref_f = np.array(ref["forces"])
        got = {t: read_ours(RUNS / case / t) for t in TIERS}
        gold = got["gold"]
        if gold is not None and (gold["forces"].size == 0 or gold["forces"].shape != ref_f.shape):
            gold = None  # interrupted (e.g. TIMEOUT) mid-SCF, no forces yet
        for t in TIERS:
            g = got[t]
            if g is None or g["forces"].shape != ref_f.shape:
                rows.append((case, t, n, None))
                continue
            # vs NCSD reference
            dE_ref = (g["e0"] - ref["e_0_energy"]) * 1000
            fr_ref = fstats(g["forces"], ref_f)
            # vs our gold
            if gold is not None and t != "gold":
                dE_gold = (g["e0"] - gold["e0"]) * 1000
                fr_gold = fstats(g["forces"], gold["forces"])
            else:
                dE_gold, fr_gold = None, (None, None, None)
            rows.append((case, t, n, dict(dE_ref=dE_ref, fr_ref=fr_ref,
                                          dE_gold=dE_gold, fr_gold=fr_gold,
                                          nscf=g["nscf"], conv=g["conv"])))

    h = (f"{'case':16s} {'tier':6s} {'N':>4s} | {'dE_ref':>9s} {'meV/at':>8s} "
         f"{'Frmse':>8s} {'Fmax':>8s} | {'dE_gold':>9s} {'meV/at':>8s} "
         f"{'Frmse':>8s} {'Fmax':>8s} | {'nscf':>4s}")
    print(h)
    print("-" * len(h))
    for case, t, n, r in rows:
        if r is None:
            print(f"{case:16s} {t:6s} {n:4d} |  (no output yet)")
            continue
        eref, mref = r["dE_ref"], r["dE_ref"] / n
        fr = r["fr_ref"]
        if r["dE_gold"] is None:
            gstr = f"{'--':>9s} {'--':>8s} {'--':>8s} {'--':>8s}"
        else:
            gstr = (f"{r['dE_gold']:9.2f} {r['dE_gold']/n:8.3f} "
                    f"{r['fr_gold'][0]:8.4f} {r['fr_gold'][2]:8.4f}")
        print(f"{case:16s} {t:6s} {n:4d} | {eref:9.2f} {mref:8.3f} "
              f"{fr[0]:8.4f} {fr[2]:8.4f} | {gstr} | {r['nscf']:4d}")

    print("\ndE_ref  : our energy(sigma->0)  -  NCSD reference e_0_energy      [meV]")
    print("dE_gold : our energy(sigma->0)  -  our fully-converged 'gold'      [meV]")
    print("Frmse/Fmax over all 3N force components                           [eV/Ang]")
    print("\ntiers:  exact = MPMorph INCAR verbatim (loose EDIFF, LREAL=T, PREC=Normal)")
    print("        tight = same basis/grid, EDIFF=1e-7")
    print("        gold  = EDIFF=1e-7 + LREAL=F + PREC=Accurate + ENCUT=700")
    print("gold-vs-tight = the systematic error baked into MP's own settings.")


if __name__ == "__main__":
    main()
