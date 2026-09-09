# NCSD amorphous-diffusivity screen

Screens the MP `amorphous_diffusivity` multi-temperature trajectories
(`../data_ncsd`, 219 compositions x {1000,1500,2000,2500} K) for

1. **how fast the Li self-diffusivity converges** within the AIMD window, and
2. **how large it is**,

to find systems where an MLIP-vs-AIMD validation of diffusion coefficients /
cross-correlations is feasible.

## Get the data

```
conda run -n dft python fetch_ncsd.py --dry-run   # 48 files, 4.0 GB
conda run -n dft python fetch_ncsd.py            # resumable; re-run to verify
ln -s $SCRATCH/data/ncsd ../data_ncsd            # if the symlink is missing
```

`fetch_ncsd.py` reads the public S3 bucket over plain REST (standard library
only -- no boto3, no credentials) and mirrors its layout under `--dest`
(default `$SCRATCH/data/ncsd`):

```
<dest>/trajectories/temperature=<T>K/chemical_system=<S>/dt=<...>.jsonl.gz
```

A file whose size already matches the bucket listing is skipped, so re-running
costs one listing request per (system, temperature) and nothing else. Downloads
land as `<name>.part` and are renamed only when complete.

By default it fetches the 17 systems in `md_targets.txt`, the set chosen for
long MLIP MD + AIMD validation of Li diffusion. **Every one is charge-balanced**
(see the section below); the file itself carries the reasoning per system.

| role | systems |
|---|---|
| trend series (paired, run both or neither) | `Li-N`/`Li-P`, `Al-Li-Na-P`/`Al-K-Li-P`, `Li-Rb-S`/`Li-Rb-Se`, `F-Li-P`/`As-F-Li`, `Li-O-P-W`/`Li-O-P-Ti`, `Li-S-Sb`/`Bi-Li-S` |
| standalone | `Cl-Li` (molten-salt benchmark), `Cu-Li-S` (two mobile cations), `Li-O-Zr` (the converging oxide), `Li-Nb-S` (fastest balanced conductor) |
| control, not an electrolyte | `Ge-Li-Sn` (Zintl phase, best-converging system in the set) |

Each trend series keeps every element but one and swaps that one for a
same-group or same-oxidation-state congener; five of the six do it at identical
stoichiometry, so the substitution is the only thing that changes. Enumerate the
alternatives with `python group_series.py` (add `--isovalent`).

**There is no `Li-O` binary in the dataset.** The Li-X binaries on offer are
Li-N, Li-P, Li-S, Li-Se, Li-Te, Br-Li, Cl-Li, I-Li and the metals -- but only
Li-N, Li-P and Li-Te are charge-balanced, so `Li-O-Zr` stands in for the Li
oxide. Asking for a system that is not in the bucket is a hard error that prints
the full available list.

Other options: `--systems Li-S,Li-Se`, `--systems all` (the whole dataset, 74 GB
without 5000 K), `--temps`, `--procs`, `--dest`. **5000 K** works but is not in
the default temperature list -- those trajectories are roughly 5x larger than
the rest (Li-S alone is 633 MB), so the 12 targets at 5000 K are about 7 GB.

## Pipeline

| file | role |
|---|---|
| `fetch_ncsd.py`  | download the trajectories from the public S3 bucket (see above) |
| `md_targets.txt` | the 12 systems chosen for long MD, one per line |
| `ncsd_traj.py`   | read a `.jsonl.gz` file, assemble the contiguous equilibrium-volume trajectory, drop the leading <=4 ps equilibration segments (2 fs timestep) |
| `diffusivity.py` | FFT MSD, 3-D Einstein D (`MSD -> 6 D t`), diffusive exponent beta, block/half/lag-window/prefix convergence, collective ("charge") MSD |
| `worker.py`      | one file -> rows (per species + `*all*`) |
| `run_screen.py`  | parallel driver -> `results/screen.{jsonl,csv}` (resumable) |
| `screen.sbatch`  | runs the above on 1 Perlmutter CPU node (~2 min) |
| `fetch_mp_reference.py` | MP published `diffusivity` + MSD_dt tables -> `results/mp_reference.csv` |
| `charge_balance.py` | classify a cell composition as a stoichiometric ionic compound (`balanced`), a composition no oxidation-state assignment neutralises (`unbalanced`), or an alloy with no anion (`no_anion`) |
| `group_series.py` | chemical-trend series: keep all but one element, swap it for a same-group congener; flags which series are fully charge-balanced and which are blocked by a non-compound congener |
| `report.py`      | -> `results/summary.md`, `results/shortlist_*.csv`, `results/figs/*.png`; `--charge-balanced-only` writes the `*_cb` variants |
| `onsager.py`     | Onsager coefficients per species: `D_self`, `D_coll`, correlation factor `f` -> `results/onsager.csv` |
| `onsager_plots.py` / `onsager_raw.py` / `onsager_species.py` | figures for one ion, the un-averaged per-origin clouds, and the cross-ion comparison |

```
sbatch screen.sbatch
python fetch_mp_reference.py      # optional cross-check
python report.py                            # all 219 compositions
python report.py --charge-balanced-only     # the 60 that are actual compounds
python onsager_plots.py --charge-balanced-only
```

## Key outputs

- `results/summary_cb.md` – narrative + tables, charge-balanced compositions only
  (**start here**); `results/summary.md` is the same over all 219
- `results/shortlist_all_cb.csv` / `results/shortlist_electrolyte_cb.csv` – systems
  whose Li channel converges at >=2 temperatures (all chemistries / insulating-anion
  only); un-suffixed versions are the unfiltered equivalents
- `results/figs_cb/`, `results/figs/li_cb/` – the same figures, filtered
- `results/screen.csv` – every (system, T, species) with all metrics

## Charge balance: read this before using any shortlist

The NCSD cells are PACKMOL packings at a target composition, and **most of those
compositions are not compounds**. Of the 219 systems: 60 balanced, 101 unbalanced,
58 alloys with no anion. `Br-Li` is `Br25Li75` = Li3Br (+0.5 e/atom), not LiBr;
`Li-S` is Li3S, not Li2S. Those cells are not charged -- VASP runs them neutral --
they are Li-rich metallic melts, so their large "Li diffusivity" is metallic
self-diffusion, not ionic conduction, and they are the systems that dominated the
top of the unfiltered shortlist. Always run the `--charge-balanced-only` variants
when picking solid-electrolyte validation targets.

`md_targets.txt` used to list seven such cells (`Br-Li`, `Li-S`, `Li-Se`,
`Li-O-Si`, `Li-S-Sn`, `Li-Se-Sn`, `Li-P-Sn`). It has been rewritten: every system
in it now passes the balance test.

## Convention note

The MPContribs `diffusivity` field is computed as `MSD_slope / 2` (1-D Einstein
relation on the 3-D MSD) and is therefore **3x** the standard `MSD_slope / 6`
used here (confirmed: ratio 3.10, IQR 2.81-3.42, n=803). Divide MP values by 3
before comparing to a properly-computed AIMD or MLIP D.
