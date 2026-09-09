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

By default it fetches the 12 systems in `md_targets.txt`, the shortlist chosen
for long MLIP MD + AIMD validation of Li diffusion:

| group | systems |
|---|---|
| binaries | `Li-N`, `Br-Li`, `Cl-Li`, `Li-S`, `Li-Se` |
| oxides | `Ga-Li-O`, `Li-O-Si` |
| ternaries | `Cu-Li-S`, `Li-S-Sn`, `Li-Se-Sn`, `Li-O-Sn`, `Li-P-Sn` |

**There is no `Li-O` binary in the dataset.** The Li-X binaries on offer are
Li-N, Li-P, Li-S, Li-Se, Li-Te, Br-Li, Cl-Li, I-Li and the metals; `Ga-Li-O` and
`Li-O-Si` stand in for the Li oxide. Asking for a system that is not in the
bucket is a hard error that prints the full available list.

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
| `report.py`      | -> `results/summary.md`, `results/shortlist_*.csv`, `results/figs/*.png` |
| `onsager.py`     | Onsager coefficients per species: `D_self`, `D_coll`, correlation factor `f` -> `results/onsager.csv` |
| `onsager_plots.py` / `onsager_raw.py` / `onsager_species.py` | figures for one ion, the un-averaged per-origin clouds, and the cross-ion comparison |

```
sbatch screen.sbatch
python fetch_mp_reference.py      # optional cross-check
python report.py
```

## Key outputs

- `results/summary.md` – narrative + tables (start here)
- `results/shortlist_all.csv` / `results/shortlist_electrolyte.csv` – systems whose
  Li channel converges at >=2 temperatures (all chemistries / insulating-anion only)
- `results/screen.csv` – every (system, T, species) with all metrics

## Convention note

The MPContribs `diffusivity` field is computed as `MSD_slope / 2` (1-D Einstein
relation on the 3-D MSD) and is therefore **3x** the standard `MSD_slope / 6`
used here (confirmed: ratio 3.10, IQR 2.81-3.42, n=803). Divide MP values by 3
before comparing to a properly-computed AIMD or MLIP D.
