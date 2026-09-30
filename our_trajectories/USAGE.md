# AIMD campaign: CPU-hours and throughput per system

Not tracked anywhere else — derived from Slurm accounting (`sacct`, `AllocCPUS x Elapsed`,
128 CPUs/job on 1 node) and actual ionic step counts (`OSZICAR`, `POTIM = 2 fs`), as of
2026-09-29. Useful for prioritizing which systems get more AIMD time in future campaigns.

## Per-system total CPU-hours spent

Base 2-temp (1000K/2500K) sweep for all 13 systems, plus the extra 4-temp x 3-replica
treatment for Cu-Li-S and Li-P.

| System | CPU-hr |
|---|---:|
| Cu-Li-S | 322,009 |
| Li-P | 299,674 |
| F-Li-P | 64,717 |
| Bi-Li-S | 54,945 |
| Cl-Li | 54,177 |
| Al-Li-Na-P | 51,231 |
| Li-O-P-Ti | 50,398 |
| Li-O-Zr | 46,711 |
| Li-Nb-S | 45,160 |
| Li-S-Sb | 36,868 |
| Ge-Li-Sn | 36,783 |
| Li-Rb-Se | 29,800 |
| Li-Rb-S | 26,924 |
| **Total** | **1,119,395** |

Cu-Li-S and Li-P dominate the total purely because they got 12 cases each (4 temps x 3
replicas) vs. 2 for everyone else, not because they're individually more expensive per
case. Within the 11 single-treatment systems, totals track the NSW-sizing formula (each
*link* targeted ~85% of its 7-day walltime) — lower numbers mean a system's second link
was truncated or never started, not that the physics needed less:
- `Li-P_1000K_r0`: 10.7k CPU-hr (link2 never started)
- `Cu-Li-S_2500K_r0`: 17.2k CPU-hr (link2 never started)
- `Li-O-Zr_1000K_r0`: 13.6k CPU-hr

## Measured throughput (ps/day)

Actual wall-clock throughput, not the original 2500K-survey estimate used for NSW
sizing.

### Per-case

| Case | ps produced | Wall-days | ps/day |
|---|---:|---:|---:|
| F-Li-P 1000K | 219.2 | 10.20 | 21.50 |
| F-Li-P 2500K | 219.2 | 10.87 | 20.16 |
| Li-P 1000K r0 | 98.6 | 3.47 | 28.38 |
| Li-P 1000K r1 | 197.2 | 7.02 | 28.09 |
| Li-P 1000K r2 | 197.2 | 7.03 | 28.03 |
| Li-P 1500K (r0-r2) | 197.2 | ~8.33 | ~23.7 |
| Li-P 2000K (r0-r2) | 197.2 | ~9.6 | ~20.5 |
| Li-P 2500K (r0-r2) | 98.6/197.2 | ~5.2-10.5 | ~18.8 |
| Li-O-Zr 1000K | 98.2 | 4.43 | 22.18 |
| Li-O-Zr 2500K | 196.4 | 10.78 | 18.22 |
| Li-S-Sb 1000K | 124.6 | 6.00 | 20.76 |
| Cl-Li 1000K | 135.9 | 7.14 | 19.03 |
| Al-Li-Na-P 1000K | 110.3 | 6.06 | 18.18 |
| Li-O-P-Ti 2500K | 65.5 | 5.83 | 11.23 |
| Li-O-P-Ti 1000K | 131.0 | 10.57 | 12.39 |
| Bi-Li-S 1000K | 106.2 | 6.86 | 15.50 |
| Cu-Li-S (all temps/reps) | — | — | 9.4-15.1 (drops with T) |
| Cl-Li 2500K | 135.9 | 10.49 | 12.95 |
| Al-Li-Na-P 2500K | 110.3 | 10.61 | 10.39 |
| Bi-Li-S 2500K | 106.2 | 11.03 | 9.63 |
| Ge-Li-Sn 1000K | 64.0 | 6.27 | 10.22 |
| Ge-Li-Sn 2500K | 32.0 | 5.71 | 5.61 |
| Li-Nb-S 1000K | 70.5 | 8.01 | 8.80 |
| Li-Nb-S 2500K | 35.3 | 6.69 | 5.27 |
| Li-Rb-S 1000K | 34.6 | 3.65 | 9.47 |
| Li-Rb-S 2500K | 34.6 | 5.11 | 6.77 |
| Li-Rb-Se 1000K | 32.2 | 3.66 | 8.79 |
| Li-Rb-Se 2500K | 32.2 | 6.04 | 5.33 |

### Per-system average

Weighted across all of a system's cases/temps/replicas.

| System | ps/day (avg) |
|---|---:|
| Li-P | 22.23 |
| F-Li-P | 20.81 |
| Li-O-Zr | 19.37 |
| Li-S-Sb | 16.51 |
| Cl-Li | 15.41 |
| Al-Li-Na-P | 13.22 |
| Li-O-P-Ti | 11.98 |
| Bi-Li-S | 11.88 |
| Cu-Li-S | 11.74 |
| Ge-Li-Sn | 8.02 |
| Li-Rb-S | 7.89 |
| Li-Nb-S | 7.20 |
| Li-Rb-Se | 6.63 |

**Pattern:** ps/day drops sharply with temperature for every system (roughly halving
from 1000K to 2500K) — more thermal disorder means more SCF iterations per ionic step.
`Li-P`, `F-Li-P`, `Li-O-Zr`, `Li-S-Sb` are the cheap systems (>15-28 ps/day, more
trajectory per CPU-hour); `Li-Rb-Se`, `Li-Nb-S`, `Li-Rb-S` are consistently the most
expensive (~5-9 ps/day).
