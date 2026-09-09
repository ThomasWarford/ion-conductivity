# NCSD amorphous-diffusivity screen  --  convergence & magnitude

219 compositions x 4 temperatures (1000-2500 K), AIMD, 2 fs timestep, 33 ps median production run.

D is the 3-D Einstein self-diffusivity (MSD -> 6 D t). The MPContribs `diffusivity` field fits MSD -> 2 D t and is therefore ~3x larger (confirmed below).


**Convergence flags (Li):** D >= 1e-5 cm^2/s; late-time MSD exponent beta in [0.90, 1.15]; block CV <= 0.25 over 4 independent sub-runs; D stable to <25% across lag windows; |D(2nd half)/D(1st half) - 1| <= 0.4; D from a <=20 ps trajectory prefix already within 25% of the full-run D (t_run_converge); |energy drift| <= 1 meV/atom/ps.

## 0. Charge balance of the cell compositions

The NCSD cells are PACKMOL packings at a target composition, and most of those compositions are not stoichiometric compounds. A composition is **balanced** here if some assignment of common oxidation states sums to zero (`charge_balance.py`; the most electronegative element is forced into its anion state). **no_anion** means no element present forms an anion at all -- an alloy, not an ionic compound.

| status | n systems | % |
|---|---|---|
| balanced | 60 | 27% |
| unbalanced | 101 | 46% |
| no_anion | 58 | 26% |
| **total** | **219** | |

The imbalance is not marginal. `q/atom` is the smallest net formal charge per atom any oxidation-state assignment can reach:

| |q|/atom | n unbalanced systems |
|---|---|
| 0 - 0.1 | 17 |
| 0.1 - 0.25 | 16 |
| 0.25 - 0.5 | 22 |
| 0.5 - inf | 46 |

**The fastest Li movers are exactly the compositions that are not compounds** -- Li-rich packings whose 'Li diffusivity' is self-diffusion in a lithium melt. Top 12 excluded systems by D_Li(1000 K):

| system | formula | reduced | status | q/atom | anion class | D_Li 1000K |
|---|---|---|---|---|---|---|
| K-Li | K20Li90 | K2Li9 | no_anion | - | metallic | 2.43e-04 |
| Li | Li100 | Li | no_anion | - | metallic | 2.36e-04 |
| Li-Na | Li90Na20 | Na2Li9 | no_anion | - | metallic | 2.36e-04 |
| Li-Rb | Li68Rb34 | RbLi2 | no_anion | - | metallic | 1.81e-04 |
| Li-Sr | Li80Sr20 | SrLi4 | no_anion | - | metallic | 1.77e-04 |
| Br-Li | Br25Li75 | Li3Br | unbalanced | +0.50 | halide | 1.71e-04 |
| I-Li | I34Li68 | Li2I | unbalanced | +0.33 | halide | 1.61e-04 |
| Li-Sn-Sr | Li84Sn7Sr14 | Sr2Li12Sn | unbalanced | +0.80 | metallic | 1.56e-04 |
| Cu-Li | Cu20Li80 | Li4Cu | no_anion | - | metallic | 1.40e-04 |
| Li-Zr | Li80Zr20 | Li4Zr | no_anion | - | metallic | 1.35e-04 |
| K-Li-Si | K10Li80Si10 | KLi8Si | unbalanced | +0.50 | network(B/C/Si/Ge) | 1.33e-04 |
| Cr-Li | Cr25Li75 | Li3Cr | no_anion | - | metallic | 1.32e-04 |

Everything below uses **all** compositions, balanced or not -- run with `--charge-balanced-only` for the filtered version.

## 1. By temperature

| T (K) | median prod (ps) | median D_Li (cm^2/s) | D_Li 10-90% range | frac diffusive (beta>0.9) | frac Li-converged | frac charge-converged |
|---|---|---|---|---|---|---|
| 1000 | 34 | 4.01e-05 | 6.2e-06 - 1.1e-04 | 0.64 | 0.37 | 0.28 |
| 1500 | 32 | 1.47e-04 | 3.6e-05 - 2.7e-04 | 0.74 | 0.48 | 0.26 |
| 2000 | 32 | 2.79e-04 | 1.2e-04 - 4.6e-04 | 0.77 | 0.46 | 0.26 |
| 2500 | 32 | 4.29e-04 | 1.7e-04 - 6.6e-04 | 0.75 | 0.50 | 0.25 |

## 2. Convergence vs magnitude, by anion chemistry (1000 K)

| anion class | n systems | median D_Li (cm^2/s) | frac Li-converged |
|---|---|---|---|
| chalcogenide | 18 | 5.27e-05 | 0.28 |
| halide | 18 | 2.68e-05 | 0.22 |
| metallic | 71 | 6.36e-05 | 0.59 |
| network(B/C/Si/Ge) | 54 | 4.17e-05 | 0.35 |
| oxide | 42 | 6.95e-06 | 0.05 |
| pnictide | 16 | 3.87e-05 | 0.50 |

## 3. Shortlist A - any chemistry, Li converged at >= 2 T  (121 systems)

Mostly Li metal / Li-alloy melts: large D, textbook convergence. Full table: results/shortlist_all.csv. Top 20 by D_Li(1000 K):

| system   | formula      | reduced   | anion_class        |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:---------|:-------------|:----------|:-------------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| K-Li     | K20Li90      | K2Li9     | metallic           |     90 |          4 |               0 |           5   |     1.0295 |         0.064  |       0.95  |   0.203 |     0.000243 |     0.000429 |     0.000757 |     0.000998 |
| Li       | Li100        | Li        | metallic           |    100 |          4 |               0 |          15   |     0.97   |         0.0725 |     nan     |   0.197 |     0.000236 |     0.000403 |     0.000641 |     0.000989 |
| Li-Na    | Li90Na20     | Na2Li9    | metallic           |     90 |          3 |               1 |           7.5 |     1.0215 |         0.074  |       1.765 |   0.174 |     0.000236 |     0.000478 |     0.000652 |     0.000791 |
| Li-Rb    | Li68Rb34     | RbLi2     | metallic           |     68 |          4 |               0 |          15   |     0.9845 |         0.0725 |       0.25  |   0.239 |     0.000181 |     0.000458 |     0.000612 |     0.001051 |
| Li-Sr    | Li80Sr20     | SrLi4     | metallic           |     80 |          4 |               0 |           5   |     1.0575 |         0.0955 |       0.6   |   0.197 |     0.000177 |     0.00032  |     0.000515 |     0.000713 |
| Br-Li    | Br25Li75     | Li3Br     | halide             |     75 |          3 |               2 |           7.5 |     1.072  |         0.067  |       0.39  |   0.235 |     0.000171 |     0.000366 |     0.00063  |     0.000892 |
| I-Li     | I34Li68      | Li2I      | halide             |     68 |          4 |               0 |          15   |     0.9225 |         0.2045 |       0.3   |   0.242 |     0.000161 |     0.000397 |     0.000565 |     0.000945 |
| Li-Sn-Sr | Li84Sn7Sr14  | Sr2Li12Sn | metallic           |     84 |          4 |               0 |           5   |     1.1145 |         0.0975 |       0.16  |   0.214 |     0.000156 |     0.000324 |     0.000597 |     0.000643 |
| Cu-Li    | Cu20Li80     | Li4Cu     | metallic           |     80 |          3 |               2 |           5   |     1.0145 |         0.0775 |       0.79  |   0.245 |     0.00014  |     0.000298 |     0.000445 |     0.000888 |
| Li-Zr    | Li80Zr20     | Li4Zr     | metallic           |     80 |          3 |               2 |          10   |     1.041  |         0.093  |       0.16  |   0.222 |     0.000135 |     0.000328 |     0.000575 |     0.000577 |
| K-Li-Si  | K10Li80Si10  | KLi8Si    | network(B/C/Si/Ge) |     80 |          2 |               1 |          20   |     1.032  |         0.0985 |       0.67  |   0.257 |     0.000133 |     0.000367 |     0.000502 |     0.000885 |
| Cr-Li    | Cr25Li75     | Li3Cr     | metallic           |     75 |          3 |               1 |           7.5 |     0.9275 |         0.113  |       0.29  |   0.203 |     0.000132 |     0.000259 |     0.000385 |     0.000572 |
| Li-Sc    | Li75Sc25     | Li3Sc     | metallic           |     75 |          4 |               0 |           7.5 |     1.0235 |         0.1235 |       1.52  |   0.223 |     0.000128 |     0.000265 |     0.000452 |     0.000609 |
| Ca-Li-Sn | Ca14Li84Sn7  | Li12Ca2Sn | metallic           |     84 |          3 |               1 |           7.5 |     0.9725 |         0.1055 |       0.75  |   0.234 |     0.000127 |     0.000281 |     0.000461 |     0.000668 |
| Hf-Li    | Hf20Li90     | Li9Hf2    | metallic           |     90 |          2 |               1 |          20   |     1.069  |         0.116  |       0.13  |   0.237 |     0.000127 |     0.000354 |     0.000486 |     0.000679 |
| B-Li     | B25Li75      | Li3B      | network(B/C/Si/Ge) |     75 |          3 |               1 |          20   |     0.9755 |         0.108  |       3.595 |   0.216 |     0.000126 |     0.000266 |     0.000461 |     0.000549 |
| Li-Ta    | Li80Ta20     | Li4Ta     | metallic           |     80 |          3 |               1 |          10   |     0.9525 |         0.115  |       0.54  |   0.22  |     0.000125 |     0.000313 |     0.000409 |     0.000618 |
| Hg-Li    | Hg20Li80     | Li4Hg     | metallic           |     80 |          2 |               2 |          15   |     0.9925 |         0.076  |       0.74  |   0.231 |     0.000124 |     0.000256 |     0.000346 |     0.000731 |
| Li-W     | Li80W20      | Li4W      | metallic           |     80 |          3 |               2 |          10   |     1.055  |         0.1075 |       0.365 |   0.256 |     0.000124 |     0.000261 |     0.000638 |     0.000657 |
| Li-Si-Sr | Li80Si10Sr10 | SrLi8Si   | network(B/C/Si/Ge) |     80 |          3 |               1 |          25   |     1.004  |         0.086  |       0.745 |   0.242 |     0.000117 |     0.000299 |     0.000356 |     0.00075  |

## 4. Shortlist B - insulating-anion (halide/oxide/S/Se/N/P) chemistry, Li converged at >= 2 T  (25 systems)

The relevant set for solid-electrolyte MLIP validation. Full table: results/shortlist_electrolyte.csv.

| system    | formula        | reduced    | anion_class   |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:----------|:---------------|:-----------|:--------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Br-Li     | Br25Li75       | Li3Br      | halide        |     75 |          3 |               2 |           7.5 |     1.072  |         0.067  |       0.39  |   0.235 |     0.000171 |     0.000366 |     0.00063  |     0.000892 |
| I-Li      | I34Li68        | Li2I       | halide        |     68 |          4 |               0 |          15   |     0.9225 |         0.2045 |       0.3   |   0.242 |     0.000161 |     0.000397 |     0.000565 |     0.000945 |
| Li-O-W    | Li80O15W5      | Li16WO3    | oxide         |     80 |          4 |               1 |           5   |     1.04   |         0.112  |       0.235 |   0.295 |     8.8e-05  |     0.000264 |     0.000452 |     0.000708 |
| Cl-Li     | Cl50Li50       | LiCl       | halide        |     50 |          3 |               0 |           5   |     1.007  |         0.1135 |       1.05  |   0.259 |     8.3e-05  |     0.00026  |     0.000345 |     0.000538 |
| Li-Se     | Li75Se25       | Li3Se      | chalcogenide  |     75 |          3 |               1 |           5   |     0.991  |         0.0635 |       2.81  |   0.282 |     7.7e-05  |     0.000218 |     0.000386 |     0.000556 |
| Li-Se-Sn  | Li84Se14Sn7    | Li12SnSe2  | chalcogenide  |     84 |          4 |               1 |           7.5 |     0.962  |         0.048  |       1.92  |   0.291 |     6.7e-05  |     0.000188 |     0.000396 |     0.00048  |
| Cl-Li-Si  | Cl32Li60Si8    | Li15Si2Cl8 | halide        |     60 |          2 |               3 |          15   |     0.9905 |         0.102  |       2.7   |   0.308 |     6.7e-05  |     0.000197 |     0.000426 |     0.000543 |
| Li-S      | Li75S25        | Li3S       | chalcogenide  |     75 |          3 |               1 |          10   |     1.0345 |         0.108  |       1.68  |   0.303 |     6.3e-05  |     0.000159 |     0.000315 |     0.000552 |
| Li-S-Sn   | Li81S18Sn9     | Li9SnS2    | chalcogenide  |     81 |          3 |               3 |          15   |     0.9315 |         0.089  |       1.565 |   0.303 |     5.6e-05  |     0.000193 |     0.000268 |     0.00052  |
| Bi-Li-S   | Bi15Li45S45    | Li3BiS3    | chalcogenide  |     45 |          3 |               1 |          25   |     0.9575 |         0.1155 |       0.7   |   0.294 |     5e-05    |     0.000158 |     0.000265 |     0.000397 |
| Li-P-Sn   | Li42P27Sn36    | Li14Sn12P9 | pnictide      |     42 |          2 |               0 |          25   |     0.9365 |         0.1835 |       0.675 |   0.312 |     4.8e-05  |     0.000148 |     0.000347 |     0.000383 |
| Li-N      | Li75N25        | Li3N       | pnictide      |     75 |          4 |               0 |           5   |     1.0195 |         0.0705 |       1.64  |   0.32  |     4.7e-05  |     0.000155 |     0.000279 |     0.00046  |
| Li-S-Sb   | Li50S40Sb10    | Li5SbS4    | chalcogenide  |     50 |          3 |               2 |          15   |     0.9255 |         0.172  |       0.385 |   0.347 |     4.7e-05  |     0.000145 |     0.000316 |     0.000544 |
| Cu-Li-S   | Cu17Li51S34    | Li3CuS2    | chalcogenide  |     51 |          2 |               1 |          10   |     1.047  |         0.105  |       0.465 |   0.353 |     3.6e-05  |     0.000135 |     0.000274 |     0.00042  |
| Li-S-Ti   | Li36S36Ti9     | Li4TiS4    | chalcogenide  |     36 |          2 |               1 |          20   |     1.0005 |         0.16   |       1.12  |   0.36  |     3.5e-05  |     0.000122 |     0.000304 |     0.00041  |
| Cu-Li-O-P | Cu5Li25O40P10  | Li5CuP2O8  | oxide         |     25 |          2 |               0 |          20   |     1.0405 |         0.1275 |       1.64  |   0.312 |     2.4e-05  |     8.7e-05  |     0.000218 |     0.000168 |
| Li-P      | Li75P25        | Li3P       | pnictide      |     75 |          3 |               1 |          25   |     1.0395 |         0.144  |       1.925 |   0.472 |     1.7e-05  |     5.4e-05  |     0.00025  |     0.000425 |
| Li-O-Sn   | Li40O40Sn20    | Li2SnO2    | oxide         |     40 |          2 |               2 |          25   |     0.965  |         0.2115 |       4.16  |   0.384 |     1.5e-05  |     6.5e-05  |     0.000134 |     0.000223 |
| Li-Rb-Se  | Li34Rb34Se34   | RbLiSe     | chalcogenide  |     34 |          2 |               3 |          15   |     0.9325 |         0.1245 |       0.585 |   0.474 |     1.3e-05  |     6.3e-05  |     0.000205 |     0.000338 |
| Al-K-Li-P | Al17K34Li17P34 | K2LiAlP2   | pnictide      |     17 |          2 |               2 |          15   |     0.809  |         0.24   |       1.695 |   0.501 |     9e-06    |     6.7e-05  |     0.000132 |     0.00035  |
| Li-O-Zn   | Li50O45Zn20    | Li10Zn4O9  | oxide         |     50 |          3 |               2 |          15   |     0.971  |         0.1155 |       0.915 |   0.487 |     8e-06    |     6e-05    |     0.000133 |     0.000229 |
| Ga-Li-O   | Ga10Li50O40    | Li5GaO4    | oxide         |     50 |          2 |               3 |          15   |     1.0175 |         0.085  |       1     |   0.526 |     7e-06    |     6.1e-05  |     0.00017  |     0.000278 |
| Li-O-V    | Li20O60V20     | LiVO3      | oxide         |     20 |          2 |               2 |          25   |     0.962  |         0.141  |       1.25  |   0.477 |     7e-06    |     6.6e-05  |     0.000144 |     0.000183 |
| As-Li-O   | As20Li20O60    | LiAsO3     | oxide         |     20 |          2 |               0 |          20   |     0.957  |         0.181  |       1.725 |   0.539 |     6e-06    |     8.3e-05  |     0.000126 |     0.000281 |
| Li-O-Zr   | Li42O49Zr14    | Li6Zr2O7   | oxide         |     42 |          3 |               1 |          20   |     0.9665 |         0.163  |       0.7   |   0.486 |     6e-06    |     3.7e-05  |     9.1e-05  |     0.00017  |

## 5. Cross-check vs MPContribs

- MP `diffusivity` / refit of their own MSD table as slope/6 : median **3.10** (IQR 2.81-3.42, n=803)  -> the published values are 3x the standard 3-D Einstein D.
- MP `diffusivity` / our independent trajectory D_einstein : median **3.10** (n=804)  -> our analysis reproduces their MSD; the 3x offset is the 1D-vs-3D Einstein prefactor (2 D t vs 6 D t), not a physics disagreement.

## 6. Cross-correlations / Onsager terms

The collective (charge) MSD converges far worse than the self MSD: only 25% of systems reach a diffusive collective regime even at 2500 K (vs 50% for self-D). Systems with `charge_conv_T >= 2` in the shortlists are the ones where validating the cross terms against AIMD is also realistic.

## 7. Top picks - any chemistry

| system   | formula     | reduced   | anion_class        |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:---------|:------------|:----------|:-------------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Co-Li-Sn | Co7Li84Sn21 | Li12CoSn3 | metallic           |     84 |          3 |               3 |           5   |     1.0695 |         0.076  |       0.39  |   0.312 |     5.1e-05  |     0.000189 |     0.0003   |     0.000463 |
| Cu-Li    | Cu20Li80    | Li4Cu     | metallic           |     80 |          3 |               2 |           5   |     1.0145 |         0.0775 |       0.79  |   0.245 |     0.00014  |     0.000298 |     0.000445 |     0.000888 |
| Br-Li    | Br25Li75    | Li3Br     | halide             |     75 |          3 |               2 |           7.5 |     1.072  |         0.067  |       0.39  |   0.235 |     0.000171 |     0.000366 |     0.00063  |     0.000892 |
| Li-Si-Y  | Li80Si10Y10 | Li8YSi    | network(B/C/Si/Ge) |     80 |          4 |               2 |           7.5 |     1.0475 |         0.055  |       0.795 |   0.254 |     0.000104 |     0.000241 |     0.000385 |     0.000664 |
| Cd-Li-Si | Cd5Li75Si20 | Li15CdSi4 | network(B/C/Si/Ge) |     75 |          4 |               2 |           5   |     1.0045 |         0.069  |       0.635 |   0.341 |     4.7e-05  |     0.000174 |     0.000368 |     0.000482 |
| Ag-Li    | Ag25Li75    | Li3Ag     | metallic           |     75 |          3 |               2 |           5   |     0.9585 |         0.096  |       0.425 |   0.232 |     0.000103 |     0.000242 |     0.000362 |     0.000545 |

## 7. Top picks - electrolyte chemistry

| system   | formula      | reduced    | anion_class   |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:---------|:-------------|:-----------|:--------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Br-Li    | Br25Li75     | Li3Br      | halide        |     75 |          3 |               2 |           7.5 |     1.072  |         0.067  |       0.39  |   0.235 |     0.000171 |     0.000366 |     0.00063  |     0.000892 |
| Cl-Li-Si | Cl32Li60Si8  | Li15Si2Cl8 | halide        |     60 |          2 |               3 |          15   |     0.9905 |         0.102  |       2.7   |   0.308 |     6.7e-05  |     0.000197 |     0.000426 |     0.000543 |
| Li-S-Sn  | Li81S18Sn9   | Li9SnS2    | chalcogenide  |     81 |          3 |               3 |          15   |     0.9315 |         0.089  |       1.565 |   0.303 |     5.6e-05  |     0.000193 |     0.000268 |     0.00052  |
| Ga-Li-O  | Ga10Li50O40  | Li5GaO4    | oxide         |     50 |          2 |               3 |          15   |     1.0175 |         0.085  |       1     |   0.526 |     7e-06    |     6.1e-05  |     0.00017  |     0.000278 |
| Li-Rb-Se | Li34Rb34Se34 | RbLiSe     | chalcogenide  |     34 |          2 |               3 |          15   |     0.9325 |         0.1245 |       0.585 |   0.474 |     1.3e-05  |     6.3e-05  |     0.000205 |     0.000338 |
| Li-Se    | Li75Se25     | Li3Se      | chalcogenide  |     75 |          3 |               1 |           5   |     0.991  |         0.0635 |       2.81  |   0.282 |     7.7e-05  |     0.000218 |     0.000386 |     0.000556 |
