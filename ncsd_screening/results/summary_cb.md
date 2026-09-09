# NCSD amorphous-diffusivity screen  --  convergence & magnitude  (CHARGE-BALANCED COMPOSITIONS ONLY)

60 compositions x 4 temperatures (1000-2500 K), AIMD, 2 fs timestep, 32 ps median production run.

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

Everything below uses the 60 balanced compositions only.

## 1. By temperature

| T (K) | median prod (ps) | median D_Li (cm^2/s) | D_Li 10-90% range | frac diffusive (beta>0.9) | frac Li-converged | frac charge-converged |
|---|---|---|---|---|---|---|
| 1000 | 33 | 9.94e-06 | 4.8e-06 - 4.8e-05 | 0.48 | 0.10 | 0.27 |
| 1500 | 32 | 6.03e-05 | 2.0e-05 - 1.5e-04 | 0.65 | 0.35 | 0.23 |
| 2000 | 31 | 1.45e-04 | 7.9e-05 - 3.0e-04 | 0.77 | 0.28 | 0.27 |
| 2500 | 31 | 2.39e-04 | 1.5e-04 - 4.2e-04 | 0.67 | 0.35 | 0.22 |

## 2. Convergence vs magnitude, by anion chemistry (1000 K)

| anion class | n systems | median D_Li (cm^2/s) | frac Li-converged |
|---|---|---|---|
| chalcogenide | 9 | 3.58e-05 | 0.22 |
| halide | 8 | 9.94e-06 | 0.12 |
| metallic | 1 | 1.96e-05 | 0.00 |
| network(B/C/Si/Ge) | 4 | 5.28e-05 | 0.25 |
| oxide | 33 | 6.92e-06 | 0.03 |
| pnictide | 5 | 1.67e-05 | 0.20 |

## 3. Shortlist A - any chemistry, Li converged at >= 2 T  (20 systems)

Mostly Li metal / Li-alloy melts: large D, textbook convergence. Full table: results/shortlist_all.csv. Top 20 by D_Li(1000 K):

| system    | formula        | reduced   | anion_class        |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:----------|:---------------|:----------|:-------------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Cl-Li     | Cl50Li50       | LiCl      | halide             |     50 |          3 |               0 |             5 |     1.007  |         0.1135 |       1.05  |   0.259 |      8.3e-05 |     0.00026  |     0.000345 |     0.000538 |
| Ge-Li-Sn  | Ge24Li72Sn12   | Li6SnGe2  | network(B/C/Si/Ge) |     72 |          4 |               2 |            15 |     1.066  |         0.095  |       0.33  |   0.295 |      6.8e-05 |     0.000158 |     0.000322 |     0.00056  |
| Li-Te     | Li68Te34       | Li2Te     | chalcogenide       |     68 |          2 |               2 |            15 |     1.045  |         0.119  |       0.83  |   0.266 |      6.3e-05 |     0.000212 |     0.000313 |     0.0004   |
| Bi-Li-S   | Bi15Li45S45    | Li3BiS3   | chalcogenide       |     45 |          3 |               1 |            25 |     0.9575 |         0.1155 |       0.7   |   0.294 |      5e-05   |     0.000158 |     0.000265 |     0.000397 |
| Li-N      | Li75N25        | Li3N      | pnictide           |     75 |          4 |               0 |             5 |     1.0195 |         0.0705 |       1.64  |   0.32  |      4.7e-05 |     0.000155 |     0.000279 |     0.00046  |
| Li-S-Sb   | Li50S40Sb10    | Li5SbS4   | chalcogenide       |     50 |          3 |               2 |            15 |     0.9255 |         0.172  |       0.385 |   0.347 |      4.7e-05 |     0.000145 |     0.000316 |     0.000544 |
| Cu-Li-S   | Cu17Li51S34    | Li3CuS2   | chalcogenide       |     51 |          2 |               1 |            10 |     1.047  |         0.105  |       0.465 |   0.353 |      3.6e-05 |     0.000135 |     0.000274 |     0.00042  |
| Li-S-Ti   | Li36S36Ti9     | Li4TiS4   | chalcogenide       |     36 |          2 |               1 |            20 |     1.0005 |         0.16   |       1.12  |   0.36  |      3.5e-05 |     0.000122 |     0.000304 |     0.00041  |
| Cu-Li-O-P | Cu5Li25O40P10  | Li5CuP2O8 | oxide              |     25 |          2 |               0 |            20 |     1.0405 |         0.1275 |       1.64  |   0.312 |      2.4e-05 |     8.7e-05  |     0.000218 |     0.000168 |
| Cu-Li-Sn  | Cu20Li60Sn20   | Li3CuSn   | metallic           |     60 |          2 |               0 |            10 |     0.941  |         0.1385 |       3.33  |   0.442 |      2e-05   |     0.000132 |     0.000301 |     0.000392 |
| Li-Si-V   | Li51Si34V17    | Li3VSi2   | network(B/C/Si/Ge) |     51 |          2 |               3 |            20 |     0.9795 |         0.116  |       0.51  |   0.418 |      1.9e-05 |     7.4e-05  |     0.000208 |     0.000356 |
| Li-P      | Li75P25        | Li3P      | pnictide           |     75 |          3 |               1 |            25 |     1.0395 |         0.144  |       1.925 |   0.472 |      1.7e-05 |     5.4e-05  |     0.00025  |     0.000425 |
| Li-O-Sn   | Li40O40Sn20    | Li2SnO2   | oxide              |     40 |          2 |               2 |            25 |     0.965  |         0.2115 |       4.16  |   0.384 |      1.5e-05 |     6.5e-05  |     0.000134 |     0.000223 |
| Li-Rb-Se  | Li34Rb34Se34   | RbLiSe    | chalcogenide       |     34 |          2 |               3 |            15 |     0.9325 |         0.1245 |       0.585 |   0.474 |      1.3e-05 |     6.3e-05  |     0.000205 |     0.000338 |
| Al-K-Li-P | Al17K34Li17P34 | K2LiAlP2  | pnictide           |     17 |          2 |               2 |            15 |     0.809  |         0.24   |       1.695 |   0.501 |      9e-06   |     6.7e-05  |     0.000132 |     0.00035  |
| Li-O-Zn   | Li50O45Zn20    | Li10Zn4O9 | oxide              |     50 |          3 |               2 |            15 |     0.971  |         0.1155 |       0.915 |   0.487 |      8e-06   |     6e-05    |     0.000133 |     0.000229 |
| Ga-Li-O   | Ga10Li50O40    | Li5GaO4   | oxide              |     50 |          2 |               3 |            15 |     1.0175 |         0.085  |       1     |   0.526 |      7e-06   |     6.1e-05  |     0.00017  |     0.000278 |
| Li-O-V    | Li20O60V20     | LiVO3     | oxide              |     20 |          2 |               2 |            25 |     0.962  |         0.141  |       1.25  |   0.477 |      7e-06   |     6.6e-05  |     0.000144 |     0.000183 |
| As-Li-O   | As20Li20O60    | LiAsO3    | oxide              |     20 |          2 |               0 |            20 |     0.957  |         0.181  |       1.725 |   0.539 |      6e-06   |     8.3e-05  |     0.000126 |     0.000281 |
| Li-O-Zr   | Li42O49Zr14    | Li6Zr2O7  | oxide              |     42 |          3 |               1 |            20 |     0.9665 |         0.163  |       0.7   |   0.486 |      6e-06   |     3.7e-05  |     9.1e-05  |     0.00017  |

## 4. Shortlist B - insulating-anion (halide/oxide/S/Se/N/P) chemistry, Li converged at >= 2 T  (16 systems)

The relevant set for solid-electrolyte MLIP validation. Full table: results/shortlist_electrolyte.csv.

| system    | formula        | reduced   | anion_class   |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:----------|:---------------|:----------|:--------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Cl-Li     | Cl50Li50       | LiCl      | halide        |     50 |          3 |               0 |             5 |     1.007  |         0.1135 |       1.05  |   0.259 |      8.3e-05 |     0.00026  |     0.000345 |     0.000538 |
| Bi-Li-S   | Bi15Li45S45    | Li3BiS3   | chalcogenide  |     45 |          3 |               1 |            25 |     0.9575 |         0.1155 |       0.7   |   0.294 |      5e-05   |     0.000158 |     0.000265 |     0.000397 |
| Li-N      | Li75N25        | Li3N      | pnictide      |     75 |          4 |               0 |             5 |     1.0195 |         0.0705 |       1.64  |   0.32  |      4.7e-05 |     0.000155 |     0.000279 |     0.00046  |
| Li-S-Sb   | Li50S40Sb10    | Li5SbS4   | chalcogenide  |     50 |          3 |               2 |            15 |     0.9255 |         0.172  |       0.385 |   0.347 |      4.7e-05 |     0.000145 |     0.000316 |     0.000544 |
| Cu-Li-S   | Cu17Li51S34    | Li3CuS2   | chalcogenide  |     51 |          2 |               1 |            10 |     1.047  |         0.105  |       0.465 |   0.353 |      3.6e-05 |     0.000135 |     0.000274 |     0.00042  |
| Li-S-Ti   | Li36S36Ti9     | Li4TiS4   | chalcogenide  |     36 |          2 |               1 |            20 |     1.0005 |         0.16   |       1.12  |   0.36  |      3.5e-05 |     0.000122 |     0.000304 |     0.00041  |
| Cu-Li-O-P | Cu5Li25O40P10  | Li5CuP2O8 | oxide         |     25 |          2 |               0 |            20 |     1.0405 |         0.1275 |       1.64  |   0.312 |      2.4e-05 |     8.7e-05  |     0.000218 |     0.000168 |
| Li-P      | Li75P25        | Li3P      | pnictide      |     75 |          3 |               1 |            25 |     1.0395 |         0.144  |       1.925 |   0.472 |      1.7e-05 |     5.4e-05  |     0.00025  |     0.000425 |
| Li-O-Sn   | Li40O40Sn20    | Li2SnO2   | oxide         |     40 |          2 |               2 |            25 |     0.965  |         0.2115 |       4.16  |   0.384 |      1.5e-05 |     6.5e-05  |     0.000134 |     0.000223 |
| Li-Rb-Se  | Li34Rb34Se34   | RbLiSe    | chalcogenide  |     34 |          2 |               3 |            15 |     0.9325 |         0.1245 |       0.585 |   0.474 |      1.3e-05 |     6.3e-05  |     0.000205 |     0.000338 |
| Al-K-Li-P | Al17K34Li17P34 | K2LiAlP2  | pnictide      |     17 |          2 |               2 |            15 |     0.809  |         0.24   |       1.695 |   0.501 |      9e-06   |     6.7e-05  |     0.000132 |     0.00035  |
| Li-O-Zn   | Li50O45Zn20    | Li10Zn4O9 | oxide         |     50 |          3 |               2 |            15 |     0.971  |         0.1155 |       0.915 |   0.487 |      8e-06   |     6e-05    |     0.000133 |     0.000229 |
| Ga-Li-O   | Ga10Li50O40    | Li5GaO4   | oxide         |     50 |          2 |               3 |            15 |     1.0175 |         0.085  |       1     |   0.526 |      7e-06   |     6.1e-05  |     0.00017  |     0.000278 |
| Li-O-V    | Li20O60V20     | LiVO3     | oxide         |     20 |          2 |               2 |            25 |     0.962  |         0.141  |       1.25  |   0.477 |      7e-06   |     6.6e-05  |     0.000144 |     0.000183 |
| As-Li-O   | As20Li20O60    | LiAsO3    | oxide         |     20 |          2 |               0 |            20 |     0.957  |         0.181  |       1.725 |   0.539 |      6e-06   |     8.3e-05  |     0.000126 |     0.000281 |
| Li-O-Zr   | Li42O49Zr14    | Li6Zr2O7  | oxide         |     42 |          3 |               1 |            20 |     0.9665 |         0.163  |       0.7   |   0.486 |      6e-06   |     3.7e-05  |     9.1e-05  |     0.00017  |

## 5. Cross-check vs MPContribs

- MP `diffusivity` / refit of their own MSD table as slope/6 : median **3.10** (IQR 2.81-3.42, n=803)  -> the published values are 3x the standard 3-D Einstein D.
- MP `diffusivity` / our independent trajectory D_einstein : median **3.21** (n=220)  -> our analysis reproduces their MSD; the 3x offset is the 1D-vs-3D Einstein prefactor (2 D t vs 6 D t), not a physics disagreement.

## 6. Cross-correlations / Onsager terms

The collective (charge) MSD converges far worse than the self MSD: only 22% of systems reach a diffusive collective regime even at 2500 K (vs 35% for self-D). Systems with `charge_conv_T >= 2` in the shortlists are the ones where validating the cross terms against AIMD is also realistic.

## 7. Top picks - any chemistry

| system   | formula      | reduced   | anion_class        |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:---------|:-------------|:----------|:-------------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Ga-Li-O  | Ga10Li50O40  | Li5GaO4   | oxide              |     50 |          2 |               3 |            15 |     1.0175 |         0.085  |       1     |   0.526 |      7e-06   |     6.1e-05  |     0.00017  |     0.000278 |
| Li-Rb-Se | Li34Rb34Se34 | RbLiSe    | chalcogenide       |     34 |          2 |               3 |            15 |     0.9325 |         0.1245 |       0.585 |   0.474 |      1.3e-05 |     6.3e-05  |     0.000205 |     0.000338 |
| Ge-Li-Sn | Ge24Li72Sn12 | Li6SnGe2  | network(B/C/Si/Ge) |     72 |          4 |               2 |            15 |     1.066  |         0.095  |       0.33  |   0.295 |      6.8e-05 |     0.000158 |     0.000322 |     0.00056  |
| Li-Si-V  | Li51Si34V17  | Li3VSi2   | network(B/C/Si/Ge) |     51 |          2 |               3 |            20 |     0.9795 |         0.116  |       0.51  |   0.418 |      1.9e-05 |     7.4e-05  |     0.000208 |     0.000356 |
| Li-Te    | Li68Te34     | Li2Te     | chalcogenide       |     68 |          2 |               2 |            15 |     1.045  |         0.119  |       0.83  |   0.266 |      6.3e-05 |     0.000212 |     0.000313 |     0.0004   |
| Li-S-Sb  | Li50S40Sb10  | Li5SbS4   | chalcogenide       |     50 |          3 |               2 |            15 |     0.9255 |         0.172  |       0.385 |   0.347 |      4.7e-05 |     0.000145 |     0.000316 |     0.000544 |

## 7. Top picks - electrolyte chemistry

| system   | formula      | reduced   | anion_class   |   n_Li |   n_conv_T |   charge_conv_T |   min_trun_ps |   med_beta |   med_block_cv |   med_haven |   Ea_eV |   D_Li_1000K |   D_Li_1500K |   D_Li_2000K |   D_Li_2500K |
|:---------|:-------------|:----------|:--------------|-------:|-----------:|----------------:|--------------:|-----------:|---------------:|------------:|--------:|-------------:|-------------:|-------------:|-------------:|
| Ga-Li-O  | Ga10Li50O40  | Li5GaO4   | oxide         |     50 |          2 |               3 |            15 |     1.0175 |         0.085  |       1     |   0.526 |      7e-06   |     6.1e-05  |     0.00017  |     0.000278 |
| Li-Rb-Se | Li34Rb34Se34 | RbLiSe    | chalcogenide  |     34 |          2 |               3 |            15 |     0.9325 |         0.1245 |       0.585 |   0.474 |      1.3e-05 |     6.3e-05  |     0.000205 |     0.000338 |
| Li-S-Sb  | Li50S40Sb10  | Li5SbS4   | chalcogenide  |     50 |          3 |               2 |            15 |     0.9255 |         0.172  |       0.385 |   0.347 |      4.7e-05 |     0.000145 |     0.000316 |     0.000544 |
| Cu-Li-S  | Cu17Li51S34  | Li3CuS2   | chalcogenide  |     51 |          2 |               1 |            10 |     1.047  |         0.105  |       0.465 |   0.353 |      3.6e-05 |     0.000135 |     0.000274 |     0.00042  |
| Li-O-Zn  | Li50O45Zn20  | Li10Zn4O9 | oxide         |     50 |          3 |               2 |            15 |     0.971  |         0.1155 |       0.915 |   0.487 |      8e-06   |     6e-05    |     0.000133 |     0.000229 |
| Li-N     | Li75N25      | Li3N      | pnictide      |     75 |          4 |               0 |             5 |     1.0195 |         0.0705 |       1.64  |   0.32  |      4.7e-05 |     0.000155 |     0.000279 |     0.00046  |
