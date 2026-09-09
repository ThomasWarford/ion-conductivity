"""Fetch MP amorphous_diffusivity reference values for the multi-temperature
systems, by chemical system (the bulk /contributions pagination only returns the
5000 K subset; the `data__chemsys__contains` filter reaches the rest).

For each (system, T <5000, element): the reported `diffusivity` (their pipeline,
= MSD slope/2 of the 3-D MSD) and, from the attached MSD_dt table, a refit as
slope/6 (standard 3-D Einstein) and slope/2 (their convention).

Writes results/mp_reference.csv  (one row per system, temperature, element).
"""
from __future__ import annotations

import glob
import json
import os
import time
import urllib.request

import numpy as np
import pandas as pd

KEY = os.environ.get("MP_API_KEY", "tCczQJW0ghZ2GyCPefsQpYNeTrslordL")
B = "https://contribs-api.materialsproject.org"
OUT = os.path.join(os.path.dirname(__file__), "results", "mp_reference.csv")
DATA = os.path.join(os.path.dirname(__file__), "..", "data_ncsd", "trajectories")


def get(url, tries=6):
    req = urllib.request.Request(url, headers={"X-API-KEY": KEY, "User-Agent": "ncsd"})
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception as e:  # noqa: BLE001
            print("  retry", i, e, flush=True)
            time.sleep(2 + 3 * i)
    raise RuntimeError(url)


def systems():
    s = sorted({os.path.basename(os.path.dirname(p)).split("chemical_system=")[-1]
                for p in glob.glob(os.path.join(DATA, "temperature=1000K",
                                                "chemical_system=*", "*.jsonl.gz"))})
    return s


def fit_table(tid):
    d = get(f"{B}/tables/{tid}/?_fields=data,columns&data_per_page=1000")
    a = np.array(d["data"], dtype=float)
    t = a[:, 0]
    m = (t >= 0.2 * t[-1]) & (t <= 0.8 * t[-1]) & (t > 0)
    out = {"msd_tmax_ps": float(t[-1])}
    if m.sum() >= 5:
        for j, c in enumerate(d["columns"][1:], 1):
            sl = np.polyfit(t[m], a[m, j], 1)[0]
            out[c.replace("_output", "")] = (sl / 6 * 1e-4, sl / 2 * 1e-4)
    return out


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    recs = []
    syss = systems()
    for k, s in enumerate(syss):
        d = get(f"{B}/contributions/?project=amorphous_diffusivity"
                f"&data__chemsys__contains={s}"
                f"&_fields=id,identifier,formula,data,tables&_limit=60")
        for c in d["data"]:
            cd = c["data"]
            if cd.get("chemsys") != s:
                continue
            T = (cd.get("temperature") or {}).get("value")
            if T is None or T >= 5000:
                continue
            props = cd.get("properties") or {}
            reported = {v["element"]: (v.get("diffusivity") or {}).get("value")
                        for v in props.values() if v.get("element")}
            fit = {}
            tb = [x for x in (c.get("tables") or []) if x.get("name") == "MSD_dt"]
            if tb:
                try:
                    fit = fit_table(tb[0]["id"])
                except Exception as e:  # noqa: BLE001
                    print("  table fail", c["identifier"], e, flush=True)
            for el, dv in reported.items():
                e6, e2 = fit.get(el, (np.nan, np.nan))
                recs.append(dict(system=s, temperature=T, element=el,
                                 identifier=c["identifier"], formula=c.get("formula"),
                                 D_reported_cm2s=dv,
                                 D_table_einstein_cm2s=e6, D_table_mpconv_cm2s=e2,
                                 msd_tmax_ps=fit.get("msd_tmax_ps")))
        if (k + 1) % 25 == 0:
            print(f"  {k + 1}/{len(syss)} systems", flush=True)

    df = pd.DataFrame(recs)
    df.to_csv(OUT, index=False)
    li = df[df.element == "Li"]
    r = (li.D_reported_cm2s / li.D_table_einstein_cm2s).replace([np.inf, -np.inf], np.nan).dropna()
    print(f"wrote {OUT}  ({len(df)} rows; {li.temperature.nunique()} temps)")
    print(f"Li  D_reported / (MSD slope/6): median {r.median():.2f}  "
          f"IQR ({r.quantile(.25):.2f}, {r.quantile(.75):.2f})  n={len(r)}")


if __name__ == "__main__":
    main()
