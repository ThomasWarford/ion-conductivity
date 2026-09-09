"""Pull the Materials Project 'amorphous_diffusivity' MPContribs table:
per-contribution reported diffusivity for each element, chemsys, temperature.
Writes mp_diffusivity.csv  (one row per contribution == one trajectory file).
"""
import csv, json, time, urllib.request

KEY = "tCczQJW0ghZ2GyCPefsQpYNeTrslordL"
BASE = "https://contribs-api.materialsproject.org/contributions/"

def get(url):
    req = urllib.request.Request(url, headers={"X-API-KEY": KEY})
    for _ in range(5):
        try:
            return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception as e:
            print("retry", e); time.sleep(3)
    raise RuntimeError(url)

rows = []
page = 1
while True:
    url = (f"{BASE}?project=amorphous_diffusivity&_limit=100&_page={page}"
           "&_fields=id,identifier,formula,data")
    d = get(url)
    for c in d["data"]:
        dat = c.get("data", {})
        props = (dat.get("properties") or {})
        elem_d = {}
        for slot, v in props.items():
            el = v.get("element")
            dv = (v.get("diffusivity") or {}).get("value")
            if el is not None:
                elem_d[el] = dv
        rows.append({
            "contribs_id": c["id"],
            "mp_id": c.get("identifier"),
            "formula": c.get("formula"),
            "chemsys": dat.get("chemsys"),
            "temperature": (dat.get("temperature") or {}).get("value"),
            "diffusivity_cm2s": json.dumps(elem_d),
        })
    print(f"page {page}/{d['total_pages']}  rows={len(rows)}")
    if page >= d["total_pages"]:
        break
    page += 1

with open("mp_diffusivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print("wrote mp_diffusivity.csv", len(rows))
