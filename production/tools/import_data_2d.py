"""Import a non-tomographic blinded data run as the page's data-2d cell.

The source is the manifest sp_validation's dashboard publisher writes: a
directory holding cosmoval.json (schema 4) and cosmoval/*.png.

The run on the page is the bound-9 plain suite of sp_validation run/tomo1101 @ 7081c2d0,
with the run-directory figures redrawn on Fir and the catalogue-level figures from the Nibi render:

    uv run --with pillow python import_data_2d.py \\
        /project/rrg-kyi/cdaley/unions/analyses/production-page/fir-data-2d \\
        --catalogue "ShapePipe SP v2.0, upper bound 9 on galaxy size relative to the PSF; 170.6M galaxies on 20,515 of the 22,901 tiles in the v2.0.0 release, 4,690 deg²" \\
        --nz "for the covariance only: the sum of the four bins of the first working DR6 BPZ+SOM n(z), not the recalibrated n(z) of CosmoStat/sp_validation#468" \\
        --code "sp_validation run/tomo1101 @ 7081c2d0 (measured 1 Nov); figures rendered 10 Nov with feat/page-figures-2d @ 489e3f3c (catalogue-level panels: run/dashboard-tomo1101 @ 54220f69)"

--stale "<what it predates>" marks a run measured on older inputs.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

from cvrun import Run

SECTION = {"data_vectors": "data_vectors", "bmodes": "bmodes", "leakage": "psf",
           "tangential": "tangential", "survey": "survey"}
COLUMN = {"cosebis_n6": "Bₙ, n ≤ 6", "cosebis_n20": "Bₙ, n ≤ 20", "xip_B": "ξ₊ᴮ", "xim_B": "ξ₋ᴮ",
          "combined": "ξ±ᴮ joint", "cl_BB": "Cℓ BB"}

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
parser.add_argument("--catalogue", required=True)
parser.add_argument("--code", required=True)
parser.add_argument("--nz", required=True)
parser.add_argument("--stale", help="what this run predates; omit for a run on current code and inputs")
args = parser.parse_args()

old = json.loads((args.source / "cosmoval.json").read_text())
meta, table = old["metadata"], old["pte_table"]
fiducial = table["rows"][0]
fails = [COLUMN[k] for k, v in fiducial["values"].items() if v < table["threshold"]]
run = Run(
    "data-2d",
    title="Blinded data, non-tomographic",
    summary=(f"At the fiducial cuts ({fiducial['theta']}, {fiducial['ell']}) "
             + (f"{' and '.join(fails)} fail" if fails else "every B-mode null passes")
             + "; the other B-mode nulls pass."
             if fails else "Every B-mode null passes at the fiducial cuts."),
    banner=f"Blinded (blind: {meta['blind_name']}): E-mode data vectors are concealed; do not interpret cosmology",
    provenance={
        "catalogue": args.catalogue,
        "n(z)": args.nz,
        "code": args.code,
        "date": datetime.fromisoformat(meta["time_utc"]).strftime("%-d %b %Y"),
        "blind": meta["blind_name"],
        "n_eff": f"{meta['n_eff']:.2f} arcmin⁻²",
        "σ_e": f"{meta['sigma_e']:.3f} per component",
    },
    status="stale" if args.stale else "current",
    stale_note=args.stale,
)
for key, value in fiducial["values"].items():
    run.headline(COLUMN[key], value, table["threshold"])
run.doc["headline_note"] = f"PTE of each B-mode null test at the fiducial cuts ({fiducial['theta']}, {fiducial['ell']})"
run.table(
    "bmodes",
    [{"key": c["key"], "label": COLUMN[c["key"]], "group": c["group"]} for c in table["columns"]],
    [{"label": r["selection"], "note": f"{r['theta']} · {r['ell']}", "values": r["values"]} for r in table["rows"]],
    caption="Probability to exceed the χ² of each B-mode data vector against zero, with the analytic covariance.",
    threshold=table["threshold"],
)
for group, title in old["groups"].items():
    run.group(group, title)
for f in old["figures"]:
    caption = f["caption"]
    run.figure(args.source / f["url"], SECTION[f["section"]], caption, size=f["size"], group=f["group"],
               alt=caption, fig_id=f["id"],
               narrow=args.source / f["narrow_url"] if f.get("narrow_url") else None)
run.save()
