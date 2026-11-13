"""Import a non-tomographic blinded data run as the page's data-2d cell.

The source is the manifest sp_validation's dashboard publisher writes: a
directory holding cosmoval.json (schema 4) and cosmoval/*.png.

The run on the page is the whole-sample suite on the all-tiles sample, sp_validation
run/alltiles-2d @ a4257516 (edge 28629e7b plus three fixes), measured and rendered on Fir:

    uv run --with pillow python import_data_2d.py \\
        /project/rrg-kyi/cdaley/unions/analyses/cosmo-val-v2-alltiles-2d/dashboard \\
        --catalogue "ShapePipe v2.0.0 calibrated sample: 191.0M galaxies on all 22,901 run-1 tiles, 5,223 deg²" \\
        --nz "for the covariance only: the sum of the four bins of the recalibrated n(z) of CosmoStat/sp_validation#468" \\
        --code "sp_validation run/alltiles-2d @ a4257516 (edge 28629e7b plus three fixes), measured 13 Nov"

--stale "<what it predates>" marks a run measured on older inputs.
"""
import argparse
import json
from datetime import datetime
from pathlib import Path

from cvrun import Run

SECTION = {"data_vectors": "data_vectors", "bmodes": "bmodes", "leakage": "psf",
           "tangential": "tangential", "survey": "survey"}

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
parser.add_argument("--catalogue", required=True)
parser.add_argument("--code", required=True)
parser.add_argument("--nz", required=True)
parser.add_argument("--stale", help="what this run predates; omit for a run on current code and inputs")
args = parser.parse_args()

old = json.loads((args.source / "cosmoval.json").read_text())
meta = old["metadata"]
run = Run(
    "data-2d",
    title="Blinded data, non-tomographic",
    summary=("The whole calibrated sample as one bin, measured blinded. The B-mode data vectors come first, in three "
             "estimators, each with its errors."),
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
# The PTE-versus-scale-cut maps stay off the page: it shows the B-mode data vectors themselves.
figures = [f for f in old["figures"] if f["group"] != "scale_cuts"]
for group, title in old["groups"].items():
    if any(f["group"] == group for f in figures):
        run.group(group, title)
for f in figures:
    caption = f["caption"]
    run.figure(args.source / f["url"], SECTION[f["section"]], caption, size=f["size"], group=f["group"],
               alt=caption, fig_id=f["id"],
               narrow=args.source / f["narrow_url"] if f.get("narrow_url") else None)
run.save()
