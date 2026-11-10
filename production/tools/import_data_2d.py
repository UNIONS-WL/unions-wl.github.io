"""Import a non-tomographic blinded data run as the page's data-2d cell.

The source is the manifest sp_validation's dashboard publisher writes: a
directory holding cosmoval.json (schema 4) and cosmoval/*.png.

The run on the page now is run-int1028, as published in commit 4b2b060:

    git archive 4b2b060 production/cosmoval.json production/cosmoval | tar -x -C /tmp/int1028
    uv run --with pillow python import_data_2d.py /tmp/int1028/production \\
        --catalogue "ShapePipe SP v2.0, upper bound 3 on galaxy size relative to the PSF, 20,515 of 22,901 tiles, 4,655 deg²" \\
        --nz "the n(z) of the v1.4.6.3 catalogue, used as a stand-in; not calibrated for v2" \\
        --code "sp_validation integration/2026-10-28 @ fbb2b3e8" \\
        --stale "Measured on 28 Oct, before the last 2,386 tiles were published, with the size bound of 3 \\
(now 9, sp_validation#448) and the v1.4.6.3 n(z) as a stand-in for the covariance."
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

# Captions that must say the E modes are blinded.
CAPTION = {"bmode_pure_eb": "Blinded ξ± decomposed into pure E and B modes. Shaded: outside the fiducial 12′–83′; "
                            "full range 1′–250′."}

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
    caption = CAPTION.get(f["id"], f["caption"])
    run.figure(args.source / f["url"], SECTION[f["section"]], caption, size=f["size"], group=f["group"],
               alt=caption, fig_id=f["id"],
               narrow=args.source / f["narrow_url"] if f.get("narrow_url") else None)
run.save()
