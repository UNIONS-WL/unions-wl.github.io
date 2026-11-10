"""Import the non-tomographic GLASS mock as the page's mock-2d cell.

No cosmo_val run on a non-tomographic v2 mock has finished, so this cell shows
the construction checks of the seed-1 catalogue-painted mock (sp_validation#459):
the figures and diagnostics.json written by analyses/glass-painted-pr/figures.sh.
When the painted suite on Nibi finishes, replace this importer's source with its
cosmo_val products.

    uv run --with pillow python import_mock_2d.py \\
        /project/rrg-kyi/cdaley/unions/analyses/glass-painted-pr/figures \\
        --code "sp_validation steward/459-rebuilt @ 6c9bac2b" --date "9 Nov 2026"
"""
import argparse
import json
from pathlib import Path

from cvrun import Run

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
parser.add_argument("--code", required=True)
parser.add_argument("--date", required=True)
args = parser.parse_args()
d = json.loads((args.source / "diagnostics.json").read_text())
rows = d["catalogue"]["rows"]

run = Run(
    "mock-2d",
    title="GLASS mock, non-tomographic",
    summary=(f"One painted mock exists and its construction is checked: the sky matches its input spectrum "
             f"(mean ratio {d['spectra']['ee_mean_ratio_ell100_8000']:.3f} over ℓ = 100–8000) and a pure-E sky "
             f"reads as pure E (BB/EE = {d['sign']['correct_bb_over_ee']:.0e}). Its null tests are not measured yet."),
    banner="GLASS mock with a known input cosmology; no data vector is shown",
    provenance={
        "catalogue": f"GLASS seed 1 painted on SP v2.0, size bound 9 ({rows:,} galaxies): real positions, weights and |e|, random orientations",
        "n(z)": "v2 non-tomographic n(z) (data/nz-v2-full-products/nz_v2_nontomo.txt)",
        "code": args.code,
        "date": args.date,
        "blind": "none applied: these figures compare the mock with its own input",
    },
    status="stale",
    stale_note=("No cosmo_val run on a non-tomographic v2 mock has finished. The painted suite on Nibi "
                "(clean, blinded, PSF-leakage and additive-B twins) stopped at 39 of 74 steps on 1 Nov and was "
                "resumed on 10 Nov. Until it lands, this cell shows how the mock is built, not its null tests."),
    lacks="Source clustering correlated with the lensing field; the PSF is the data's own at each galaxy.",
)
run.figure(args.source / "known_answer_spectra.png", "known",
           "The full-sky convergence and E-mode spectra of the seed-1 GLASS sky against the input theory. "
           "The ratio stays at 1 up to ℓ = 8000; a pixel lookup at nside 4096 would follow the dotted curve.")
run.figure(args.source / "catalogue_noise.png", "known",
           "The painted catalogue against its inputs, on one million galaxies: the |e| distribution is unchanged by "
           "the random rotation, the saved shear equals the sky at each position, and the mean painted "
           "ellipticity follows the true shear with slope 1.")
run.figure(args.source / "sign_convention.png", "known",
           "A pure-E sky measured with the pipeline's pseudo-Cℓ estimator. With the catalogue's polarisation "
           "convention it reads as pure E; with the opposite sign it reads as mixed E and B.", size="half")
run.figure(args.source / "reduced_shear_effect.png", "known",
           "The effect of painting reduced shear g = γ/(1−κ) instead of γ on the E-mode spectrum, "
           "and the B-mode power it creates.", size="half")
run.save()
