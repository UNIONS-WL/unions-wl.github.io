"""Import a blinded tomographic cosmo_val data run into the page's data-tomo cell.

Reads <run>/provenance.json, <run>/report/summary.json, <run>/report/figures/*.png
and the fiducial θ cut from <run>/report/figures.py. Copies only the report's
figures: the blind stays in place.

The run on the page (sp_validation run/tomo-nz468 @ 7081c2d0 on the all-tiles sample, the
bins and n(z) of CosmoStat/sp_validation#468), from production/tools:

    uv run --with pillow python import_data_tomo.py \\
        /project/rrg-kyi/cdaley/unions/analyses/cosmo-val-v2-tomo-alltiles \\
        --catalogue "ShapePipe v2.0.0 calibrated sample, on all 22,901 run-1 tiles" \\
        --nz "recalibrated DR6 BPZ+SOM four-bin n(z) and bins of CosmoStat/sp_validation#468, on a dz = 0.05 grid; validation and null tests only"

--stale-note "<what it predates>" marks a run measured on older inputs.
"""
import argparse
import json
import re
from pathlib import Path

from cvrun import Run


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--catalogue", required=True, help="catalogue description; tiles and areas are appended")
    ap.add_argument("--nz", required=True, help="n(z) description; per-bin mean redshifts are appended")
    ap.add_argument("--stale-note", help="marks the run stale, saying what it predates")
    args = ap.parse_args()

    report = args.run_dir / "report"
    prov = json.loads((args.run_dir / "provenance.json").read_text())
    sm = json.loads((report / "summary.json").read_text())
    cut = "–".join(f"{float(x):g}" for x in re.search(r"^CUT = \(([\d.]+), ([\d.]+)\)",
                                                       (report / "figures.py").read_text(), re.M).groups())
    tomo = prov["survey_tomo"]
    n_modes = sm["joint"]["cosebis"][1] // len(sm["pairs"])
    n_bands = sm["joint"]["cl_bb"][1] // len(sm["pairs"])

    zs = ", ".join(f"{z:.3f}" for z in sm["nz_mean"])
    neff = ", ".join(f"{tomo['samples'][f'bin{i}']['n_eff']:.2f}" for i in range(1, 5))
    provenance = {
        "catalogue": (f"{args.catalogue}; {prov['survey']['rows'] / 1e6:.1f}M galaxies, "
                      f"{prov['survey']['area_deg2']:,.0f} deg²; tomographic footprint (tiles with z-band photometry) "
                      f"{tomo['tiles']['tomographic']:,} tiles, {tomo['area_deg2']:,.0f} deg²"),
        "n(z)": f"{args.nz}; ⟨z⟩ = {zs}",
        "code": f"sp_validation {prov['git_branch']} @ {prov['git_sha'][:8]}",
        "date": prov["time_utc"][:10],
        "blind": prov["blind_name"],
        "n_eff per bin": f"{neff} arcmin⁻² on the tomographic footprint",
    }

    summary = ("The four tomographic bins and n(z) of CosmoStat/sp_validation#468, measured blinded. The B-mode data "
               "vectors come first, for all ten bin pairs in three estimators, each with its errors.")

    run = Run("data-tomo", title="UNIONS v2, four tomographic bins (blinded data)", summary=summary,
              banner=f"Blinded (blind: {prov['blind_name']}): E-mode data vectors are concealed; do not interpret cosmology",
              provenance=provenance, status="stale" if args.stale_note else "current", stale_note=args.stale_note)

    fig = lambda name: report / "figures" / f"{name}.png"
    cc = sm["cov_check"]

    run.figure(fig("xi_B"), "bmodes",
               f"Pure-mode B-mode correlation functions θξ₊ᴮ (circles) and θξ₋ᴮ (squares) for each pair of tomographic "
               f"bins, with errors from the analytic covariance. Shaded scales lie outside the {cut}′ range of the test.")
    run.figure(fig("cosebis_B"), "bmodes",
               f"COSEBIs B modes Bₙ (n ≤ {n_modes}, {cut}′) for each pair of tomographic bins, in units of their own "
               f"standard deviation, over bands at ±1σ and ±2σ. "
               f"The modes are correlated, so single points say little.")
    run.figure(fig("cl_BB"), "bmodes",
               f"Pseudo-Cℓ BB (filled circles) and EB (open diamonds) for each pair of tomographic bins, in units of "
               f"their own standard deviation from a Gaussian NaMaster covariance, over bands at ±1σ and ±2σ. The ℓ axis "
               f"is linear in √ℓ, on which the {n_bands} bandpowers are evenly spaced; BB and EB are offset slightly in ℓ "
               f"to keep them apart.")
    run.figure(fig("xi_E"), "data_vectors",
               "Blinded shear correlation functions θξ₊ (filled) and θξ₋ (open) for each pair of tomographic bins, "
               "with errors from the analytic covariance. The blinding shifts these data by an unknown cosmology, so "
               "they are shown alone, without a prediction to compare against.")
    run.figure(fig("cl_EE"), "data_vectors",
               "Blinded pseudo-Cℓ EE power spectra ℓCℓ for each pair of tomographic bins, on one y range and an axis "
               "linear in √ℓ, with errors from the Gaussian NaMaster covariance; shown alone for the same reason.")

    run.figure(fig("rho"), "psf",
               "ρ statistics of the PSF stars, which measure the PSF model residuals; all four bins share them.",
               size="half")
    run.figure(fig("tau"), "psf",
               "Correlations between galaxy shapes and PSF quantities for each tomographic bin (one row per bin): τ₀, "
               "and τ₂ and τ₅ multiplied by θ, with each bin's fitted leakage model as a line.")
    run.figure(fig("leakage"), "psf",
               "PSF-leakage coefficients α, β and η per bin, from a fit of the τ statistics to the ρ templates with ρ held fixed. "
               "The τ estimator biases α low by about 25% on a simulation with known leakage "
               "(CosmoStat/shear_psf_leakage#58), so read these α as lower limits.")

    run.figure(fig("nz"), "survey", f"The four-bin redshift distributions of this run: {args.nz}.", size="half")
    run.group("covariances", "Covariances")
    run.figure(fig("cov_jk"), "survey",
               f"Analytic (OneCovariance) ξ± variance divided by the data's 100-patch jackknife variance, per bin pair "
               f"(auto-correlations coloured, cross-correlations grey); the band is the jackknife's own ±14% scatter. "
               f"The median ratio is {cc['median_small']:.2f} below 10′ (required within 15% of 1) and "
               f"{cc['median_cut']:.2f} at {cut}′.", group="covariances")
    run.figure(fig("cov_xi"), "survey",
               "Correlation matrix of the analytic covariance of the tomographic ξ±, in data-file order.",
               size="half", group="covariances")
    run.figure(fig("cov_cl"), "survey",
               f"Correlation matrix of the Gaussian covariance of pseudo-Cℓ BB: ten bin pairs × {n_bands} bandpowers.",
               size="half", group="covariances")
    run.save()


if __name__ == "__main__":
    main()
