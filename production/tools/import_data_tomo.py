"""Import a blinded tomographic cosmo_val data run into the page's data-tomo cell.

Reads <run>/provenance.json, <run>/report/summary.json, <run>/report/figures/*.png
and the fiducial θ cut from <run>/report/figures.py. Copies only the report's
figures and B-mode PTEs: the blind stays in place.

The 31 Oct run (sp_validation run/tomo1030b @ 2e25c76a), from production/tools:

    uv run --with pillow python import_data_tomo.py \\
        /project/rrg-kyi/cdaley/unions/analyses/mock-ensemble-tomo/inputs/data-snapshot-31oct \\
        --catalogue "ShapePipe SP v2.0, upper bound 9 on galaxy size relative to the PSF" \\
        --nz "first working DR6 BPZ+SOM four-bin n(z), validation and null tests only (not the recalibrated n(z) of CosmoStat/sp_validation#468)" \\
        --stale-note "This run predates the recalibrated n(z) bins of CosmoStat/sp_validation#468; the remeasurement on the new bins is running."

Omit --stale-note for a current run.
"""
import argparse
import json
import re
from pathlib import Path

from cvrun import Run

THRESHOLD = 0.05
# key, label, column group, headline label
STATS = [
    ("cosebis", "COSEBIs B", "COSEBIs", "COSEBIs B"),
    ("xip_B", "ξ₊ᴮ", "pure E/B", "ξ₊ᴮ"),
    ("xim_B", "ξ₋ᴮ", "pure E/B", "ξ₋ᴮ"),
    ("combined", "ξ±ᴮ", "pure E/B", "ξ±ᴮ combined"),
    ("cl_bb", "Cℓ BB", "pseudo-Cℓ", "Cℓ BB"),
    ("cl_eb", "Cℓ EB", "pseudo-Cℓ", "Cℓ EB"),
    ("small", "2–20′ ξ±ᴮ", "small scales", "2–20′ ξ±ᴮ"),
]
HEADLINE_ORDER = ["cl_bb", "cl_eb", "cosebis", "xip_B", "xim_B", "combined", "small"]


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
    joint = {k: v[2] for k, v in sm["joint"].items()}
    names = {k: h for k, _, _, h in STATS}
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

    fail = [names[k] for k in HEADLINE_ORDER if joint[k] < THRESHOLD]
    ok = [names[k] for k in HEADLINE_ORDER if joint[k] >= THRESHOLD]
    listing = lambda xs: " and ".join(", ".join(xs).rsplit(", ", 1))
    summary = (f"Jointly over all ten bin pairs, {listing(ok)} pass; {listing(fail)} fail (PTE < {THRESHOLD})." if fail else
               f"Jointly over all ten bin pairs, every B-mode null passes (PTE ≥ {THRESHOLD}).")

    run = Run("data-tomo", title="UNIONS v2, four tomographic bins (blinded data)", summary=summary,
              banner=f"Blinded (blind: {prov['blind_name']}): E-mode data vectors are concealed; do not interpret cosmology",
              provenance=provenance, status="stale" if args.stale_note else "current", stale_note=args.stale_note)

    for k in HEADLINE_ORDER:
        run.headline(names[k], joint[k], THRESHOLD)
    run.doc["headline_note"] = "PTE of each B-mode null test, jointly over the ten bin pairs"

    rows = [{"label": f"bins {p.replace('x', '×')}", "values": {k: v.get(k) for k, *_ in STATS}}
            for p, v in sm["pairs"].items()]
    rows.append({"label": "joint", "note": "all ten pairs, full cross-pair covariance",
                 "values": {k: joint.get(k) for k, *_ in STATS}})
    run.table("bmodes", [{"key": k, "label": lab, "group": g} for k, lab, g, _ in STATS], rows, threshold=THRESHOLD,
              caption=(f"Probability to exceed (PTE) the measured χ² if the B modes are zero, for each bin pair and "
                       f"jointly over all ten pairs. Scale cuts: COSEBIs modes n ≤ {n_modes} and pure-mode ξ±ᴮ at {cut}′; "
                       f"pseudo-Cℓ over all {n_bands} bandpowers; the small-scale test is pure-mode ξ±ᴮ at 2–20′. "
                       f"Real-space covariance: {sm['xi_cov']}; pseudo-Cℓ covariance: Gaussian (NaMaster). "
                       f"Values below {THRESHOLD} are failures."))

    fig = lambda name: report / "figures" / f"{name}.png"
    cc = sm["cov_check"]

    run.figure(fig("pte_map"), "bmodes",
               "PTE of each B-mode null test for every pair of the four tomographic bins. Blue cells are below 0.05; "
               "dark red cells are above 0.95, which would point to an overestimated covariance. Each panel title "
               "gives the joint PTE over all ten pairs.")
    run.figure(fig("xi_B"), "bmodes",
               f"Pure-mode B-mode correlation functions ξ₊ᴮ (circles) and ξ₋ᴮ (squares) for each bin pair, with errors "
               f"from the analytic covariance. Shaded scales lie outside the {cut}′ range of the test; each panel prints "
               f"its PTEs.")
    run.figure(fig("cosebis_B"), "bmodes",
               f"COSEBIs B modes Bₙ (n ≤ {n_modes}, {cut}′) for each bin pair, in units of their own standard deviation, "
               f"with bands at ±1σ and ±2σ. The modes are correlated, so judge by the PTE rather than by single points.")
    run.figure(fig("cl_BB"), "bmodes",
               f"Pseudo-Cℓ BB (circles) and EB (squares) for each bin pair, with a Gaussian NaMaster covariance built on "
               f"each bin's own map and mask. The PTEs use all {n_bands} bandpowers.")
    run.group("galaxy_mask", "Masking nearby galaxies")
    run.figure(fig("messier_xiB"), "bmodes",
               "Whole-sample ξ±ᴮ as measured (points) and after masking six objects that contaminated the PSF star "
               "sample (M33, M101, M13, M92, Draco, Ursa Minor; crosses): no point moves by more than 0.06σ. This test "
               "comes from an earlier non-tomographic run on 138M galaxies, not from this run.", group="galaxy_mask")

    run.figure(fig("xi_E"), "data_vectors",
               "Blinded θξ₊ and θξ₋ for each bin pair against a Planck 2018 prediction without intrinsic alignments. "
               "The band shows how far the blinding can move the data (S₈ ± 0.075, Ωₘ ± 0.1); it can raise ξ₊ above "
               "about 50′ by up to a factor of two in the high-redshift pairs, so such a rise is not by itself a systematic.")
    run.figure(fig("cl_EE"), "data_vectors",
               "Blinded pseudo-Cℓ EE for each bin pair against the fiducial prediction convolved with each band's window.")
    run.figure(fig("amp_consistency"), "data_vectors",
               "Blinded consistency test: one template amplitude per bin pair, fitted jointly to ξ± at 5–100′ and shown "
               "relative to the common amplitude; the χ² tests whether all ten pairs share one amplitude. Grey: template "
               "without intrinsic alignments; blue: with NLA alignments (A_IA = 1); red: NLA after subtracting the "
               "PSF-leakage term from the data. The test was designed after seeing the blinded data.", size="half")

    run.figure(fig("rho"), "psf",
               "ρ statistics of the PSF stars, which measure the PSF model residuals; all four bins share them.",
               size="half")
    run.figure(fig("tau"), "psf",
               "τ₀, τ₂ and τ₅ (correlations between galaxy shapes and PSF quantities) for each tomographic bin, with "
               "each bin's fitted leakage model as lines.")
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
