"""Import the sampled tomographic GLASS mock ensemble as the mock-tomo run.

The ensemble: GLASS mocks sampled on the v2 tomographic footprint with the data's
n(z), per-bin n_eff and sigma_e, measured by the same cosmo_val estimators as the
blinded data, each mock carrying the data run's analytic covariance. The cell asks
whether, on mocks with a known input and no systematics, the B-mode null tests come
back null and the analytic covariance describes the mock scatter.

Inputs are a placement directory written by the recipe's report.py (placement.json
and figures/) and the ensemble's run root (provenance.json, glass_v2_tomo.yaml and
the mock SACCs). Numbers and figures come from the placement directory, so rerunning
on a newer one updates the cell:

    cd production/tools && uv run --with pillow --with matplotlib --with numpy --with sacc \\
        python import_mock_tomo.py \\
        /project/rrg-kyi/cdaley/unions/analyses/mock-ensemble-sampled/placement-n14 \\
        /project/rrg-kyi/cdaley/unions/analyses/mock-ensemble-sampled

The placement report also draws the blinded data's E modes against the mock bands
(band_*_E.png); those are left out. The data_vectors figures here are drawn from the
mock SACCs alone, which carry the throwaway mock blind.
"""
import argparse
import datetime
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sacc
import yaml

from cvrun import Run

MOCK_BLIND = "mock-ensemble-tomo"
PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]
# Joint tests in page order: placement key -> label.
TESTS = {"cl_bb": "Cℓ BB", "cl_eb": "Cℓ EB", "cosebis": "COSEBIs B", "xip_B": "ξ+ B (12–83′)",
         "xim_B": "ξ− B (12–83′)", "combined": "ξ± B (12–83′)", "small": "ξ± B (2–20′)"}
# Placement element keys for the variance map, B from "elements", E from "e_modes".
VARIANCE_ROWS = {"cl_bb": "Cℓ BB", "cl_eb": "Cℓ EB", "cosebis_B": "COSEBIs B", "xip_B": "ξ+ B",
                 "xim_B": "ξ− B", "cl_ee": "Cℓ EE", "cosebis_E": "COSEBIs E", "xip_E": "ξ+ E",
                 "xim_E": "ξ− E"}


def key(p):
    return f"{p[0] + 1}x{p[1] + 1}"


def label(p):
    return f"{p[0] + 1}×{p[1] + 1}"


def git(checkout, *args):
    return subprocess.check_output(["git", "-C", checkout, *args], text=True).strip()


def fig_pte_ecdf(stats, n, path):
    """ECDF of the mocks' analytic PTEs per joint test against the uniform line."""
    fig, axes = plt.subplots(2, 4, figsize=(8, 4.2), sharex=True, sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    for ax, (test, name) in zip(axes.flat, TESTS.items()):
        j = stats[test]["joint"]
        p = np.sort(j["mock_analytic_pte"])
        ax.plot([0, 1], [0, 1], color="0.6", lw=1)
        ax.step(np.concatenate([[0], p, [1]]), np.concatenate([[0], np.arange(1, n + 1) / n, [1]]),
                where="post", color="#2c5d8f", lw=1.5)
        ax.set_title(f"{name}\ndof {j['dof']}, KS p = {j['mock_pte_ks_p']:.2f}", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    axes.flat[-1].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("analytic PTE", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of mocks", fontsize=8)
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def fig_variance(result, path):
    """Mean ensemble/analytic variance per statistic and bin pair."""
    a = np.array([[np.mean(result["elements"][k][key(p)]["diagonal_variance_ratio"])
                   if k in result["elements"] else
                   np.mean(result["e_modes"][k][key(p)]["variance_ratio"])
                   for p in PAIRS] for k in VARIANCE_ROWS])
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    fig.patch.set_facecolor("white")
    im = ax.imshow(a, cmap="PuOr_r", vmin=0.6, vmax=1.4, aspect="auto")
    for i, j in np.ndindex(a.shape):
        ax.text(j, i, f"{a[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(len(PAIRS)), [label(p) for p in PAIRS], fontsize=8)
    ax.set_yticks(range(len(VARIANCE_ROWS)), list(VARIANCE_ROWS.values()), fontsize=8)
    ax.set_xlabel("bin pair", fontsize=8)
    fig.colorbar(im, ax=ax, label="ensemble variance / analytic variance")
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def load_mocks(paths):
    mocks = [sacc.Sacc.load_fits(p) for p in paths]
    assert all(m.metadata.get("blind") == MOCK_BLIND for m in mocks)
    return mocks


def fig_bands(mocks, kinds, tag, scale, xlabel, ylabel, path):
    """Mock-only E-mode bands per bin pair: median and 68% / 95% of the ensemble."""
    fig, axes = plt.subplots(4, 4, figsize=(8, 7.5), layout="constrained")
    fig.patch.set_facecolor("white")
    for i, j in np.ndindex(4, 4):
        if j > i:
            axes[i, j].set_visible(False)
    for p in PAIRS:
        ax = axes[p[1], p[0]]
        for dtype, colour, name in kinds:
            r = mocks[0].indices(dtype, (f"source_{p[0]}", f"source_{p[1]}"))
            x = np.array([mocks[0].data[k].tags[tag] for k in r])
            v = np.array([m.mean[r] for m in mocks]) * scale(x)
            q = np.quantile(v, [.025, .16, .5, .84, .975], axis=0)
            ax.fill_between(x, q[0], q[4], color=colour, alpha=.15, lw=0)
            ax.fill_between(x, q[1], q[3], color=colour, alpha=.35, lw=0)
            ax.plot(x, q[2], color=colour, lw=1, label=name)
        ax.set_xscale("log")
        ax.axhline(0, color="0.6", lw=.5)
        ax.set_title(label(p), fontsize=8)
        ax.tick_params(labelsize=7)
    for i in range(4):
        axes[i, 0].set_ylabel(ylabel, fontsize=8)
        axes[3, i].set_xlabel(xlabel, fontsize=8)
    if len(kinds) > 1:
        axes[0, 0].legend(fontsize=7)
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("placement", type=Path, help="report.py output: placement.json and figures/")
    ap.add_argument("root", type=Path, help="ensemble run root: provenance.json, glass_v2_tomo.yaml")
    args = ap.parse_args()

    result = json.loads((args.placement / "placement.json").read_text())
    prov = json.loads((args.root / "provenance.json").read_text())
    glass = yaml.safe_load((args.root / "glass_v2_tomo.yaml").read_text())
    data_prov = json.loads((Path(prov["data_run"]) / "provenance.json").read_text())
    stats, n = result["statistics"], result["n_mocks"]
    paths = result["inputs"]["mocks"]
    newest = max(Path(p).stat().st_mtime for p in paths)
    date = datetime.datetime.fromtimestamp(newest).strftime("%Y-%m-%d")
    cv = prov["cosmo_val_checkout"]
    branch = git(cv, "rev-parse", "--abbrev-ref", "HEAD")
    s_lo = min(stats[k]["joint"]["s"] for k in TESTS)
    s_hi = max(stats[k]["joint"]["s"] for k in TESTS)
    ks_fail = [TESTS[k] for k in TESTS if stats[k]["joint"]["mock_pte_ks_p"] < 0.05]
    high = [TESTS[k] for k in TESTS if stats[k]["joint"]["s"] - 1 > 2 * stats[k]["joint"]["s_bootstrap_error"]]

    run = Run(
        "mock-tomo",
        title="GLASS mocks, tomographic",
        summary=(f"On {n} sampled GLASS mocks the mean χ²/dof of the seven joint B-mode tests against the "
                 f"analytic covariance lies between {s_lo:.2f} and {s_hi:.2f}"
                 + (f", more than 2σ above 1 for {' and '.join(high)}" if high else "")
                 + "; the mock PTEs are consistent with uniform"
                 + (f" except for {' and '.join(ks_fail)} (KS p < 0.05)." if ks_fail else " for every test.")),
        banner="GLASS mocks under a throwaway blind (mock-ensemble-tomo); the input cosmology is known",
        provenance={
            "catalogue": f"GLASS sampled mock on the v2 tomographic footprint, {n} seeds, nside {glass['nside']}",
            "n(z)": (f"first working DR6 BPZ+SOM n(z) ({Path(glass['nz_path']).name}, {glass['nbins']} bins), "
                     "the one the 31 Oct blinded data run used; not the recalibrated n(z) of CosmoStat/sp_validation#468"),
            "code": (f"generator sp_validation make_glass_sim.py @ {prov['generator_sha'][:8]}; "
                     f"measurement {branch} @ {prov['cosmo_val_sha'][:8]}, with the estimators and analytic covariance of "
                     f"the blinded tomographic data run as rerun at {data_prov['git_sha'][:8]}"),
            "date": date,
            "blind": f"{prov['mock_blind']}, a throwaway mock blind",
        },
        status="stale",
        stale_note=("These mocks were generated with the first working n(z) and its bins, before the recalibration of "
                    "CosmoStat/sp_validation#468. They match the 31 Oct data run, not the remeasurement on the new bins."),
        lacks=("Sampled mocks have unit weights, Gaussian shape noise at the data's per-bin σ_e and n_eff, "
               "no PSF leakage and no additive B mode, so this run has no PSF section."),
    )

    for test, name in TESTS.items():
        run.headline(name, stats[test]["joint"]["mock_pte_ks_p"])
    run.doc["headline_note"] = f"KS p-value that the {n} mock PTEs of each joint B-mode test are uniform"

    rows = []
    for test, name in TESTS.items():
        j = stats[test]["joint"]
        rows.append({"label": name, "values": {
            "median": float(np.median(j["mock_analytic_pte"])), "ks": j["mock_pte_ks_p"]}})
    run.table("bmodes", [{"key": "median", "label": "median mock PTE"},
                         {"key": "ks", "label": "KS p, mock PTEs uniform"}], rows,
              caption=(f"Joint B-mode null tests over all ten bin pairs on {n} mocks, each against the analytic "
                       "covariance. A PTE well below 0.05 means the mocks fail their own null test, or their "
                       "PTEs are not uniform."))
    rows = []
    for test, name in TESTS.items():
        j = stats[test]["joint"]
        rows.append({"label": name, "values": {
            "dof": j["dof"], "s": j["s"], "s_err": j["s_bootstrap_error"],
            "low": int(np.sum(np.asarray(j["mock_analytic_pte"]) < 0.05))}})
    run.table("bmodes", [{"key": "dof", "label": "dof"},
                         {"key": "s", "label": "mean mock χ²/dof"},
                         {"key": "s_err", "label": "± (bootstrap)"},
                         {"key": "low", "label": f"mocks with PTE < 0.05 (of {n})"}], rows,
              caption=("Mean mock χ²/dof is 1 when the analytic covariance describes the scatter of the mocks; "
                       f"about {0.05 * n:.1f} mocks per test should have a PTE below 0.05 by chance."),
              threshold=0)

    out = args.placement / "figures"
    tmp = Path(tempfile.mkdtemp())

    fig_pte_ecdf(stats, n, tmp / "pte_ecdf.png")
    run.figure(tmp / "pte_ecdf.png", "known",
               f"Cumulative distribution of the {n} mocks' PTEs for each joint B-mode test, against the "
               "uniform distribution (grey line) expected when the analytic covariance is right.")
    fig_variance(result, tmp / "variance_map.png")
    lo, hi = result["sampling_variance_bands_gaussian"]["0.68"]
    run.figure(tmp / "variance_map.png", "known",
               f"Variance of the {n} mocks divided by the analytic variance, averaged over each statistic's "
               f"elements, per bin pair. E rows use ξ± at 12–83′ and Cℓ EE at ℓ 300–1600. With {n} Gaussian mocks a "
               f"single element's ratio has a 68% range of {lo:.2f}–{hi:.2f}; these averages scatter less.")
    run.figure(out / "covariance_xi_scale_split.png", "known",
               "Mock variance over analytic variance for each ξ+ B and ξ− B angular bin and bin pair; "
               "the grey band is the 95% scatter expected from a finite number of Gaussian mocks.")
    run.group("covariance", "Covariance diagnostics per test")
    for test, name in TESTS.items():
        run.figure(out / f"covariance_{test}.png", "known",
                   f"{name}: mock variance over analytic variance per element and bin pair (top), and per "
                   "eigenmode of the analytic covariance in ascending eigenvalue (bottom); grey is the 95% "
                   "scatter expected from the finite number of mocks.", group="covariance")

    caption = ("Pointwise 68% and 95% ranges of the {n} mocks for {name} per bin pair, with the blinded "
               "UNIONS v2 data as black points; the blind changes only E modes, so these B-mode points are "
               "the data's own.")
    for stat, name in (("cl_bb", "Cℓ BB"), ("cosebis_B", "COSEBIs B_n"), ("xim_B", "ξ− B")):
        run.figure(out / f"band_{stat}.png", "bmodes", caption.format(n=n, name=name))
    run.group("bands", "More B-mode bands")
    for stat, name in (("xip_B", "ξ+ B"), ("cl_eb", "Cℓ EB")):
        run.figure(out / f"band_{stat}.png", "bmodes", caption.format(n=n, name=name), group="bands")

    mocks = load_mocks(paths)
    fig_bands(mocks, [("galaxy_shear_xi_plus", "#2c5d8f", "ξ+"), ("galaxy_shear_xi_minus", "#c0392b", "ξ−")],
              "theta", lambda x: x * 1e4, "θ [arcmin]", "θ ξ± × 10⁴ [arcmin]", tmp / "mock_xi_E.png")
    run.figure(tmp / "mock_xi_E.png", "data_vectors",
               f"Median and 68% / 95% ranges of θ ξ+ and θ ξ− over the {n} mocks per bin pair, against θ in "
               "arcmin. The mock vectors carry a throwaway blind, so they are shifted from the input cosmology "
               "by an unknown smooth amount.")
    fig_bands(mocks, [("galaxy_shear_cl_ee", "#2c5d8f", "Cℓ EE")], "ell", lambda x: x * 1e7,
              "ℓ", "ℓ Cℓ × 10⁷", tmp / "mock_cl_ee.png")
    run.figure(tmp / "mock_cl_ee.png", "data_vectors",
               f"Median and 68% / 95% ranges of ℓ Cℓ EE over the {n} mocks per bin pair, against ℓ. The mock "
               "vectors carry a throwaway blind, so they are shifted from the input cosmology by an unknown "
               "smooth amount.")

    shutil.rmtree(tmp)
    run.save()


if __name__ == "__main__":
    main()
