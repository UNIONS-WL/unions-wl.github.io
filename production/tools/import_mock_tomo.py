"""Import a sampled tomographic GLASS mock ensemble as the mock-tomo run.

The ensemble: GLASS mocks sampled on the tomographic footprint with the data's n(z),
per-bin n_eff and sigma_e, measured by the same cosmo_val estimators as the blinded
data, each mock carrying the data run's analytic covariance. The cell asks whether, on
mocks with a known input and no systematics, the B-mode null tests come back null and
the analytic covariance describes the mock scatter, and where the blinded data's B
modes fall inside the ensemble.

Inputs are a placement directory written by the recipe's report.py (placement.json)
and the ensemble's run root (provenance.json, glass_v2_tomo.yaml and the mock SACCs).
Every figure is drawn here from placement.json and the mock SACCs. The run on the page
now is the ensemble matched to the bins and n(z) of CosmoStat/sp_validation#468, with
65 of its 70 seeds:

    cd production/tools && uv run --with pillow --with matplotlib --with numpy --with sacc --with pyyaml \\
        python import_mock_tomo.py \\
        /project/rrg-kyi/cdaley/unions/analyses/mock-ensemble-sampled-nz468/placement-prelim \\
        /project/rrg-kyi/cdaley/unions/analyses/mock-ensemble-sampled-nz468 --seeds 70

The cell is labelled preliminary while the placement holds fewer mocks than --seeds.

Blinding: data enter only through B-mode and EB statistics (placement z-scores, null
PTEs, the Cℓ^BB additive-template amplitude); the blind changes only E modes. The E-mode
figures are drawn from the mock SACCs alone, which carry the throwaway mock blind.
"""
import argparse
import datetime
import json
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
# B-mode element statistics: label, x-axis label, x-axis kind.
ELEMENTS = {"cl_bb": ("Cℓ BB", "multipole ℓ", "ell"), "cosebis_B": ("COSEBIs Bₙ", "mode n", "lin"),
            "xim_B": ("ξ− B", "θ [arcmin]", "log"), "xip_B": ("ξ+ B", "θ [arcmin]", "log"),
            "cl_eb": ("Cℓ EB", "multipole ℓ", "ell")}
ADDITIVE = {"clbb_template_common": "all ten pairs", "clbb_template_autos": "four auto-pairs",
            "clbb_template_crosses": "six cross-pairs"}
OFFSETS = {"mean_offset_cl_bb": "Cℓ BB", "mean_offset_xip_B": "ξ+ B (12–83′)",
           "mean_offset_xim_B": "ξ− B (12–83′)", "mean_offset_combined": "ξ± B (12–83′)",
           "mean_offset_small": "ξ± B (2–20′)"}
BLUE, RED, GREY = "#2c5d8f", "#c0392b", "0.55"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})


def key(p):
    return f"{p[0] + 1}x{p[1] + 1}"


def label(p):
    return f"{p[0] + 1}×{p[1] + 1}"


def git(checkout, *args):
    return subprocess.check_output(["git", "-C", checkout, *args], text=True).strip()


def sqrt_ell(ax):
    """Square-root ℓ axis."""
    ax.set_xscale("function", functions=(lambda x: np.sqrt(np.clip(x, 0, None)), lambda x: x ** 2))
    ax.set_xticks([100, 400, 900, 1600])
    ax.set_xticks([i * 100 for i in range(1, 21)], minor=True)
    ax.set_xlim(0, 2100)


def set_x(ax, kind):
    if kind == "ell":
        sqrt_ell(ax)
    elif kind == "log":
        ax.set_xscale("log")


def triangle(figsize=(8, 7.2)):
    """4×4 lower triangle of bin-pair panels; returns fig and {pair: ax}."""
    fig, axes = plt.subplots(4, 4, figsize=figsize, sharex=True, sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    for i, j in np.ndindex(4, 4):
        if j > i:
            axes[i, j].set_visible(False)
    for i in range(4):
        axes[i, i].tick_params(labelbottom=True)
    return fig, axes, {p: axes[p[1], p[0]] for p in PAIRS}


def pair_tag(ax, p):
    ax.set_title(f"bins {label(p)}", pad=3)


def save(fig, path):
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)


def fig_pte_ecdf(stats, n, path):
    """ECDF of the mocks' analytic PTEs per joint test against the uniform line."""
    fig, axes = plt.subplots(2, 4, figsize=(8, 4.4), sharex=True, sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    for ax, (test, name) in zip(axes.flat, TESTS.items()):
        j = stats[test]["joint"]
        p = np.sort(j["mock_analytic_pte"])
        ax.plot([0, 1], [0, 1], color=GREY, lw=1)
        ax.step(np.concatenate([[0], p, [1]]), np.concatenate([[0], np.arange(1, n + 1) / n, [1]]),
                where="post", color=BLUE, lw=1.5)
        ax.set_title(f"{name}\nKS p = {j['mock_pte_ks_p']:.2g}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
    axes.flat[-1].set_visible(False)
    axes[0, 3].tick_params(labelbottom=True)
    axes[0, 3].set_xlabel("PTE")
    for ax in axes[1, :3]:
        ax.set_xlabel("PTE")
    for ax in axes[:, 0]:
        ax.set_ylabel("fraction of mocks")
    save(fig, path)


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
        ax.text(j, i, f"{a[i, j]:.2f}", ha="center", va="center", fontsize=7.5)
    ax.set_xticks(range(len(PAIRS)), [label(p) for p in PAIRS])
    ax.set_yticks(range(len(VARIANCE_ROWS)), list(VARIANCE_ROWS.values()))
    ax.set_xlabel("bin pair")
    fig.colorbar(im, ax=ax, label="mock variance / analytic variance")
    save(fig, path)


def fig_variance_scale(result, band, path):
    """Mock/analytic variance per element against scale: each bin pair thin, their mean bold."""
    fig, axes = plt.subplots(2, 3, figsize=(8, 5), sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    axes.flat[-1].set_visible(False)
    for ax, (stat, (name, xlabel, kind)) in zip(axes.flat, ELEMENTS.items()):
        el = result["elements"][stat]
        x = np.array(el["1x1"]["coordinate"])
        r = np.array([el[key(p)]["diagonal_variance_ratio"] for p in PAIRS])
        ax.axhspan(*band, color="0.9", lw=0)
        ax.axhline(1, color="k", lw=0.6)
        for row in r:
            ax.plot(x, row, color=GREY, lw=0.5, alpha=0.6)
        ax.plot(x, r.mean(axis=0), color=BLUE, lw=1.8, marker="o", ms=2.5)
        set_x(ax, kind)
        ax.set_title(name)
        ax.set_xlabel(xlabel)
    for ax in axes[:, 0]:
        ax.set_ylabel("mock / analytic variance")
    save(fig, path)


def fig_eigen(stats, band, path):
    """Mock variance along each eigenmode of the analytic covariance, per joint test."""
    fig, axes = plt.subplots(2, 4, figsize=(10, 4.8), sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    for ax, (test, name) in zip(axes.flat, TESTS.items()):
        v = np.array(stats[test]["joint"]["eigenmode_variance"])
        ax.axhspan(*band, color="0.9", lw=0)
        ax.axhline(1, color="k", lw=0.6)
        ax.plot(np.arange(len(v)), v, ".", color=BLUE, ms=2.5)
        ax.set_title(name)
    axes.flat[-1].set_visible(False)
    for ax in list(axes[1, :3]) + [axes[0, 3]]:
        ax.set_xlabel("eigenmode (ascending eigenvalue)")
    for ax in axes[:, 0]:
        ax.set_ylabel("mock / analytic variance")
    save(fig, path)


def fig_z(result, stat, path):
    """Blinded-data B-mode elements in units of the mock scatter, per bin pair."""
    name, xlabel, kind = ELEMENTS[stat]
    fig, axes, ax_of = triangle()
    el = result["elements"][stat]
    for p, ax in ax_of.items():
        x = np.array(el[key(p)]["coordinate"])
        z = np.array(el[key(p)]["data_z"])
        ax.axhspan(-2, 2, color="#c9d9ea", lw=0)
        ax.axhspan(-1, 1, color="#94b4d4", lw=0)
        ax.axhline(0, color="k", lw=0.5)
        out = np.abs(z) > 2
        ax.plot(x[~out], z[~out], "o", color="k", ms=2.5)
        ax.plot(x[out], z[out], "o", color=RED, ms=3)
        set_x(ax, kind)
        pair_tag(ax, p)
    zmax = max(np.abs(el[key(p)]["data_z"]).max() for p in PAIRS)
    axes[0, 0].set_ylim(-max(4.5, zmax + 0.6), max(4.5, zmax + 0.6))
    for i in range(4):
        axes[i, 0].set_ylabel("data offset [mock σ]")
        axes[3, i].set_xlabel(xlabel)
    for i in range(3):
        axes[i, i].set_xlabel(xlabel)
    save(fig, path)


def fig_additive(add, n, path):
    """Cℓ^BB additive-template amplitude: mock histogram and data, for all, auto and cross pairs."""
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), sharey=True, layout="constrained")
    fig.patch.set_facecolor("white")
    for ax, (k, name) in zip(axes, ADDITIVE.items()):
        a = add[k]
        ax.hist(a["mocks"], bins=np.linspace(-1.5, 1.5, 25), color=BLUE, alpha=0.6)
        ax.axvline(a["data"], color=RED, lw=2)
        ax.set_title(f"{name}\ndata {a['data']:.2f} ± {a['sigma']:.2f}, {a['rank_above']} of {n} mocks above")
        ax.set_xlabel("template amplitude")
    axes[0].set_ylabel("mocks")
    save(fig, path)


def load_mocks(paths):
    mocks = [sacc.Sacc.load_fits(p) for p in paths]
    assert all(m.metadata.get("blind") == MOCK_BLIND for m in mocks)
    return mocks


def fig_bands(mocks, kinds, tag, scale, kind, xlabel, ylabel, path):
    """Mock-only E-mode bands per bin pair: median and 68% / 95% of the ensemble."""
    fig, axes, ax_of = triangle()
    for p, ax in ax_of.items():
        for dtype, colour, name in kinds:
            r = mocks[0].indices(dtype, (f"source_{p[0]}", f"source_{p[1]}"))
            x = np.array([mocks[0].data[k].tags[tag] for k in r])
            v = np.array([m.mean[r] for m in mocks]) * scale(x)
            q = np.quantile(v, [.025, .16, .5, .84, .975], axis=0)
            ax.fill_between(x, q[0], q[4], color=colour, alpha=.15, lw=0)
            ax.fill_between(x, q[1], q[3], color=colour, alpha=.35, lw=0)
            ax.plot(x, q[2], color=colour, lw=1, label=name)
        set_x(ax, kind)
        ax.axhline(0, color=GREY, lw=.5)
        pair_tag(ax, p)
    for i in range(4):
        axes[i, 0].set_ylabel(ylabel)
        axes[3, i].set_xlabel(xlabel)
    for i in range(3):
        axes[i, i].set_xlabel(xlabel)
    if len(kinds) > 1:
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper right", bbox_to_anchor=(0.97, 0.97))
    save(fig, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("placement", type=Path, help="report.py output: placement.json")
    ap.add_argument("root", type=Path, help="ensemble run root: provenance.json, glass_v2_tomo.yaml")
    ap.add_argument("--seeds", type=int, required=True, help="seeds the ensemble is built to")
    args = ap.parse_args()

    result = json.loads((args.placement / "placement.json").read_text())
    prov = json.loads((args.root / "provenance.json").read_text())
    glass = yaml.safe_load((args.root / "glass_v2_tomo.yaml").read_text())
    data_prov = json.loads((Path(prov["data_run"]) / "provenance.json").read_text())
    stats, n, add = result["statistics"], result["n_mocks"], result["additive"]
    paths = result["inputs"]["mocks"]
    newest = max(Path(p).stat().st_mtime for p in paths)
    date = datetime.datetime.fromtimestamp(newest).strftime("%Y-%m-%d")
    cv = prov["cosmo_val_checkout"]
    branch = git(cv, "rev-parse", "--abbrev-ref", "HEAD")
    preliminary = n < args.seeds
    s = {k: stats[k]["joint"]["s"] for k in TESTS}
    ks_fail = [TESTS[k] for k in TESTS if stats[k]["joint"]["mock_pte_ks_p"] < 0.05]
    t = add["clbb_template_common"]
    lo, hi = result["sampling_variance_bands_gaussian"]["0.68"]
    band95 = result["sampling_variance_bands_gaussian"]["0.95"]

    run = Run(
        "mock-tomo",
        title="GLASS mocks, tomographic" + (f" (preliminary, {n} of {args.seeds} seeds)" if preliminary else ""),
        summary=((f"Preliminary: {n} of {args.seeds} seeds; the rest arrive today. " if preliminary else "")
                 + f"On {n} sampled GLASS mocks matched to the data cell's bins and n(z), the analytic covariance "
                 f"describes the mock scatter of Cℓ BB and Cℓ EB (mean χ²/dof {s['cl_bb']:.2f}, {s['cl_eb']:.2f}) "
                 f"but is {100 * (min(s[k] for k in ('cosebis', 'xip_B', 'xim_B', 'combined', 'small')) - 1):.0f}–"
                 f"{100 * (max(s[k] for k in ('cosebis', 'xip_B', 'xim_B', 'combined', 'small')) - 1):.0f}% low for "
                 "the COSEBIs and ξ± B tests"
                 + (f", so their mock PTEs are not uniform ({', '.join(ks_fail)}: KS p < 0.05)" if ks_fail else "")
                 + f". The blinded data's Cℓ BB additive-template amplitude, {t['data']:.2f} ± {t['sigma']:.2f}, "
                 f"is exceeded by {t['rank_above']} of {n} mocks ({t['z']:.1f}σ above the mock mean): the data's "
                 "B-mode excess is additive, not a property of mocks without systematics."),
        banner=("GLASS mocks under a throwaway blind; the input cosmology is known. The data B modes shown "
                "against them are blinded data, which the blind leaves unchanged"),
        provenance={
            "catalogue": (f"GLASS sampled mocks on the v2 tomographic footprint (4,065 deg²), {n} of {args.seeds} "
                          f"seeds, nside {glass['nside']}, each bin's n_eff and σ_e from the data"),
            "n(z)": ("recalibrated DR6 BPZ+SOM four-bin n(z) and bins of CosmoStat/sp_validation#468, "
                     "the same as the data cell"),
            "code": (f"generator sp_validation make_glass_sim.py @ {prov['generator_sha'][:8]}; "
                     f"measurement {branch} @ {prov['cosmo_val_sha'][:8]}, with the estimators and analytic "
                     f"covariance of the blinded #468 data run ({data_prov['git_branch']} @ {data_prov['git_sha'][:8]})"),
            "date": date,
            "blind": f"{prov['mock_blind']}, a throwaway mock blind; data B modes from the run under blind unions",
        },
        lacks=("Sampled mocks have unit weights, Gaussian shape noise at the data's per-bin σ_e and n_eff, "
               "no PSF leakage and no additive B mode, so this run has no PSF section."),
    )

    for test, name in TESTS.items():
        run.headline(name, stats[test]["joint"]["mock_pte_ks_p"])
    run.doc["headline_note"] = (f"KS p-value that the {n} mock PTEs of each joint B-mode test are uniform, "
                                "as they are when the analytic covariance is right")

    # Known answers: mock-only calibration.
    rows = []
    for test, name in TESTS.items():
        j = stats[test]["joint"]
        rows.append({"label": name, "values": {
            "dof": j["dof"], "s": j["s"], "s_err": j["s_bootstrap_error"],
            "low": int(np.sum(np.asarray(j["mock_analytic_pte"]) < 0.05)),
            "median": float(np.median(j["mock_analytic_pte"]))}})
    run.table("known", [{"key": "dof", "label": "elements"},
                        {"key": "s", "label": "mean mock χ²/dof"},
                        {"key": "s_err", "label": "± (bootstrap)"},
                        {"key": "median", "label": "median mock PTE"},
                        {"key": "low", "label": f"mocks with PTE < 0.05 (of {n})"}], rows,
              caption=("Each joint B-mode test runs over all ten bin pairs on every mock, against the analytic "
                       "covariance. A PTE (probability to exceed) is the chance that pure noise gives a χ² at least "
                       "this large. Mean mock χ²/dof is 1 and the median PTE 0.5 when the analytic covariance "
                       f"describes the mock scatter; about {0.05 * n:.1f} mocks per test fall below 0.05 by chance."),
              threshold=0)

    tmp = Path(tempfile.mkdtemp())
    fig_pte_ecdf(stats, n, tmp / "pte_ecdf.png")
    run.figure(tmp / "pte_ecdf.png", "known",
               f"Cumulative distribution of the {n} mocks' PTEs for each joint B-mode test, against the uniform "
               "distribution (grey diagonal) expected when the analytic covariance is right. A curve above the "
               "diagonal means too many low PTEs: the mocks scatter more than the analytic covariance says.")
    fig_variance(result, tmp / "variance_map.png")
    run.figure(tmp / "variance_map.png", "known",
               f"Variance of the {n} mocks divided by the analytic variance, averaged over each statistic's "
               f"elements, per bin pair. E rows use ξ± at 12–83′ and Cℓ EE at ℓ 300–1600. With {n} Gaussian mocks a "
               f"single element's ratio has a 68% range of {lo:.2f}–{hi:.2f}; these averages scatter less.")
    fig_variance_scale(result, band95, tmp / "variance_scale.png")
    run.figure(tmp / "variance_scale.png", "known",
               "Mock variance over analytic variance for each B-mode element against scale: grey lines are the ten "
               "bin pairs, blue their mean. The shaded band is the 95% range expected for one element from a "
               f"finite set of {n} Gaussian mocks.")
    fig_eigen(stats, band95, tmp / "eigenmodes.png")
    run.group("covariance", "Covariance along eigenmodes")
    run.figure(tmp / "eigenmodes.png", "known",
               "Mock variance along each eigenmode of the analytic covariance of each joint test, in ascending "
               "eigenvalue; the shaded band is the 95% range expected from the finite number of mocks.",
               group="covariance")

    # B-mode nulls: the blinded data placed in the ensemble.
    rows = []
    for test, name in TESTS.items():
        j = stats[test]["joint"]
        rows.append({"label": name, "values": {
            "analytic": j["data_analytic_pte"], "recal": j["recalibrated_pte"], "rank": j["empirical_pte"]}})
    run.table("bmodes", [{"key": "analytic", "label": "data PTE, analytic covariance"},
                         {"key": "recal", "label": "data PTE, covariance rescaled to the mocks"},
                         {"key": "rank", "label": f"fraction of {n} mocks with larger χ²"}], rows,
              caption=("Joint B-mode null tests of the blinded data over all ten bin pairs. The rescaled PTE "
                       "multiplies the analytic covariance by the mocks' mean χ²/dof; the last column counts mocks "
                       "directly."))
    rows = [{"label": f"Cℓ BB template amplitude, {name}", "values": {
        "data": add[k]["data"], "z": add[k]["z"], "above": (add[k]["rank_above"] + 1) / (n + 1)}}
        for k, name in ADDITIVE.items()]
    rows += [{"label": f"mean data/σ, {name}", "values": {
        "data": add[k]["data"], "z": add[k]["z"], "above": (add[k]["rank_above"] + 1) / (n + 1)}}
        for k, name in OFFSETS.items()]
    run.table("bmodes", [{"key": "data", "label": "data"},
                         {"key": "z", "label": "σ above mock mean"},
                         {"key": "above", "label": f"(mocks above + 1) / {n + 1}"}], rows,
              caption=("χ² ignores sign; an additive B mode (PSF leakage, a shared shear offset) pushes B modes "
                       "the same way across scales and pairs. The template amplitude A is the least-squares fit "
                       "of C_BB = A × 2.6×10⁻¹¹ exp[−(ℓ/2000)²] to Cℓ BB with the analytic covariance; "
                       "the mean data/σ averages each element over its analytic standard deviation. The mocks have "
                       "no additive B, so they give each number's null distribution."))
    fig_additive(add, n, tmp / "additive.png")
    run.figure(tmp / "additive.png", "bmodes",
               f"Cℓ BB additive-template amplitude on the {n} mocks (histogram) and on the blinded data (red line), "
               "fitted over all pairs, the auto-pairs and the cross-pairs.")

    z_caption = ("Blinded UNIONS data, {name}: each element minus the mock mean, in units of the mock standard "
                 "deviation, per bin pair (i×j is tomographic bin i with bin j). Shaded bands are ±1σ and ±2σ of "
                 "the {n} mocks; points beyond 2σ are red. The blind changes only E modes, so these are the "
                 "data's own B modes.")
    for stat in ("cl_bb", "cosebis_B", "xim_B"):
        fig_z(result, stat, tmp / f"z_{stat}.png")
        run.figure(tmp / f"z_{stat}.png", "bmodes", z_caption.format(name=ELEMENTS[stat][0], n=n))
    run.group("bands", "More B-mode placements")
    for stat in ("xip_B", "cl_eb"):
        fig_z(result, stat, tmp / f"z_{stat}.png")
        run.figure(tmp / f"z_{stat}.png", "bmodes", z_caption.format(name=ELEMENTS[stat][0], n=n), group="bands")

    # Data vectors: mocks only.
    mocks = load_mocks(paths)
    fig_bands(mocks, [("galaxy_shear_xi_plus", BLUE, "ξ+"), ("galaxy_shear_xi_minus", RED, "ξ−")],
              "theta", lambda x: x * 1e4, "log", "θ [arcmin]", "θ ξ± × 10⁴ [arcmin]", tmp / "mock_xi_E.png")
    blind_note = ("The mock vectors carry a throwaway blind, so they are shifted from the input cosmology by an "
                  "unknown smooth amount; no data are shown.")
    run.figure(tmp / "mock_xi_E.png", "data_vectors",
               f"Median and 68% / 95% ranges of θ ξ+ and θ ξ− over the {n} mocks per bin pair. " + blind_note)
    fig_bands(mocks, [("galaxy_shear_cl_ee", BLUE, "Cℓ EE")], "ell", lambda x: x * 1e7, "ell",
              "multipole ℓ", "ℓ Cℓ × 10⁷", tmp / "mock_cl_ee.png")
    run.figure(tmp / "mock_cl_ee.png", "data_vectors",
               f"Median and 68% / 95% ranges of ℓ Cℓ EE over the {n} mocks per bin pair. " + blind_note)

    run.save()


if __name__ == "__main__":
    main()
