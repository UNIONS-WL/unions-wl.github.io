# /// script
# dependencies = ["numpy", "astropy", "matplotlib"]
# ///
"""An Easter egg on the plate: the PSF stars around M33, production/eggs/m33.png.

    srun -n 1 -c 2 --mem=8G uv run production/tools/egg_m33.py \\
        /project/rrg-kyi/cdaley/unions/data/v2.x/v2.0.0/catalogs/psf_catalog_ngmix.fits production/eggs/m33.png
"""
import sys
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.io import fits
RA0, DEC0, HALF = 23.462, 30.660, 2.0
d = fits.getdata(sys.argv[1], 1)
ra, dec = np.asarray(d["RA"], float), np.asarray(d["Dec"], float)
x = -(ra - RA0) * np.cos(np.radians(DEC0)); y = dec - DEC0     # east to the left
m = (abs(x) < HALF) & (abs(y) < HALF)
x, y = x[m], y[m]
r = np.hypot(x, y)
print("stars in box", m.sum(), "; nearest star to M33 [deg]", r.min())
INK, MUTED, RAG = "#1B1611", "#6E6452", "#FBF6E9"
plt.rcParams.update({"font.family": "serif", "font.size": 11, "text.color": INK, "axes.labelcolor": INK,
                     "axes.edgecolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
cmap = LinearSegmentedColormap.from_list("wash", [RAG, "#8CC2B0", "#3A6CB4", "#2E2A45"])
fig, ax = plt.subplots(figsize=(6.2, 5.4), facecolor=RAG)
n = 48
H, xe, ye = np.histogram2d(x, y, bins=n, range=[[-HALF, HALF]] * 2)
area = (2 * HALF / n * 60) ** 2
im = ax.pcolormesh(xe, ye, (H / area).T, cmap=cmap, rasterized=True)
ax.set_aspect("equal"); ax.set_facecolor(RAG)
ax.set_xlabel("offset east ← → west of M33 [deg]"); ax.set_ylabel("offset north of M33 [deg]")
cb = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.9); cb.set_label("PSF stars per arcmin²"); cb.outline.set_edgecolor(MUTED)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=170)
