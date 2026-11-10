# /// script
# dependencies = ["numpy", "healpy", "astropy", "contourpy", "scipy", "shapely"]
# ///
"""Draw the plate that opens the production page: production/plate.svg.

A north-polar chart (azimuthal equidistant, sky view: RA increases clockwise) of the
UNIONS lensing footprint as two watercolour washes (the whole calibrated sample, and
the part with z-band photometry that carries tomographic bins), the Milky Way as a
dark wash, a graticule, the bright northern stars, and the analysis triangle drawn
over it as a constellation. status.html inlines the file and fills in the run stars.

    srun -n 1 -c 2 --mem=8G uv run production/tools/plate.py \\
        /project/rrg-kyi/cdaley/unions/data/v2.x/cosmo-val-v2-tomo/footprint-nside4096-ring.fits \\
        /project/rrg-kyi/cdaley/unions/analyses/cosmo-val-v2-tomo-nz468/catalogue/footprint-tomo-nside4096-ring.fits
"""
import sys
from pathlib import Path

import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from contourpy import contour_generator
from scipy.ndimage import gaussian_filter
from shapely.geometry import Polygon

SIZE, C, R = 680, 340.0, 300.0   # viewBox side, centre, radius of the rim
DEC_RIM, RA_DOWN = 24.0, 190.0   # declination at the rim; right ascension pointing down
SS = 2                           # raster samples per SVG unit


def xy(ra, dec):
    r = R * (90 - np.asarray(dec, float)) / (90 - DEC_RIM)
    a = np.radians(np.asarray(ra, float) - RA_DOWN + 90)
    return C + r * np.cos(a), C + r * np.sin(a)


def raster_radec():
    g = (np.arange(SIZE * SS) + 0.5) / SS
    x, y = np.meshgrid(g, g)
    r = np.hypot(x - C, y - C)
    dec = 90 - r / R * (90 - DEC_RIM)
    ra = (np.degrees(np.arctan2(y - C, x - C)) - 90 + RA_DOWN) % 360
    return ra, dec, r <= R


def loops(z, level, tol=0.3, min_area=3.0):
    """Closed contour loops of a raster at a level, as one even-odd SVG path."""
    d = []
    for line in contour_generator(z=z).lines(level):
        if len(line) < 4:
            continue
        poly = Polygon(line / SS).simplify(tol)
        if poly.is_empty or poly.area < min_area or poly.geom_type != "Polygon":
            continue
        pts = np.asarray(poly.exterior.coords)[:-1]
        d.append("M" + "L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + "Z")
    return "".join(d)


def footprint(path, ra, dec, inside):
    m = hp.read_map(path, dtype=np.uint8)
    hit = m[hp.ang2pix(hp.get_nside(m), ra, dec, lonlat=True)].astype(float)
    hit[~inside] = 0
    return gaussian_filter(hit, 1.6 * SS)


# Bright northern stars, J2000 degrees: name, RA, Dec, visual magnitude.
STARS = [
    ("Polaris", 37.95, 89.26, 2.0), ("Kochab", 222.68, 74.16, 2.1), ("Pherkad", 230.18, 71.83, 3.0),
    ("zeta UMi", 236.01, 77.79, 4.3), ("eta UMi", 244.38, 75.76, 5.0), ("eps UMi", 251.49, 82.04, 4.2),
    ("Yildun", 263.05, 86.59, 4.4),
    ("Dubhe", 165.93, 61.75, 1.8), ("Merak", 165.46, 56.38, 2.4), ("Phecda", 178.46, 53.69, 2.4),
    ("Megrez", 183.86, 57.03, 3.3), ("Alioth", 193.51, 55.96, 1.8), ("Mizar", 200.98, 54.93, 2.2),
    ("Alkaid", 206.89, 49.31, 1.9),
    ("Caph", 2.29, 59.15, 2.3), ("Schedar", 10.13, 56.54, 2.2), ("gamma Cas", 14.18, 60.72, 2.2),
    ("Ruchbah", 21.45, 60.24, 2.7), ("Segin", 28.60, 63.67, 3.4),
    ("Vega", 279.23, 38.78, 0.0), ("Deneb", 310.36, 45.28, 1.3), ("Capella", 79.17, 46.00, 0.1),
    ("Castor", 113.65, 31.89, 1.6), ("Pollux", 116.33, 28.03, 1.1), ("Alpheratz", 2.10, 29.09, 2.1),
    ("Mirfak", 51.08, 49.86, 1.8), ("Algol", 47.04, 40.96, 2.1), ("Alphecca", 233.67, 26.71, 2.2),
    ("Eltanin", 269.15, 51.49, 2.2), ("Rastaban", 262.61, 52.30, 2.8), ("Alderamin", 319.64, 62.59, 2.5),
    ("Menkalinan", 89.88, 44.95, 1.9), ("Almach", 30.97, 42.33, 2.1), ("Mirach", 17.43, 35.62, 2.1),
    ("Cor Caroli", 194.01, 38.32, 2.9), ("Thuban", 211.10, 64.38, 3.7), ("Sadr", 305.56, 40.26, 2.2),
    ("Albireo", 292.68, 27.96, 3.1), ("Elnath", 81.57, 28.61, 1.7),
]
NAMED = {"Vega": (8, 4), "Deneb": (8, 4), "Capella": (8, 4), "Castor": (8, 4), "Alkaid": (-8, 14)}
ASTERISMS = [
    ["Dubhe", "Merak", "Phecda", "Megrez", "Dubhe"], ["Megrez", "Alioth", "Mizar", "Alkaid"],
    ["Caph", "Schedar", "gamma Cas", "Ruchbah", "Segin"],
    ["Polaris", "Yildun", "eps UMi", "zeta UMi", "Kochab", "Pherkad", "eta UMi", "zeta UMi"],
]
GALAXIES = [("M31", 10.68, 41.27), ("M33", 23.46, 30.66), ("M81", 148.89, 69.07), ("Coma", 194.95, 27.98)]


def main():
    full_path, tomo_path = sys.argv[1:3]
    ra, dec, inside = raster_radec()
    full = loops(footprint(full_path, ra, dec, inside), 0.5)
    tomo = loops(footprint(tomo_path, ra, dec, inside), 0.5)

    gal = SkyCoord(ra * u.deg, dec * u.deg).galactic
    l, b = gal.l.wrap_at(180 * u.deg).deg, gal.b.deg
    glow = np.exp(-(b / (7 + 4 * np.cos(np.radians(l - 70)))) ** 2) * inside
    way = [loops(gaussian_filter(glow, 2 * SS), lv, tol=0.6, min_area=30) for lv in (0.25, 0.6)]

    star = {n: xy(a, d) for n, a, d, _ in STARS}
    out = []
    w = out.append
    w(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE} {SIZE}" role="img" '
      'aria-label="North-polar chart of the UNIONS lensing footprint, with the analysis drawn over it as a '
      'constellation: pipeline at the pole, mocks and data on the survey, and the four validation runs as stars.">')
    w('''<defs>
<filter id="pl-wash" x="-5%" y="-5%" width="110%" height="110%">
  <feTurbulence type="fractalNoise" baseFrequency="0.018" numOctaves="3" seed="7" result="warp"/>
  <feDisplacementMap in="SourceGraphic" in2="warp" scale="11" xChannelSelector="R" yChannelSelector="G" result="edge"/>
  <feTurbulence type="fractalNoise" baseFrequency="0.035 0.05" numOctaves="4" seed="3" result="grain"/>
  <feColorMatrix in="grain" type="matrix" values="0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  1.5 0 0 0 0.12" result="mottle"/>
  <feComposite in="edge" in2="mottle" operator="in" result="pigment"/>
  <feMorphology in="edge" operator="erode" radius="1.6" result="core"/>
  <feComposite in="edge" in2="core" operator="out" result="rim"/>
  <feGaussianBlur in="rim" stdDeviation="0.9" result="pool"/>
  <feMerge><feMergeNode in="pigment"/><feMergeNode in="pool"/></feMerge>
</filter>
<filter id="pl-wash2" x="-5%" y="-5%" width="110%" height="110%">
  <feTurbulence type="fractalNoise" baseFrequency="0.022" numOctaves="3" seed="21" result="warp"/>
  <feDisplacementMap in="SourceGraphic" in2="warp" scale="12" xChannelSelector="R" yChannelSelector="G" result="edge"/>
  <feTurbulence type="fractalNoise" baseFrequency="0.016 0.022" numOctaves="4" seed="13" result="grain"/>
  <feColorMatrix in="grain" type="matrix" values="0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  2.4 0 0 0 -0.32" result="mottle"/>
  <feComposite in="edge" in2="mottle" operator="in" result="pigment"/>
  <feGaussianBlur in="pigment" stdDeviation="0.7"/>
</filter>
<filter id="pl-bleed" x="-10%" y="-10%" width="120%" height="120%">
  <feGaussianBlur stdDeviation="9"/>
</filter>
<filter id="pl-cloud" x="-5%" y="-5%" width="110%" height="110%">
  <feTurbulence type="fractalNoise" baseFrequency="0.011" numOctaves="4" seed="11" result="warp"/>
  <feDisplacementMap in="SourceGraphic" in2="warp" scale="26" xChannelSelector="R" yChannelSelector="G" result="edge"/>
  <feGaussianBlur in="edge" stdDeviation="5" result="soft"/>
  <feTurbulence type="fractalNoise" baseFrequency="0.02 0.028" numOctaves="5" seed="5" result="grain"/>
  <feColorMatrix in="grain" type="matrix" values="0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  2.2 0 0 0 -0.25" result="mottle"/>
  <feComposite in="soft" in2="mottle" operator="in"/>
</filter>
<filter id="pl-paper" x="0" y="0" width="100%" height="100%">
  <feTurbulence type="fractalNoise" baseFrequency="0.9" numOctaves="2" seed="2" result="fibre"/>
  <feColorMatrix in="fibre" type="matrix" values="0 0 0 0 0.35  0 0 0 0 0.27  0 0 0 0 0.16  0 0 0 0.55 -0.2"/>
</filter>
<filter id="pl-ink" x="-2%" y="-2%" width="104%" height="104%">
  <feTurbulence type="fractalNoise" baseFrequency="0.06" numOctaves="2" seed="9" result="warp"/>
  <feDisplacementMap in="SourceGraphic" in2="warp" scale="1.6" xChannelSelector="R" yChannelSelector="G"/>
</filter>
<radialGradient id="pl-vignette" cx="50%" cy="50%" r="50%">
  <stop offset="0.72" stop-color="#8A6A3A" stop-opacity="0"/><stop offset="1" stop-color="#8A6A3A" stop-opacity="0.17"/>
</radialGradient>
<clipPath id="pl-disc"><circle cx="340" cy="340" r="300"/></clipPath>
<path id="pl-star" d="M0,-1 L0.23,-0.23 L1,0 L0.23,0.23 L0,1 L-0.23,0.23 L-1,0 L-0.23,-0.23Z"/>
</defs>''')
    w(f'<circle cx="{C}" cy="{C}" r="{R}" fill="#FBF6E9"/>')
    w('<g clip-path="url(#pl-disc)">')
    w(f'<path d="{way[0]}" fill="#4B4660" fill-opacity="0.34" fill-rule="evenodd" filter="url(#pl-cloud)"/>')
    w(f'<path d="{way[1]}" fill="#2E2A45" fill-opacity="0.38" fill-rule="evenodd" filter="url(#pl-cloud)"/>')
    w(f'<path d="{full}" fill="#57A898" fill-opacity="0.22" fill-rule="evenodd" filter="url(#pl-bleed)"/>')
    w(f'<path d="{full}" fill="#4FA391" fill-opacity="0.7" fill-rule="evenodd" filter="url(#pl-wash)"/>')
    w(f'<path d="{tomo}" fill="#3A6CB4" fill-opacity="0.36" fill-rule="evenodd" filter="url(#pl-wash2)"/>')
    w(f'<rect width="{SIZE}" height="{SIZE}" filter="url(#pl-paper)"/>')
    w(f'<circle cx="{C}" cy="{C}" r="{R}" fill="url(#pl-vignette)"/>')
    w('</g>')

    # Graticule and rim, in pen.
    w('<g class="pl-grat" filter="url(#pl-ink)">')
    for d in (30, 45, 60, 75):
        w(f'<circle cx="{C}" cy="{C}" r="{R * (90 - d) / (90 - DEC_RIM):.1f}"/>')
    for h in range(0, 24, 2):
        (x0, y0), (x1, y1) = xy(15 * h, 82), xy(15 * h, DEC_RIM)
        w(f'<path d="M{x0:.1f},{y0:.1f}L{x1:.1f},{y1:.1f}"/>')
    w('</g>')
    w(f'<circle class="pl-rim" cx="{C}" cy="{C}" r="{R}"/><circle class="pl-rim2" cx="{C}" cy="{C}" r="{R + 22}"/>')
    for h in range(24):
        for k in range(4):
            a = 15 * (h + k / 4)
            (x0, y0), (x1, y1) = xy(a, DEC_RIM), xy(a, DEC_RIM - (1.6 if k else 3.2) * (90 - DEC_RIM) / R * 2)
            w(f'<path class="pl-tick" d="M{x0:.1f},{y0:.1f}L{x1:.1f},{y1:.1f}"/>')
        x, y = xy(15 * h, DEC_RIM - 13.5 * (90 - DEC_RIM) / R)
        w(f'<text class="pl-hour" x="{x:.1f}" y="{y + 3.5:.1f}" text-anchor="middle">{h}<tspan dy="-3" font-size="6.5">h</tspan></text>')
    for d in (30, 45, 60, 75):
        x, y = xy(RA_DOWN - 90 + 3, d)
        w(f'<text class="pl-dec" x="{x + 3:.1f}" y="{y - 3:.1f}">+{d}°</text>')

    # The fixed stars.
    w('<g class="pl-sky">')
    for chain in ASTERISMS:
        w('<path class="pl-figure" d="M' + "L".join(f"{star[n][0]:.1f},{star[n][1]:.1f}" for n in chain) + '"/>')
    for n, a, d, mag in STARS:
        if n == "Polaris":
            continue
        x, y = star[n]
        s = 5.2 - 1.05 * mag
        if mag < 2.3:
            w(f'<use href="#pl-star" transform="translate({x:.1f},{y:.1f}) scale({max(s, 2.6):.1f})"/>')
        else:
            w(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{max(s / 2.4, 0.8):.1f}"/>')
        if n in NAMED:
            dx, dy = NAMED[n]
            w(f'<text class="pl-name" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{"end" if dx < 0 else "start"}">{n}</text>')
    for n, a, d in GALAXIES:
        x, y = xy(a, d)
        w(f'<ellipse class="pl-galaxy" cx="{x:.1f}" cy="{y:.1f}" rx="4.2" ry="2.2" transform="rotate(-30 {x:.1f} {y:.1f})"/>'
          f'<text class="pl-name" x="{x + 7:.1f}" y="{y + 3:.1f}">{n}</text>')
    for label, (a, d), rot in (("Cassiopeia", (14, 66), 0),):
        x, y = xy(a, d)
        w(f'<text class="pl-const" x="{x:.1f}" y="{y:.1f}" text-anchor="middle">{label}</text>')
    x, y = xy(*((SkyCoord(l=0 * u.deg, b=90 * u.deg, frame="galactic").icrs.ra.deg,
                 SkyCoord(l=0 * u.deg, b=90 * u.deg, frame="galactic").icrs.dec.deg)))
    w(f'<path class="pl-tick" d="M{x - 4:.1f},{y:.1f}h8M{x:.1f},{y - 4:.1f}v8"/>'
      f'<text class="pl-name" x="{x - 7:.1f}" y="{y + 3:.1f}" text-anchor="end">Galactic pole</text>')
    w('</g>')
    for text, (a, d), rot in (("Milky Way", (322, 47), 68), ("Milky Way", (70, 38), -58)):
        x, y = xy(a, d)
        w(f'<text class="pl-way" x="{x:.1f}" y="{y:.1f}" text-anchor="middle" transform="rotate({rot} {x:.1f} {y:.1f})">{text}</text>')

    # The analysis, in red ink: pipeline at the pole, mocks and data on the survey.
    P, M, D = star["Polaris"], xy(234, 45), xy(146, 45)
    w('<g class="pl-analysis">')
    w(f'<path class="pl-edge" filter="url(#pl-ink)" d="M{P[0]:.1f},{P[1]:.1f}L{M[0]:.1f},{M[1]:.1f}L{D[0]:.1f},{D[1]:.1f}Z"/>')

    def along(p, q, text, dy, cls="pl-loop"):
        (x0, y0), (x1, y1) = p, q
        if x1 < x0:
            (x0, y0), (x1, y1) = (x1, y1), (x0, y0)
        ang = np.degrees(np.arctan2(y1 - y0, x1 - x0))
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        w(f'<text class="{cls}" x="{mx:.1f}" y="{my:.1f}" dy="{dy}" text-anchor="middle" '
          f'transform="rotate({ang:.1f} {mx:.1f} {my:.1f})">{text}</text>')

    along(P, M, "recovers known inputs", -6)
    along(P, D, "measured blinded", -6)
    along(M, D, "do the data sit inside the mocks?", -7)

    def vertex(p, name, href, anchor, dx, dy, sub):
        x, y = p
        w(f'<a href="{href}"><use class="pl-vertex" href="#pl-star" transform="translate({x:.1f},{y:.1f}) scale(11)"/>'
          f'<circle class="pl-hit" cx="{x:.1f}" cy="{y:.1f}" r="16"/>'
          f'<text class="pl-vname" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}">{name}</text>'
          f'<text class="pl-vsub" x="{x + dx:.1f}" y="{y + dy + 14:.1f}" text-anchor="{anchor}">{sub}</text></a>')

    vertex(P, "PIPELINE", "#compute", "start", 14, -12, '<tspan id="plate-next">compute and fair share</tspan>')
    vertex(M, "MOCKS", "#mock-tomo", "end", -13, -8, "known cosmology")
    vertex(D, "DATA", "#data-tomo", "start", 13, -8, "blinded")

    # The four runs: stars hanging from mocks and data; status.html sets their state and text.
    for run, base, (dx, dy), anchor, label in (
            ("mock-tomo", M, (4, 34), "start", "tomographic"), ("mock-2d", M, (30, 88), "start", "non-tomographic"),
            ("data-tomo", D, (-18, 61), "end", "tomographic"), ("data-2d", D, (-46, 115), "end", "non-tomographic")):
        x, y = base[0] + dx, base[1] + dy
        tx = x + (-12 if anchor == "end" else 12)
        w(f'<a href="#{run}" class="pl-run" data-run="{run}">'
          f'<path class="pl-thread" d="M{base[0]:.1f},{base[1]:.1f}L{x:.1f},{y:.1f}"/>'
          f'<use class="pl-runstar" href="#pl-star" transform="translate({x:.1f},{y:.1f}) scale(8.5)"/>'
          f'<circle class="pl-hit" cx="{x:.1f}" cy="{y:.1f}" r="14"/>'
          f'<text class="pl-rname" x="{tx:.1f}" y="{y - 1:.1f}" text-anchor="{anchor}">{label}</text>'
          f'<text class="pl-rsay" x="{tx:.1f}" y="{y + 13:.1f}" text-anchor="{anchor}"></text></a>')
    w('</g></svg>')

    dest = Path(__file__).resolve().parents[1] / "plate.svg"
    dest.write_text("\n".join(out) + "\n")
    print(dest, f"{dest.stat().st_size / 1e3:.0f} kB")


main()
