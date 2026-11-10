"""Write one cosmo_val run onto the production page.

A run is a directory production/cosmoval/<run-id>/ holding run.json and its
figures. production/cosmoval.json is the index the page reads: every run.json
in display order. The importers in this directory build a Run and call save().

Run ids are the four cells of the page's switcher:
    data-tomo, data-2d, mock-tomo, mock-2d

Sections, in page order (a run lists only the ones it has figures or tables for):
    known         Known answers (mocks: input recovery, ensemble calibration)
    bmodes        B-mode nulls
    data_vectors  Data vectors
    psf           PSF leakage
    tangential    Tangential shear tests
    survey        Survey, n(z) and covariance
"""
import hashlib
import json
import shutil
from pathlib import Path

from PIL import Image

PRODUCTION = Path(__file__).resolve().parents[1]
SCHEMA = 5
ORDER = ["data-tomo", "data-2d", "mock-tomo", "mock-2d"]
SECTIONS = {
    "known": "Known answers",
    "bmodes": "B-mode nulls",
    "data_vectors": "Data vectors",
    "psf": "PSF leakage",
    "tangential": "Tangential shear tests",
    "survey": "Survey, n(z) and covariance",
}
PROVENANCE = ["catalogue", "n(z)", "code", "date", "blind"]
MAX_WIDTH = 2000


class Run:
    def __init__(self, run_id, *, title, summary, banner, provenance, status="current",
                 stale_note=None, lacks=None):
        """provenance: ordered {label: value}; must carry catalogue, n(z), code, date, blind.
        status: 'current' (latest code and inputs) or 'stale' (stale_note says what it predates).
        lacks: for mocks, one sentence on what the mock does not contain."""
        assert run_id in ORDER, run_id
        missing = [k for k in PROVENANCE if k not in provenance]
        assert not missing, f"provenance lacks {missing}"
        assert status in ("current", "stale") and (status == "current" or stale_note)
        self.dir = PRODUCTION / "cosmoval" / run_id
        if self.dir.exists():
            shutil.rmtree(self.dir)
        self.dir.mkdir(parents=True)
        kind, dim = run_id.split("-")
        self.doc = {
            "id": run_id, "kind": kind, "tomographic": dim == "tomo", "title": title,
            "status": status, "stale_note": stale_note, "summary": summary, "banner": banner,
            "lacks": lacks, "provenance": [[k, str(v)] for k, v in provenance.items()],
            "headline": [], "headline_note": None, "groups": {}, "tables": [], "figures": [],
        }

    def headline(self, label, pte, threshold=0.05):
        """One joint null test for the run's verdict chips; set doc["headline_note"] to say what the number is."""
        self.doc["headline"].append({"label": label, "pte": pte, "fail": pte < threshold})

    def group(self, group_id, title):
        """A drawer of secondary figures inside a section."""
        self.doc["groups"][group_id] = title

    def table(self, section, columns, rows, *, caption, threshold=0.05):
        """columns: [{key, label, group?}]; rows: [{label, note?, values: {key: number|None}}].
        Values below threshold are drawn as failures."""
        assert section in SECTIONS, section
        self.doc["tables"].append({"section": section, "caption": caption, "threshold": threshold,
                                   "columns": columns, "rows": rows})

    def _copy(self, src, stem):
        image = Image.open(src)
        if image.width > MAX_WIDTH:
            image = image.resize((MAX_WIDTH, round(image.height * MAX_WIDTH / image.width)), Image.LANCZOS)
        tmp = self.dir / "tmp.png"
        image.save(tmp, optimize=True)
        name = f"{stem}-{hashlib.sha256(tmp.read_bytes()).hexdigest()[:12]}.png"
        tmp.rename(self.dir / name)
        return f"cosmoval/{self.doc['id']}/{name}", image.width, image.height

    def figure(self, src, section, caption, *, size="full", group=None, alt=None, fig_id=None, narrow=None):
        """Copy a PNG into the run (downscaled to MAX_WIDTH, content-hashed name).
        narrow: an alternative layout of the same figure for phone widths."""
        assert section in SECTIONS, section
        assert size in ("full", "half")
        assert group is None or group in self.doc["groups"], group
        src = Path(src)
        fig_id = fig_id or src.stem
        url, width, height = self._copy(src, fig_id)
        entry = {"id": fig_id, "section": section, "group": group, "size": size, "caption": caption,
                 "alt": alt or caption, "url": url, "width": width, "height": height}
        if narrow:
            entry["narrow_url"], entry["narrow_width"], entry["narrow_height"] = self._copy(narrow, fig_id + "-narrow")
        self.doc["figures"].append(entry)

    def save(self):
        used = {f["section"] for f in self.doc["figures"]} | {t["section"] for t in self.doc["tables"]}
        self.doc["sections"] = [{"id": s, "title": t} for s, t in SECTIONS.items() if s in used]
        (self.dir / "run.json").write_text(json.dumps(self.doc, indent=1, ensure_ascii=False) + "\n")
        build_index()


def build_index():
    """Rewrite production/cosmoval.json from the run directories present."""
    runs = [json.loads(p.read_text()) for r in ORDER
            if (p := PRODUCTION / "cosmoval" / r / "run.json").exists()]
    body = json.dumps(runs, sort_keys=True, ensure_ascii=False)
    index = {"schema_version": SCHEMA, "content_hash": hashlib.sha256(body.encode()).hexdigest(), "runs": runs}
    (PRODUCTION / "cosmoval.json").write_text(json.dumps(index, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    build_index()
