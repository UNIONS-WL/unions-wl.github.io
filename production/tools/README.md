# Refreshing the production page

`production/status.html` reads four kinds of file. Each has one command.

## A cosmo_val run

The page shows four runs: data or GLASS mock, tomographic or not.
Each lives in `production/cosmoval/<run-id>/` (`run.json` and its figures), and `production/cosmoval.json` is the index of all of them.
An importer replaces one run and rewrites the index:

| run | command (from `production/tools`) | source |
|---|---|---|
| `data-tomo` | `uv run --with pillow python import_data_tomo.py <run dir> --catalogue … --nz … [--stale-note …]` | a tomographic run directory with `report/figures/`, `report/summary.json`, `provenance.json` |
| `data-2d` | `uv run --with pillow python import_data_2d.py <dir> --catalogue … --code … [--stale …]` | the schema-4 `cosmoval.json` and `cosmoval/*.png` that sp_validation's dashboard publisher writes |
| `mock-tomo` | `uv run --with pillow --with matplotlib --with numpy --with sacc --with pyyaml python import_mock_tomo.py <placement dir> <ensemble root> --seeds <N>` | the sampled GLASS ensemble and its placement report |
| `mock-2d` | `uv run --with pillow python import_mock_2d.py <figures dir> --code … --date …` | the painted-mock construction figures |

The docstring of each importer holds the exact command for the run on the page now.
Then commit `production/cosmoval*` and push.
A run measured before the current n(z), catalogue or code is imported with a stale note, and the page marks it "older inputs".

Blinding: importers copy figures and PTEs only.
Data figures must come from concealed E-mode vectors; nothing derived from an unblinded vector or a blind record goes on the page.

`cvrun.py` holds the run format (sections, tables, verdict chips, provenance labels).

## The plate

`plate.py` draws `production/plate.svg`, the chart that opens the page, from the lensing footprint and the tomographic footprint (HEALPix masks).
Its docstring holds the command; rerun it when either footprint changes.
The page fills in the four run stars from `cosmoval.json`.

## Fair share

`fairshare-push.sh <cluster> <account>…` samples the LevelFS of each account on the cluster it runs on, appends it to `production/fairshare/<cluster>.json` and pushes.
Run it from a clone used for nothing else:

    fairshare-push.sh fir rrg-kyi def-sfabbro def-waerbeke def-mjhudson
    fairshare-push.sh nibi def-sfabbro def-waerbeke def-mjhudson

The page recommends an account only from samples newer than three hours, so each cluster needs this on a timer.

## ShapePipe v2 run 1

`status.json` and `tile_history.json` are the record of the finished run and no longer change.
