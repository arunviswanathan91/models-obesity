# PDAC Obesity Atlas

Static results browser with on-demand Python figures. Serve this directory over HTTP, for example `python -m http.server 8000 --directory web` from the repository root.

The app reads only the curated `atlas_web_datasets` and `atlas_web_rows` API tables. The client uses a publishable key, never a service-role key. Collections without publication approval are not returned by row-level security.

The figure worker loads pinned Pyodide 314.0.7 and Seaborn 0.13.2. Python initializes only when Draw figure is selected. No fitting runs. `plotting.py` preserves the original study heatmap function and colour palette. SVG and PNG downloads are generated from the filtered aggregate rows. Tables also export CSV.

## Lab deployment

The lab deployment workflow checks out a pinned commit of this repository and copies `web/` to `docs/models-obesity/` after the lab site build and checks. The atlas is served at `/models-obesity/`. To deploy a new tested website version, update that pinned commit in the lab workflow using an account or integration with write access.

This repository remains the canonical source. The lab workflow uses a reviewed version rather than automatically taking untested changes.

## Interpretation

Filters keep distinct effect components, units and contrasts separate. Numerical failures are excluded by default. Heatmap circles and filled forest-plot dots indicate pointwise HDIs excluding zero; they are not multiplicity-controlled discovery markers. The browser includes a dedicated interpretation section.
