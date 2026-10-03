# PDAC Obesity Atlas

Static results browser with on-demand Python figures. Serve this directory over HTTP, for example `python -m http.server 8000 --directory web` from the repository root.

The app reads only the curated `atlas_web_datasets` and `atlas_web_rows` API tables. The client uses a publishable key, never a service-role key. Collections without publication approval are not returned by row-level security.

The figure worker loads pinned Pyodide 314.0.7 and Seaborn 0.13.2. Python initializes only when Draw figure is selected. No fitting runs. `plotting.py` preserves the original study heatmap function and colour palette. SVG and PNG downloads are generated from the filtered aggregate rows. Tables also export CSV.

## Lab deployment

The lab site already copies its `public/` directory into its GitHub Pages build. Copy the six runtime files from this directory (index.html, style.css, app.js, config.js, plot-worker.js, plotting.py) into `public/models-obesity/` in the lab website repository. Commit through an account/integration with write access. The existing lab deployment workflow then serves the atlas at `/models-obesity/`. No deployment-workflow edit is needed.

This repository remains the canonical source; repeat the export when the website changes.

## Interpretation

Filters keep distinct effect components, units and contrasts separate. Numerical failures are excluded by default. Heatmap circles and filled forest-plot dots indicate pointwise HDIs excluding zero; they are not multiplicity-controlled discovery markers. The browser includes a dedicated interpretation section.
