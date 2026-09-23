# slug-transient-storage

Breakthrough curves from combined NaCl + nutrient slug additions in
Brazilian karst streams (dry season 2023), processed for transient
storage zone estimation by inverse modelling and for nutrient spiralling
analysis.

## Pipeline

| Vignette | Does | Writes |
|---|---|---|
| `vignettes/01_hobo_qaqc.Rmd` | Instrument QAQC: import HOBO exports, harmonise headers, recalculate SpC, anchor to tracer addition, import Shiny-selected slug windows, flag window quality | `data_derived/btc_conservative.csv`, `data_derived/events_hobo.csv` |
| `vignettes/02_integration.Rmd` | Pair nutrient grab samples to the logger series, reconcile handheld probe vs logger, NaCl from SpC, CaCO₃ co-precipitation correction for SRP | `data_derived/btc_nutrients.csv`, `data_derived/events.csv`, `data_derived/master_tsm.csv` |
| `vignettes/03_btc_review.Rmd` | Visual review only, no processing. Reads the derived tables and renders interactive Plotly panels per event | `docs/03_btc_review.html` |
| `vignettes/04_tsm_uptake.Rmd` | R-native transient-storage model (no OTIS/OTIS-P): two-stage calibration (Runkel 2007) of every event's hydraulics (Stage 1, conservative NaCl BTC) and nutrient uptake (Stage 2), Bonanno et al. (2022)-style identifiability scan, four-simulation mass-balance partition of uptake between main channel and storage zone | `data_derived/tsm_uptake_summary.csv` |

`R/app_slug_interval.R` is a Shiny app used once, between 01's export of
`hobo_with_meta.rds` and its import of the trimmed windows. Its output is
version-controlled, so the vignettes knit without re-running it.

`R/app_tsm_review.R` is a separate, ongoing-use Shiny app for event-by-event
manual curation of the Stage-1/Stage-2 fits: which conservative-tracer
source to trust (logger vs. handheld probe), which grab points to exclude
as field/instrument errors, and eyeballing every fit + identifiability plot
before accepting it. Run with `shiny::runApp("R/app_tsm_review.R")`
(requires `install.packages(c("shiny", "DT"))` once). Every accepted
decision is appended to `data_derived/tsm_manual_overrides.csv` -- nothing
is changed silently. The batch pipeline (`vignettes/04_tsm_uptake.Rmd`,
`R/run_tsm_uptake.R`) and this app currently run independently; the app's
curated decisions are not yet read back into the batch script.

## Data model

Tables joined on `event_id` (`stream_YYYYMMDD_{slug_label|station|single}`):

- **`events.csv`** — one row per event: reach geometry, injected masses,
  addition time, backgrounds, water chemistry, QAQC flags, adopted
  discharge/velocity.
- **`btc_conservative.csv`** — high-resolution logger series (SpC/NaCl over
  time), one row per reading. Dense input to the Stage-1 inverse model for
  logger-sourced events.
- **`btc_nutrients.csv`** — nutrient grab samples in long format, one row
  per sample per solute.
- **`master_tsm.csv`** — the consolidated long table: nutrient
  concentrations (raw, background-corrected, and SRP co-precipitation
  corrected), background levels, and the paired handheld-probe NaCl grabs
  (`nacl_mgL_grab`), one row per event × solute × grab timestamp. Sparse
  input to Stage 1 (probe-sourced events only) and Stage 2 of the inverse
  model.
- **`tsm_uptake_summary.csv`** — one row per event × solute × correction:
  every fitted Stage-1/Stage-2 parameter, the Runkel (2007) uptake
  partition, and the standard Sw/vf/U nutrient-spiralling metrics, from
  `vignettes/04_tsm_uptake.Rmd` / `R/run_tsm_uptake.R`. This is the
  pipeline's current analysis output; see that vignette (and the header
  comments in `R/tsm_model.R`, `tsm_calibrate.R`, `tsm_partition.R`) for
  the full methodology.
- **`tsm_manual_overrides.csv`** — audit log of every decision accepted in
  `R/app_tsm_review.R` (source override, excluded points, resulting
  parameters, RMSE, timestamp). Not yet consumed by the batch pipeline.
- **`coprec_mass_balance.csv`** — the CaCO₃ co-precipitation correction's
  working table (Ca mass balance underlying `p_coprec_ugL` in
  `master_tsm.csv`).
- **`export_excel/tsm_review_workbook_all_campaigns.xlsx`** — a read-only,
  human-friendly companion export, one tab per campaign (`event_id`),
  combining that campaign's hydraulics/metadata (from `events.csv`), its
  conservative-tracer (NaCl) series -- whichever of logger or hand-probe was
  actually used to fit that event's transient-storage hydraulics -- and its
  nutrient grabs, side by side without needing to join the CSVs in R.
  Regenerated from `events.csv`/`btc_conservative.csv`/`master_tsm.csv` by
  `build_review_workbook.py`; edit those source tables, never this file,
  if a correction is needed.

Data dictionaries for each table are in `data_derived/dictionary/`
(`dict_<table_name>.csv`, columns: `column, unit, description`).

## Reproducing

```r
# from the project root, with the .Rproj open
rmarkdown::render("vignettes/01_hobo_qaqc.Rmd", output_dir = "docs")
rmarkdown::render("vignettes/02_integration.Rmd", output_dir = "docs")
rmarkdown::render("vignettes/03_btc_review.Rmd", output_dir = "docs")
rmarkdown::render("vignettes/04_tsm_uptake.Rmd", output_dir = "docs")
```

Paths are resolved with `here::here()` from the project root.

## Conventions

- All timestamps are Brasília time (`America/Sao_Paulo`, UTC−3), the
  logger's native clock. Field notebooks recorded Mato Grosso local time
  (`America/Cuiaba`, UTC−4); conversion happens at import and is
  documented at each occurrence.
- QAQC decisions are recorded as prose notes in the vignettes and as
  explicit flag columns in the derived tables, not as silent filtering.
  The same convention extends to `R/app_tsm_review.R`: every manual
  override is logged to `tsm_manual_overrides.csv`, never applied silently.
- Raw HOBOware binaries are archived separately; see `data/README.md`.

## Citation

See `CITATION.cff`.
