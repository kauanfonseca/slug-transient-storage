# slug-transient-storage

Breakthrough curves from combined NaCl + nutrient slug additions in
Brazilian karst streams (dry season 2023), processed for transient
storage zone estimation by inverse modelling and for nutrient spiralling
analysis.

## Pipeline

| Vignette | Does | Writes |
|---|---|---|
| `vignettes/01_hobo_qaqc.Rmd` | Instrument QAQC: import HOBO exports, harmonise headers, recalculate SpC, anchor to tracer addition, import Shiny-selected slug windows, flag window quality | `data_derived/btc_conservative.csv`, `data_derived/events_hobo.csv` |
| `vignettes/02_integration.Rmd` | Pair nutrient grab samples to the logger series, reconcile handheld probe vs logger, NaCl from SpC, CaCO₃ co-precipitation correction for SRP | `data_derived/btc_nutrients.csv`, `data_derived/events.csv` |
| `vignettes/03_btc_review.Rmd` | Visual review only, no processing. Reads the derived tables and renders interactive Plotly panels per event | `docs/03_btc_review.html` |

`R/app_slug_interval.R` is a Shiny app used once, between 01's export of
`hobo_with_meta.rds` and its import of the trimmed windows. Its output is
version-controlled, so the vignettes knit without re-running it.

## Data model

Three tables joined on `event_id` (`stream_YYYYMMDD_{slug_label|station|single}`):

- **`events.csv`** — one row per event: reach geometry, injected masses,
  addition time, backgrounds, water chemistry, QAQC flags.
- **`btc_conservative.csv`** — high-resolution logger series, one row per
  reading. Input to the inverse model.
- **`btc_nutrients.csv`** — nutrient grab samples in long format, one row
  per sample per solute.

Data dictionaries for each are in `data_derived/dictionary/`.

## Reproducing

```r
# from the project root, with the .Rproj open
rmarkdown::render("vignettes/01_hobo_qaqc.Rmd", output_dir = "docs")
rmarkdown::render("vignettes/02_integration.Rmd", output_dir = "docs")
rmarkdown::render("vignettes/03_btc_review.Rmd", output_dir = "docs")
```

Paths are resolved with `here::here()` from the project root.

## Conventions

- All timestamps are Brasília time (`America/Sao_Paulo`, UTC−3), the
  logger's native clock. Field notebooks recorded Mato Grosso local time
  (`America/Cuiaba`, UTC−4); conversion happens at import and is
  documented at each occurrence.
- QAQC decisions are recorded as prose notes in the vignettes and as
  explicit flag columns in the derived tables, not as silent filtering.
- Raw HOBOware binaries are archived separately; see `data/README.md`.

## Citation

See `CITATION.cff`.
