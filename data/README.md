# Data provenance

| Folder | Content | Tracked in git |
|---|---|---|
| `hobo_raw/` | HOBO conductivity logger exports, one `.xlsx` per stream-date. Direct input to `vignettes/01_hobo_qaqc.Rmd`. | yes |
| `nutrients/` | Revised NH4-N and SRP grab-sample tables. Input to `vignettes/02_integration.Rmd`. | yes |
| `field/` | Field sheet (`slugs_kauan_clean.csv`), NaCl-SpC calibration, channel widths, upstream geochemistry. | yes |
| `slug_trimmed/` | Per-event slug windows exported by `R/app_slug_interval.R`. Regenerable, tracked so the pipeline runs without re-running the app. | yes |

## Not in this repository

Proprietary HOBOware binaries (`.hobo`, `.hproj`) and the full multi-sensor
deployments (DO, light, metabolism) are archived separately:

> Fonseca, K. (2026). Raw HOBO logger deployments, karst streams dry season 2023.
> Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

The `.xlsx` files in `hobo_raw/` are the HOBOware exports of those binaries and
are sufficient to reproduce every result in this repository.
