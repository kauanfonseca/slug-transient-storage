################################################################################
# run_tsm_uptake.R
#
# Batch driver: for every slug event in the pipeline, calibrates the
# transient-storage model against the conservative (NaCl) breakthrough curve,
# then calibrates first-order uptake against each available nutrient (NH4-N,
# SRP) breakthrough curve, partitions uptake between the main channel and the
# storage zone (Runkel 2007), and computes the standard Sw/vf/U uptake
# metrics -- all in R, without OTIS/OTIS-P.
#
# Inputs (unchanged outputs of the existing pipeline, see README.md /
# handoff.md): data_derived/events.csv, data_derived/btc_conservative.csv,
# data_derived/master_tsm.csv, data/nutrients/nutrient_addition_{nitrogen,
# phosphate}_data.csv.
#
# Output: data_derived/tsm_uptake_summary.csv, one row per
# event_id x solute x correction ("raw" for NH4-N and uncorrected SRP;
# "coprec_corrected" for SRP run on conc_corr_coprec_ugL, i.e. with the
# House 1990 co-precipitation correction added back in, approximating the
# BIOTIC-only assimilation fraction -- see the header of tsm_partition.R and
# the project handoff's circularity caveat, section 4.6/7, before using the
# coprec_corrected numbers to argue anything about the co-precipitation
# mechanism itself).
#
# Usage:
#   source("R/tsm_model.R"); source("R/tsm_calibrate.R"); source("R/tsm_partition.R")
#   source("R/run_tsm_uptake.R")
#   summary_tbl <- run_tsm_uptake_all()
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

M_N <- 14.007   # g/mol, elemental nitrogen
M_P <- 30.974   # g/mol, elemental phosphorus

#' Events whose hydraulics should be borrowed from a different event's
#' conservative-tracer fit, because the event itself has no usable logger
#' series. Per Kauan's own instruction (project handoff, Section 7): for
#' SR_20231011_single (the nutrient slug day, no logger), use the transient-
#' storage parameters fitted on the SR_20231009_downstream high-frequency
#' curve (same 72 m reach), but keep SR_20231011_single's OWN discharge
#' (the nutrient-slug day's value) for everything downstream of the
#' hydraulic fit.
HYDRAULICS_BORROWED_FROM <- c(SR_20231011_single = "SR_20231009_downstream")

#' Fit Stage-1 hydraulics for every event that has a usable conservative
#' tracer logger series (has_logger == TRUE and not flag_discharge_invalid).
fit_all_hydraulics <- function(events, btc_conservative, n_lhs = 250, n_cells = 40, seed = 1) {
  set.seed(seed)
  ids <- events %>%
    filter(has_logger, !isTRUE(flag_discharge_invalid), !is.na(discharge_Ls)) %>%
    pull(event_id)

  results <- map(ids, function(eid) {
    e <- events %>% filter(event_id == eid)
    sub <- btc_conservative %>%
      filter(event_id == eid, time_since_release_s >= 0) %>%
      arrange(time_since_release_s)
    if (nrow(sub) < 10) return(NULL)

    L <- e$reach_length_m
    Q <- e$discharge_Ls / 1000
    A <- Q / e$water_velocity_ms
    mass_g <- e$nacl_mass_g

    fit <- tryCatch(
      fit_hydraulics(sub$time_since_release_s, pmax(sub$nacl_mgL, 0), L, Q, A,
                      mass_g, n_lhs = n_lhs, n_cells = n_cells),
      error = function(err) NULL
    )
    if (is.null(fit)) return(NULL)

    list(event_id = eid, L = L, Q = Q, A = A, v = Q / A,
         width = e$reach_mean_width_m, D = fit$par[["D"]],
         alpha = fit$par[["alpha"]], As = fit$par[["As"]], rmse = fit$rmse,
         lhs = fit$lhs)
  })
  names(results) <- ids
  results[!vapply(results, is.null, logical(1))]
}

#' Elemental nutrient mass injected (mg), joined from the raw addition
#' spreadsheets so the molar-mass conversion always matches the salt that was
#' actually weighed out for that event (K2HPO4, NH4Cl, ...), rather than an
#' assumed constant.
lookup_injected_mass_mg <- function(events_row, nitrogen_raw, phosphate_raw, solute) {
  ymd <- format(as.Date(events_row$date), "%Y%m%d")
  if (solute == "NH4-N") {
    r <- nitrogen_raw %>% filter(stream == events_row$stream, date == ymd) %>% slice(1)
    if (nrow(r) == 0) return(NA_real_)
    r$added_mass_NH4Cl_g * (M_N / r$molar_mass_NH4Cl_nutrient) * 1000
  } else if (solute == "SRP") {
    r <- phosphate_raw %>% filter(stream == events_row$stream, date == ymd) %>% slice(1)
    if (nrow(r) == 0) return(NA_real_)
    r$added_mass_PO4_g * (M_P / r$molar_mass_P_nutrient) * 1000
  } else {
    NA_real_
  }
}

#' Run the full Stage-2 + partition + metrics pipeline for one event x solute
#' x concentration-column combination.
run_one_uptake <- function(event_id, solute, conc_col, hyd, events, master_tsm,
                            nitrogen_raw, phosphate_raw, n_lhs = 250, n_cells = 40) {

  e <- events %>% filter(event_id == !!event_id)
  nut <- master_tsm %>%
    filter(event_id == !!event_id, solute == !!solute, !is.na(.data[[conc_col]])) %>%
    arrange(time_since_release_s)
  if (nrow(nut) < 6) {
    return(tibble(event_id = event_id, solute = solute, conc_col = conc_col,
                   fit_status = "too_few_grabs"))
  }

  mass_mg <- lookup_injected_mass_mg(e, nitrogen_raw, phosphate_raw, solute)
  if (!is.finite(mass_mg)) {
    return(tibble(event_id = event_id, solute = solute, conc_col = conc_col,
                   fit_status = "no_addition_record"))
  }

  L <- hyd$L; Q <- hyd$Q; A <- hyd$A; v <- hyd$v; width <- hyd$width
  hydraulics <- c(D = hyd$D, alpha = hyd$alpha, As = hyd$As)

  obs_conc <- pmax(nut[[conc_col]], 0)
  fitU <- tryCatch(
    fit_uptake(nut$time_since_release_s, obs_conc, L, Q, A, hydraulics, mass_mg,
               n_lhs = n_lhs, n_cells = n_cells),
    error = function(err) NULL
  )
  if (is.null(fitU)) {
    return(tibble(event_id = event_id, solute = solute, conc_col = conc_col,
                   fit_status = "uptake_fit_failed"))
  }

  part <- tryCatch(
    partition_uptake(L, Q, A, hyd$D, hyd$alpha, hyd$As,
                      fitU$par[["lambda"]], fitU$par[["lambda_s"]], mass_mg,
                      n_cells = n_cells),
    error = function(err) NULL
  )

  depth_main <- A / width
  depth_storage <- hyd$As / width
  Camb <- suppressWarnings(mean(nut$background_ugL, na.rm = TRUE))

  met <- uptake_metrics(v, depth_main, depth_storage, fitU$par[["lambda"]],
                         fitU$par[["lambda_s"]], hyd$alpha, A, hyd$As, Camb)

  # --- identifiability / reliability flags -----------------------------
  # A fitted coefficient sitting at (or within 5% in log10 space of) the
  # search bound is a classic sign of a non-identifiable or boundary
  # solution (Bonanno et al. 2022): the optimiser is being pushed as far as
  # it is allowed to go, not converging on an interior optimum. Likewise, a
  # trivially small total uptake means lambda/lambda_s are both ~0 and the
  # main-channel/storage-zone SPLIT of that near-zero removal is meaningless
  # noise, even though the arithmetic still returns a percentage.
  lam_bounds <- c(-7, -1); lams_bounds <- c(-7, -1)
  near_bound <- function(val, bounds, tol = 0.05) {
    lv <- log10(val)
    (lv - bounds[1]) < tol * diff(bounds) || (bounds[2] - lv) < tol * diff(bounds)
  }
  lambda_at_bound   <- near_bound(fitU$par[["lambda"]], lam_bounds)
  lambda_s_at_bound <- near_bound(fitU$par[["lambda_s"]], lams_bounds)
  uptake_significant <- is.finite(part$pct_total_uptake) && part$pct_total_uptake > 2

  tibble(
    event_id = event_id, solute = solute, conc_col = conc_col,
    fit_status = "ok", n_grabs = nrow(nut), mass_injected_mg = mass_mg,
    background_ugL = Camb,
    D_m2s = hyd$D, alpha_1s = hyd$alpha, As_m2 = hyd$As, A_m2 = A,
    As_over_A = hyd$As / A, hydraulic_rmse = hyd$rmse,
    lambda_1s = fitU$par[["lambda"]], lambda_s_1s = fitU$par[["lambda_s"]],
    uptake_rmse = fitU$rmse,
    pct_total_uptake = part$pct_total_uptake,
    pct_mainchannel = if (uptake_significant) part$pct_mainchannel else NA_real_,
    pct_storagezone = if (uptake_significant) part$pct_storagezone else NA_real_,
    Sw_mainchannel_m = if (!lambda_at_bound) met$Sw_mainchannel_m else NA_real_,
    vf_mainchannel_mps = met$vf_mainchannel_mps,
    U_mainchannel_ugm2h = met$U_mainchannel * 3600,
    vf_storagezone_mps = met$vf_storagezone_mps, U_storagezone_ugm2h = met$U_storagezone * 3600,
    Sw_total_m = if (uptake_significant) met$Sw_total_m else NA_real_,
    vf_total_mps = met$vf_total_mps,
    U_total_ugm2h = met$U_total * 3600,
    lambda_at_bound = lambda_at_bound, lambda_s_at_bound = lambda_s_at_bound,
    uptake_significant = uptake_significant
  )
}

#' Full batch pipeline across all events and solutes
run_tsm_uptake_all <- function(data_dir = "data_derived", nutrients_dir = "data/nutrients",
                                n_lhs_hydraulics = 250, n_lhs_uptake = 250, n_cells = 40,
                                seed = 1) {

  events <- read_csv(file.path(data_dir, "events.csv"), show_col_types = FALSE)
  btc_conservative <- read_csv(file.path(data_dir, "btc_conservative.csv"), show_col_types = FALSE)
  master_tsm <- read_csv(file.path(data_dir, "master_tsm.csv"), show_col_types = FALSE)
  nitrogen_raw <- read_csv(file.path(nutrients_dir, "nutrient_addition_nitrogen_data.csv"), show_col_types = FALSE) %>%
    mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_NH4Cl_g, molar_mass_NH4Cl_nutrient)
  phosphate_raw <- read_csv(file.path(nutrients_dir, "nutrient_addition_phosphate_data.csv"), show_col_types = FALSE) %>%
    mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_PO4_g, molar_mass_P_nutrient)

  message("Stage 1 -- fitting transient-storage hydraulics on conservative (NaCl) BTCs...")
  hyd_fits <- fit_all_hydraulics(events, btc_conservative, n_lhs = n_lhs_hydraulics,
                                  n_cells = n_cells, seed = seed)

  # borrow hydraulics for events with no usable logger series (see
  # HYDRAULICS_BORROWED_FROM above), keeping the borrowing event's own Q/A/v
  for (eid in names(HYDRAULICS_BORROWED_FROM)) {
    src <- HYDRAULICS_BORROWED_FROM[[eid]]
    if (!(eid %in% names(hyd_fits)) && src %in% names(hyd_fits) && eid %in% events$event_id) {
      e <- events %>% filter(event_id == eid)
      if (is.na(e$discharge_Ls)) next
      src_fit <- hyd_fits[[src]]
      Q <- e$discharge_Ls / 1000
      A <- Q / e$water_velocity_ms
      hyd_fits[[eid]] <- list(event_id = eid, L = e$reach_length_m, Q = Q, A = A, v = Q / A,
                               width = e$reach_mean_width_m, D = src_fit$D,
                               alpha = src_fit$alpha, As = src_fit$As / src_fit$A * A,
                               rmse = NA_real_, lhs = NULL,
                               hydraulics_borrowed_from = src)
      message(sprintf("  %s: no usable logger; hydraulics borrowed from %s", eid, src))
    }
  }

  message(sprintf("Stage 1 done: %d events with fitted/borrowed hydraulics.", length(hyd_fits)))

  jobs <- master_tsm %>%
    distinct(event_id, solute) %>%
    filter(event_id %in% names(hyd_fits))

  jobs <- jobs %>%
    mutate(conc_col = "conc_corr_ugL") %>%
    bind_rows(
      jobs %>% filter(solute == "SRP") %>% mutate(conc_col = "conc_corr_coprec_ugL")
    )

  message(sprintf("Stage 2 -- fitting uptake for %d event x solute x correction combinations...",
                   nrow(jobs)))

  results <- pmap(jobs, function(event_id, solute, conc_col) {
    message(sprintf("  %s | %s | %s", event_id, solute, conc_col))
    run_one_uptake(event_id, solute, conc_col, hyd_fits[[event_id]], events, master_tsm,
                    nitrogen_raw, phosphate_raw, n_lhs = n_lhs_uptake, n_cells = n_cells)
  })

  out <- bind_rows(results) %>%
    mutate(correction = if_else(conc_col == "conc_corr_coprec_ugL",
                                  "coprec_corrected_biotic_only", "raw")) %>%
    left_join(
      tibble(event_id = names(hyd_fits),
             hydraulics_borrowed_from = map_chr(hyd_fits, function(x)
               if (!is.null(x$hydraulics_borrowed_from)) x$hydraulics_borrowed_from else NA_character_)),
      by = "event_id"
    ) %>%
    left_join(events %>% select(event_id, stream, flag_discharge_invalid, flag_mixing_incomplete,
                                  flag_no_baseline, flag_weak_signal),
               by = "event_id") %>%
    left_join(master_tsm %>% distinct(event_id, flag_broken_pairing), by = "event_id") %>%
    relocate(stream, .after = event_id)

  attr(out, "hydraulic_fits") <- hyd_fits
  out
}
