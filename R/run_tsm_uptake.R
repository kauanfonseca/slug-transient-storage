################################################################################
# run_tsm_uptake.R
#
# Batch driver -- THE METHOD. For every slug event in the pipeline:
#
#   Stage 1 (hydraulics, conservative NaCl curve), applied identically to
#   every curve:
#     a. choose the conservative series: logger, or hand-held-probe grabs
#        where events$discharge_source == "probe", unless hard-coded in
#        CONSERVATIVE_SERIES (e.g. "logger_rescaled": logger shape with the
#        mass scaled to the logger's own recovery);
#     b. probe grabs: fill the unlogged pre-arrival gap with zeros
#        (fill_pre_arrival_gap()); every series: zero the pre-arrival
#        baseline drift (zero_pre_arrival());
#     c. fit TWO candidate models with the exact solver (tsm_model.R),
#        Q fixed from dilution gauging in both:
#          "v_fixed"  : D, alpha, As/A      (v = events$water_velocity_ms)
#          "v_fitted" : D, alpha, As/A, v   (A = Q / v fitted, as OTIS-P)
#     d. select one automatically: "v_fitted" if it lowers AICc by more
#        than AICC_MIN_GAIN AND its v is not at the edge of the search
#        range; otherwise "v_fixed". Everything downstream uses the
#        selected model's D, alpha, As, A and v.
#
#   Why two candidates and a rule, not always v_fitted: events$
#   water_velocity_ms is reach length / time-to-peak, which in a transient-
#   storage model underestimates the main-channel velocity (storage delays
#   the peak). Fitting v fixes that on dense logger curves (CB_20230907:
#   fit error after arrival 0.90 -> 0.63 mg/L, see diag_tsm_shape.R), but
#   on sparse probe-grab curves v and As/A may not be separable. The AICc
#   rule keeps the extra parameter only when the data support it, and the
#   choice is recorded per event, not assumed.
#
#   Stage 2 (uptake): with the selected hydraulics fixed, calibrate first-
#   order uptake against each available nutrient (NH4-N, SRP) breakthrough
#   curve, partition uptake between main channel and storage zone (Runkel
#   2007), and compute the standard Sw/vf/U uptake metrics -- all in R,
#   without OTIS/OTIS-P. Stage 2 and the metrics take A, v from the Stage-1
#   result (hyd$A, hyd$v), never from events$water_velocity_ms.
#
# Inputs (unchanged outputs of the existing pipeline, see README.md /
# handoff.md): data_derived/events.csv, data_derived/btc_conservative.csv,
# data_derived/master_tsm.csv, data/nutrients/nutrient_addition_{nitrogen,
# phosphate}_data.csv.
#
# Outputs (the project's final TSM results):
#   run_tsm_uptake_all() returns the uptake summary (one row per event x
#   solute x correction), with attributes
#     "hydraulic_fits"    list, one element per event (both candidates)
#     "hydraulics_table"  tibble, one row per event: both candidate fits,
#                         AICc, the selected model and QA flags
#   With write_outputs = TRUE it also writes
#     data_derived/tsm_uptake_summary.csv
#     data_derived/tsm_hydraulics_summary.csv
#
# "coprec_corrected" SRP rows use conc_corr_coprec_ugL (House 1990 co-
# precipitation correction added back, approximating BIOTIC-only
# assimilation) -- see the header of tsm_partition.R and the circularity
# caveat (section09 Sec. 9.7) before using those numbers to argue anything
# about the co-precipitation mechanism itself.
#
# Usage:
#   source("R/tsm_model.R"); source("R/tsm_calibrate.R"); source("R/tsm_partition.R")
#   source("R/run_tsm_uptake.R")
#   summary_tbl <- run_tsm_uptake_all(write_outputs = TRUE)
#
# Run time: two Stage-1 fits per event (~20-40 s each with the exact
# solver), so roughly 10-15 min for Stage 1 on the current event set.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(here)
})

M_N <- 14.007   # g/mol, elemental nitrogen
M_P <- 30.974   # g/mol, elemental phosphorus

#' Method settings for Stage 1 (kept in one place so they are reported with
#' the results, not buried in function defaults)
TSM_METHOD <- list(
  model_version    = "exact_v1",   # exact solver + two-candidate selection
  zero_pre_arrival = TRUE,         # zero pre-arrival baseline drift
  onset_frac       = 0.05,         # onset = last point < 5% of peak before it
  v_factor         = 2,            # v search: v_input/2 .. v_input*2
  n_starts         = 3,            # Nelder-Mead starts (each restarted once)
  aicc_min_gain    = 2,            # v_fitted kept only if AICc drops by > 2
  bound_tol        = 0.02          # "at bound" = within 2% of log10 range
)

#' Events whose hydraulics should be borrowed from a different event's
#' conservative-tracer fit, because the event itself has no usable logger
#' series. Per Kauan's own instruction (project handoff, Section 7): for
#' SR_20231011_single (the nutrient slug day, no logger), use the transient-
#' storage parameters fitted on the SR_20231009_downstream high-frequency
#' curve (same 72 m reach), but keep SR_20231011_single's OWN discharge
#' (the nutrient-slug day's value) for everything downstream of the
#' hydraulic fit.
#'
#' RA_20230906_N (2026-09-25, Kauan): reviewed the logger_rescaled fit for
#' both N and P side by side -- P's came out acceptable, N's independent
#' fit did not. Same reach (135 m) and same day, so N now borrows P's D,
#' alpha, As/A and velocity correction instead of being fit on its own; N
#' keeps its own Q/discharge for everything downstream of the hydraulic fit,
#' same mechanism as SR_20231011_single below.
HYDRAULICS_BORROWED_FROM <- c(SR_20231011_single = "SR_20231009_downstream",
                              RA_20230906_N = "RA_20230906_P")

#' Hard-coded conservative-series choice for Stage 1 (overrides the default
#' "logger unless events$discharge_source == 'probe'"). Values:
#'   "logger"          logger curve, Q fixed
#'   "logger_rescaled" logger SHAPE: injected NaCl mass scaled by the
#'                     logger's own recovery, Q * int(C dt) / mass, computed
#'                     on the same zeroed curve the fit sees. Q stays the
#'                     events value (from the probe), v is fitted from the
#'                     logger's timing.
#'   "probe_grab"      hand-held-probe grabs
#'
#' RA_20230906_N / _P (2026-09, diag_ra_logger.R): the logger recovered only
#' 35% (N) and 61% (P) of the salt, and a DIFFERENT fraction for each slug on
#' the same logger and day -- consistent with incomplete lateral mixing at
#' the logger, not a conversion error. With Q fixed the logger fits invent a
#' storage zone 1.6-7x the channel (As/A = 7.2, 1.6). The probe-grab fits
#' disagree with each other (v 0.051 vs 0.038, As/A 0.67 vs 0.37). The
#' rescaled logger fits agree (v 0.041 vs 0.036, As/A 0.31 vs 0.26), so the
#' logger shape is used and its area is not. N is no longer fit
#' independently here -- see HYDRAULICS_BORROWED_FROM above.
CONSERVATIVE_SERIES <- c(RA_20230906_P = "logger_rescaled")

# small helpers ---------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

aicc <- function(rmse, n, k) {
  rss <- n * rmse^2
  n * log(rss / n) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)
}

#' Which fitted parameters sit at (within tol of) a search bound, as a
#' "+"-joined string ("" when none)
params_at_bound <- function(fit, tol = TSM_METHOD$bound_tol) {
  lv <- c(D = log10(fit$par[["D"]]), alpha = log10(fit$par[["alpha"]]),
          As_ratio = log10(fit$par[["As"]] / fit$par[["A"]]),
          v = log10(fit$par[["v"]]))
  hit <- vapply(names(fit$bounds_log10), function(nm) {
    b <- fit$bounds_log10[[nm]]
    (lv[[nm]] - b[1]) < tol * diff(b) || (b[2] - lv[[nm]]) < tol * diff(b)
  }, logical(1))
  paste(names(hit)[hit], collapse = "+")
}

#' Stage 1 for every event with a usable conservative-tracer series:
#' fit both candidate models and select one (see file header).
#'
#' Which series to use is not always "the logger, because has_logger is
#' TRUE": `events$discharge_source` already records, event by event, when
#' the pipeline itself decided the LOGGER's own conductivity curve is not
#' trustworthy (contaminated start anchor, no clean baseline, weak signal --
#' see the project handoff, Section 4.3/4.8) and fell back to the
#' hand-held probe (YSI/Hanna) trapezoidal estimate instead. Whenever
#' discharge_source == "probe", the conservative series used for
#' calibration is the paired hand-held NaCl series (`nacl_mgL_grab` in
#' master_tsm -- the conductivity grabs taken alongside the nutrient
#' samples, already background-corrected and converted via NACL_SLOPE in
#' 02_integration.Rmd). That series is as sparse as the nutrient grabs
#' (~20-30 points), which is exactly where the AICc rule is expected to
#' keep v fixed more often.
fit_all_hydraulics <- function(events, btc_conservative, master_tsm,
                               n_lhs = 250, n_cells = 40, seed = 1,
                               fill_pre_arrival = TRUE, method = TSM_METHOD) {
  # Each event gets its OWN seed, derived from `seed` + its event_id (see
  # job_seed(), tsm_calibrate.R), so a given event's fit does not depend on
  # how many events were fitted before it or in what order.
  ids <- events %>%
    filter(!(flag_discharge_invalid %in% TRUE),   # (isTRUE() on a column
           !is.na(discharge_Ls),                  #  never filtered anything)
           !is.na(water_velocity_ms),
           has_logger | discharge_source == "probe" |
             event_id %in% names(CONSERVATIVE_SERIES),
           !(event_id %in% names(HYDRAULICS_BORROWED_FROM))) %>%
    pull(event_id)
  
  results <- map(ids, function(eid) {
    e <- events %>% filter(event_id == eid)
    series <- if (eid %in% names(CONSERVATIVE_SERIES)) {
      CONSERVATIVE_SERIES[[eid]]
    } else if (isTRUE(e$discharge_source == "probe")) {
      "probe_grab"
    } else {
      "logger"
    }
    stopifnot(series %in% c("logger", "logger_rescaled", "probe_grab"))
    use_probe <- series == "probe_grab"
    
    if (use_probe) {
      sub <- master_tsm %>%
        filter(event_id == eid, !is.na(nacl_mgL_grab), time_since_release_s >= 0) %>%
        distinct(time_since_release_s, nacl_mgL_grab) %>%
        arrange(time_since_release_s) %>%
        rename(nacl_mgL = nacl_mgL_grab)
      source_used <- "probe_grab"
    } else {
      sub <- btc_conservative %>%
        filter(event_id == eid, time_since_release_s >= 0) %>%
        arrange(time_since_release_s) %>%
        distinct(time_since_release_s, .keep_all = TRUE)
      source_used <- series
    }
    if (nrow(sub) < 8) {
      message(sprintf("  %s: fewer than 8 conservative points -- skipped", eid))
      return(NULL)
    }
    
    L <- e$reach_length_m
    Q <- e$discharge_Ls / 1000
    v_input <- e$water_velocity_ms
    A_input <- Q / v_input
    mass_g  <- e$nacl_mass_g
    
    # logger recovery on the zeroed, clamped curve the fit will see
    logger_recovery <- NA_real_
    if (!use_probe) {
      zp <- zero_pre_arrival(sub$time_since_release_s, pmax(sub$nacl_mgL, 0), method$onset_frac)
      m0 <- sum(diff(sub$time_since_release_s) * (head(zp$value, -1) + tail(zp$value, -1)) / 2)
      logger_recovery <- Q * m0 / mass_g
    }
    if (series == "logger_rescaled") {
      if (!is.finite(logger_recovery) || logger_recovery <= 0) {
        message(sprintf("  %s: logger recovery not usable -- skipped", eid))
        return(NULL)
      }
      mass_g <- mass_g * logger_recovery
    }
    
    # Pre-arrival gap fill: only for the hand-held-probe grabs, never the
    # auto-logging conductivity logger (see fill_pre_arrival_gap()).
    obs_conc <- pmax(sub$nacl_mgL, 0)
    gf <- if (use_probe && fill_pre_arrival) {
      fill_pre_arrival_gap(sub$time_since_release_s, obs_conc)
    } else {
      list(time = sub$time_since_release_s, value = obs_conc,
           kind = rep("observed", nrow(sub)), n_added = 0L, dt_used = NA_real_)
    }
    
    fit_one <- function(fit_velocity) {
      tryCatch(
        fit_hydraulics(gf$time, gf$value, L, Q, A_input, mass_g,
                       n_lhs = n_lhs, n_cells = n_cells,
                       seed = job_seed(seed, eid),
                       fit_velocity = fit_velocity, v_factor = method$v_factor,
                       zero_pre_arrival = method$zero_pre_arrival,
                       onset_frac = method$onset_frac,
                       n_starts = method$n_starts),
        error = function(err) {
          message(sprintf("  %s: fit (fit_velocity=%s) failed: %s", eid, fit_velocity,
                          conditionMessage(err)))
          NULL
        }
      )
    }
    f_fix <- fit_one(FALSE)
    f_vf  <- fit_one(TRUE)
    if (is.null(f_fix) && is.null(f_vf)) return(NULL)
    
    # ---- model selection ----------------------------------------------------
    aicc_fix <- if (!is.null(f_fix)) aicc(f_fix$rmse, f_fix$n_fit, 3) else NA_real_
    aicc_vf  <- if (!is.null(f_vf))  aicc(f_vf$rmse,  f_vf$n_fit,  4) else NA_real_
    d_aicc   <- aicc_fix - aicc_vf          # > 0 favours v_fitted
    v_at_bound <- !is.null(f_vf) && grepl("(^|\\+)v($|\\+)", params_at_bound(f_vf))
    
    choose_vf <- !is.null(f_vf) &&
      (is.null(f_fix) || (is.finite(d_aicc) && d_aicc > method$aicc_min_gain && !v_at_bound))
    sel <- if (choose_vf) f_vf else f_fix
    reason <- if (is.null(f_fix)) "v_fixed fit failed"
    else if (is.null(f_vf)) "v_fitted fit failed"
    else if (choose_vf) sprintf("dAICc = %.1f > %g", d_aicc, method$aicc_min_gain)
    else if (v_at_bound) "v_fitted: v at search bound"
    else sprintf("dAICc = %.1f <= %g", d_aicc, method$aicc_min_gain)
    
    message(sprintf("  %-24s %-10s -> %-8s (%s) | v %.4f -> %.4f",
                    eid, source_used, if (choose_vf) "v_fitted" else "v_fixed", reason,
                    v_input, sel$par[["v"]]))
    
    list(event_id = eid, L = L, Q = Q, A = sel$par[["A"]], v = sel$par[["v"]],
         v_input = v_input, width = e$reach_mean_width_m,
         D = sel$par[["D"]], alpha = sel$par[["alpha"]], As = sel$par[["As"]],
         rmse = sel$rmse, lhs = sel$lhs,
         hydraulic_model = if (choose_vf) "v_fitted" else "v_fixed",
         selection_reason = reason,
         aicc_v_fixed = aicc_fix, aicc_v_fitted = aicc_vf, d_aicc = d_aicc,
         params_at_bound = params_at_bound(sel),
         conservative_source = source_used,
         logger_recovery = logger_recovery, nacl_mass_fit_g = mass_g,
         n_gap_filled = gf$n_added, gap_fill_dt_s = gf$dt_used,
         n_pre_zeroed = sel$n_pre_zeroed, t_onset = sel$t_onset, n_fit = sel$n_fit,
         model_version = method$model_version,
         candidates = list(v_fixed = f_fix, v_fitted = f_vf))
  })
  names(results) <- ids
  results[!vapply(results, is.null, logical(1))]
}

#' One row per event: both candidate fits, the selection, QA flags.
#' This is the Stage-1 result table for the methods / results section.
tsm_hydraulics_table <- function(hyd_fits) {
  map_dfr(hyd_fits, function(h) {
    cand <- function(nm, what) {
      f <- h$candidates[[nm]]
      if (is.null(f)) return(NA_real_)
      switch(what, D = f$par[["D"]], alpha = f$par[["alpha"]],
             As_ratio = f$par[["As"]] / f$par[["A"]], v = f$par[["v"]], rmse = f$rmse)
    }
    tibble(
      event_id = h$event_id, model_version = h$model_version %||% NA_character_,
      conservative_source = h$conservative_source,
      logger_recovery = h$logger_recovery %||% NA_real_,
      nacl_mass_fit_g = h$nacl_mass_fit_g %||% NA_real_,
      hydraulics_borrowed_from = h$hydraulics_borrowed_from %||% NA_character_,
      hydraulic_model = h$hydraulic_model, selection_reason = h$selection_reason %||% NA_character_,
      L_m = h$L, Q_m3s = h$Q, v_input_ms = h$v_input %||% NA_real_,
      v_ms = h$v, A_m2 = h$A, D_m2s = h$D, alpha_1s = h$alpha, As_m2 = h$As,
      As_over_A = h$As / h$A, hydraulic_rmse = h$rmse,
      params_at_bound = h$params_at_bound %||% NA_character_,
      n_fit = h$n_fit %||% NA_integer_, n_pre_zeroed = h$n_pre_zeroed %||% NA_integer_,
      t_onset_s = h$t_onset %||% NA_real_, n_gap_filled = h$n_gap_filled %||% NA_integer_,
      aicc_v_fixed = h$aicc_v_fixed %||% NA_real_, aicc_v_fitted = h$aicc_v_fitted %||% NA_real_,
      d_aicc = h$d_aicc %||% NA_real_,
      vfixed_D = cand("v_fixed", "D"), vfixed_alpha = cand("v_fixed", "alpha"),
      vfixed_As_over_A = cand("v_fixed", "As_ratio"), vfixed_rmse = cand("v_fixed", "rmse"),
      vfitted_D = cand("v_fitted", "D"), vfitted_alpha = cand("v_fitted", "alpha"),
      vfitted_As_over_A = cand("v_fitted", "As_ratio"), vfitted_v = cand("v_fitted", "v"),
      vfitted_rmse = cand("v_fitted", "rmse")
    )
  })
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
#' x concentration-column combination, on the SELECTED Stage-1 hydraulics
#' (`hyd`: one element of fit_all_hydraulics()). A and v come from `hyd`.
#'
#' @param seed base seed for this job's Stage-2 LHS scan; the actual seed is
#'   job_seed(seed, "event_id|solute|conc_col") (see job_seed()).
run_one_uptake <- function(event_id, solute, conc_col, hyd, events, master_tsm,
                           nitrogen_raw, phosphate_raw, n_lhs = 250, n_cells = 40,
                           excluded_times = NULL, fill_pre_arrival = TRUE, seed = 1) {
  
  e <- events %>% filter(event_id == !!event_id)
  nut <- master_tsm %>%
    filter(event_id == !!event_id, solute == !!solute, !is.na(.data[[conc_col]])) %>%
    arrange(time_since_release_s)
  # `excluded_times`: grab timestamps (time_since_release_s) to drop before
  # fitting -- always an explicit list passed in by the caller.
  if (!is.null(excluded_times) && length(excluded_times) > 0) {
    nut <- nut %>% filter(!time_since_release_s %in% excluded_times)
  }
  base_row <- tibble(event_id = event_id, solute = solute, conc_col = conc_col,
                     hydraulic_model = hyd$hydraulic_model %||% NA_character_)
  if (nrow(nut) < 6) return(base_row %>% mutate(fit_status = "too_few_grabs"))
  
  mass_mg <- lookup_injected_mass_mg(e, nitrogen_raw, phosphate_raw, solute)
  if (!is.finite(mass_mg)) return(base_row %>% mutate(fit_status = "no_addition_record"))
  
  L <- hyd$L; Q <- hyd$Q; A <- hyd$A; v <- hyd$v; width <- hyd$width
  hydraulics <- c(D = hyd$D, alpha = hyd$alpha, As = hyd$As)
  
  obs_conc <- pmax(nut[[conc_col]], 0)
  # Nutrient grabs are always sparse, hand-taken samples, so the pre-arrival
  # gap fill applies here by default -- see fill_pre_arrival_gap().
  gf <- if (fill_pre_arrival) {
    fill_pre_arrival_gap(nut$time_since_release_s, obs_conc)
  } else {
    list(time = nut$time_since_release_s, value = obs_conc,
         kind = rep("observed", nrow(nut)), n_added = 0L, dt_used = NA_real_)
  }
  
  fitU <- tryCatch(
    fit_uptake(gf$time, gf$value, L, Q, A, hydraulics, mass_mg,
               n_lhs = n_lhs, n_cells = n_cells,
               seed = job_seed(seed, paste(event_id, solute, conc_col, sep = "|"))),
    error = function(err) NULL
  )
  if (is.null(fitU)) return(base_row %>% mutate(fit_status = "uptake_fit_failed"))
  
  part <- tryCatch(
    partition_uptake(L, Q, A, hyd$D, hyd$alpha, hyd$As,
                     fitU$par[["lambda"]], fitU$par[["lambda_s"]], mass_mg,
                     n_cells = n_cells),
    error = function(err) NULL
  )
  if (is.null(part)) {   # previously crashed on part$... when NULL
    part <- list(pct_total_uptake = NA_real_, pct_mainchannel = NA_real_,
                 pct_storagezone = NA_real_)
  }
  
  depth_main <- A / width
  depth_storage <- hyd$As / width
  Camb <- suppressWarnings(mean(nut$background_ugL, na.rm = TRUE))
  
  met <- uptake_metrics(v, depth_main, depth_storage, fitU$par[["lambda"]],
                        fitU$par[["lambda_s"]], hyd$alpha, A, hyd$As, Camb)
  
  # --- identifiability / reliability flags -----------------------------
  # A fitted coefficient at (or within 5% in log10 space of) the search
  # bound is a sign of a non-identifiable or boundary solution (Bonanno et
  # al. 2022). A trivially small total uptake makes the main-channel/
  # storage-zone SPLIT meaningless noise.
  lam_bounds <- c(-7, -1); lams_bounds <- c(-7, -1)
  near_bound <- function(val, bounds, tol = 0.05) {
    lv <- log10(val)
    (lv - bounds[1]) < tol * diff(bounds) || (bounds[2] - lv) < tol * diff(bounds)
  }
  lambda_at_bound   <- near_bound(fitU$par[["lambda"]], lam_bounds)
  lambda_s_at_bound <- near_bound(fitU$par[["lambda_s"]], lams_bounds)
  uptake_significant <- isTRUE(is.finite(part$pct_total_uptake) && part$pct_total_uptake > 2)
  
  tibble(
    event_id = event_id, solute = solute, conc_col = conc_col,
    fit_status = "ok", n_grabs = nrow(nut), mass_injected_mg = mass_mg,
    background_ugL = Camb,
    n_gap_filled_uptake = gf$n_added, gap_fill_dt_uptake_s = gf$dt_used,
    # ---- Stage-1 hydraulics used (selected model) ----
    hydraulic_model = hyd$hydraulic_model %||% NA_character_,
    conservative_source = hyd$conservative_source %||% NA_character_,
    v_input_ms = hyd$v_input %||% NA_real_, v_ms = v,
    n_gap_filled_hydraulics = hyd$n_gap_filled %||% NA_integer_,
    n_pre_zeroed_hydraulics = hyd$n_pre_zeroed %||% NA_integer_,
    hydraulic_params_at_bound = hyd$params_at_bound %||% NA_character_,
    D_m2s = hyd$D, alpha_1s = hyd$alpha, As_m2 = hyd$As, A_m2 = A,
    As_over_A = hyd$As / A, hydraulic_rmse = hyd$rmse,
    # ---- Stage-2 uptake ----
    lambda_1s = fitU$par[["lambda"]], lambda_s_1s = fitU$par[["lambda_s"]],
    uptake_rmse = fitU$rmse,
    pct_total_uptake = part$pct_total_uptake,
    pct_mainchannel = if (uptake_significant) part$pct_mainchannel else NA_real_,
    pct_storagezone = if (uptake_significant) part$pct_storagezone else NA_real_,
    Sw_mainchannel_m = if (!lambda_at_bound) met$Sw_mainchannel_m else NA_real_,
    vf_mainchannel_mps = met$vf_mainchannel_mps,
    # U = vf [m/s] * Camb [numerically ug/L == mg/m3] = mg/(m2.s); *3600 -> mg/(m2.h)
    U_mainchannel_mgm2h = met$U_mainchannel * 3600,
    vf_storagezone_mps = met$vf_storagezone_mps, U_storagezone_mgm2h = met$U_storagezone * 3600,
    Sw_total_m = if (uptake_significant) met$Sw_total_m else NA_real_,
    vf_total_mps = met$vf_total_mps,
    U_total_mgm2h = met$U_total * 3600,
    lambda_at_bound = lambda_at_bound, lambda_s_at_bound = lambda_s_at_bound,
    uptake_significant = uptake_significant
  )
}

#' Full batch pipeline across all events and solutes -- the final result
run_tsm_uptake_all <- function(data_dir = here::here("data_derived"),
                               nutrients_dir = here::here("data", "nutrients"),
                               n_lhs_hydraulics = 250, n_lhs_uptake = 250, n_cells = 40,
                               seed = 1, method = TSM_METHOD, write_outputs = FALSE) {
  # here::here() anchors paths to the project root regardless of the working
  # directory (knitting from vignettes/ sets it to vignettes/).
  
  events <- read_csv(file.path(data_dir, "events.csv"), show_col_types = FALSE)
  btc_conservative <- read_csv(file.path(data_dir, "btc_conservative.csv"), show_col_types = FALSE)
  master_tsm <- read_csv(file.path(data_dir, "master_tsm.csv"), show_col_types = FALSE)
  nitrogen_raw <- read_csv(file.path(nutrients_dir, "nutrient_addition_nitrogen_data.csv"), show_col_types = FALSE) %>%
    mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_NH4Cl_g, molar_mass_NH4Cl_nutrient)
  phosphate_raw <- read_csv(file.path(nutrients_dir, "nutrient_addition_phosphate_data.csv"), show_col_types = FALSE) %>%
    mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_PO4_g, molar_mass_P_nutrient)
  
  message("Stage 1 -- hydraulics: fitting v_fixed and v_fitted per event, selecting by AICc...")
  hyd_fits <- fit_all_hydraulics(events, btc_conservative, master_tsm, n_lhs = n_lhs_hydraulics,
                                 n_cells = n_cells, seed = seed, method = method)
  
  # Borrow hydraulics for events with no usable conservative series (see
  # HYDRAULICS_BORROWED_FROM). Carried over from the source event: D, alpha,
  # As/A, and the source's velocity CORRECTION (fitted v / input v), applied
  # to the borrowing event's own input velocity -- so the borrowing event
  # keeps its own Q and its own measured velocity, adjusted by the same
  # peak-time bias the source fit found on the same reach.
  for (eid in names(HYDRAULICS_BORROWED_FROM)) {
    src <- HYDRAULICS_BORROWED_FROM[[eid]]
    if (!(eid %in% names(hyd_fits)) && src %in% names(hyd_fits) && eid %in% events$event_id) {
      e <- events %>% filter(event_id == eid)
      if (is.na(e$discharge_Ls) || is.na(e$water_velocity_ms)) next
      s <- hyd_fits[[src]]
      Q <- e$discharge_Ls / 1000
      v <- e$water_velocity_ms * (s$v / s$v_input)
      A <- Q / v
      hyd_fits[[eid]] <- list(
        event_id = eid, L = e$reach_length_m, Q = Q, A = A, v = v,
        v_input = e$water_velocity_ms, width = e$reach_mean_width_m,
        D = s$D, alpha = s$alpha, As = s$As / s$A * A,
        rmse = NA_real_, lhs = NULL,
        hydraulic_model = s$hydraulic_model, selection_reason = paste("borrowed from", src),
        params_at_bound = s$params_at_bound, hydraulics_borrowed_from = src,
        conservative_source = "borrowed", model_version = method$model_version)
      message(sprintf("  %s: no usable conservative series; hydraulics borrowed from %s (%s)",
                      eid, src, s$hydraulic_model))
    }
  }
  
  hyd_tbl <- tsm_hydraulics_table(hyd_fits)
  message(sprintf("Stage 1 done: %d events (%d v_fitted, %d v_fixed, %d borrowed).",
                  nrow(hyd_tbl), sum(hyd_tbl$hydraulic_model == "v_fitted" & is.na(hyd_tbl$hydraulics_borrowed_from)),
                  sum(hyd_tbl$hydraulic_model == "v_fixed" & is.na(hyd_tbl$hydraulics_borrowed_from)),
                  sum(!is.na(hyd_tbl$hydraulics_borrowed_from))))
  
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
                   nitrogen_raw, phosphate_raw, n_lhs = n_lhs_uptake, n_cells = n_cells,
                   seed = seed)
  })
  
  out <- bind_rows(results) %>%
    mutate(correction = if_else(conc_col == "conc_corr_coprec_ugL",
                                "coprec_corrected_biotic_only", "raw"),
           model_version = method$model_version) %>%
    left_join(hyd_tbl %>% select(event_id, hydraulics_borrowed_from, d_aicc),
              by = "event_id") %>%
    left_join(events %>% select(event_id, stream, flag_discharge_invalid, flag_mixing_incomplete,
                                flag_no_baseline, flag_weak_signal),
              by = "event_id") %>%
    left_join(master_tsm %>% distinct(event_id, flag_broken_pairing), by = "event_id") %>%
    relocate(stream, .after = event_id)
  
  attr(out, "hydraulic_fits")   <- hyd_fits
  attr(out, "hydraulics_table") <- hyd_tbl
  attr(out, "method")           <- method
  
  if (write_outputs) {
    write_csv(out, file.path(data_dir, "tsm_uptake_summary.csv"))
    write_csv(hyd_tbl, file.path(data_dir, "tsm_hydraulics_summary.csv"))
    message("written: tsm_uptake_summary.csv, tsm_hydraulics_summary.csv in ", data_dir)
  }
  out
}