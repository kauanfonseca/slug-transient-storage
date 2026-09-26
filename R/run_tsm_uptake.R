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
#   Before Stage 2 fitting, each event/solute's nutrient grabs go through:
#     a. hard-coded exclusions (NUTRIENT_EXCLUDED_TIMES / bad samples,
#        NUTRIENT_VALUE_THRESHOLD / spike values) -- see their comments;
#     b. a general clock-offset search (find_nutrient_time_shift(),
#        TSM_METHOD$nutrient_time_shift) that aligns the grab times against
#        the SELECTED Stage-1 conservative shape, correcting a nutrient-vs-
#        logger clock desync found in nearly every event (2026-09-25, see
#        claude/tsm_uptake_pipeline_status.md). The shift itself is not
#        hard-coded -- only the search is -- and is reported per event x
#        solute as nutrient_time_shift_s in the output.
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
  bound_tol        = 0.02,         # "at bound" = within 2% of log10 range
  nutrient_time_shift = TRUE,      # search + apply a per-event/solute clock-offset correction before Stage 2 (see find_nutrient_time_shift())
  nutrient_shift_range_s = seq(-900, 900, by = 15),
  storage_profile      = TRUE,     # profile lambda_s -> ranges for the channel/storage split (see profile_storage_uptake())
  storage_profile_tol  = 1.05,     # keep (lambda, lambda_s) pairs with sqrt-RMSE <= tol x best
  storage_profile_grid = seq(-7, -1, by = 0.25),  # log10(lambda_s) grid, same bounds as fit_uptake()
  storage_detect_pct   = 1         # storage uptake "detected" when the profile's lower bound of pct_storagezone exceeds this
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
#' 2026-09-26 (Kauan): SR_20231011_single no longer borrows -- the list is
#' now empty and the mechanism is kept for future use. The borrowed
#' parameters (D 0.006, As/A 0.39, from the higher-flow 09 Oct day) predict
#' a NaCl peak of ~630 mg/L at 14 min on 11 Oct, against the day's own 28
#' probe grabs peaking at ~150 mg/L at 17 min (RMSE 4.06 on those grabs).
#' Fitted on the grabs themselves (probe_grab, like any event with
#' discharge_source == "probe"), the AICc rule selects v_fitted (D 0.043,
#' alpha 3.2e-3, As/A 0.61, nothing at a bound, RMSE 0.31). The grabs are
#' dense (20 s) from 590 to 1250 s and then only one at 2600 s, so the tail
#' is loosely constrained. See diag_sr1011.R.
HYDRAULICS_BORROWED_FROM <- setNames(character(0), character(0))

#' Events that share ANOTHER event's Stage-1 hydraulics outright -- same
#' point, same day, one TSM estimate for both slugs -- as opposed to
#' HYDRAULICS_BORROWED_FROM above, which corrects for a genuinely different
#' day's flow via a velocity ratio and keeps the borrowing event's own Q.
#'
#' RA_20230906_N (2026-09-25, Kauan): reviewed the logger_rescaled fits for
#' N and P side by side -- P's came out acceptable, N's independent fit did
#' not (and applying P's parameters to N's own curve, checked separately,
#' still didn't rise cleanly). Since N and P are the same reach and the same
#' day, use P's Stage-1 result -- Q, A, v, D, alpha, As, all of it -- as the
#' hydraulics for N's Stage-2 uptake fit too, unmodified. N's own
#' water_velocity_ms/discharge_Ls are NOT used for anything past this point;
#' N is not independently fit or velocity-corrected, it shares P's number.
HYDRAULICS_SHARED_WITH <- c(RA_20230906_N = "RA_20230906_P")

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

#' Time-to-peak (s) used for the Stage-1 input velocity v_input = L / t_peak,
#' overriding events$water_velocity_ms. v_input is the fixed v of the
#' "v_fixed" candidate and the centre of the "v_fitted" search (v_input/2 ..
#' v_input*2), so a wrong t_peak makes both candidates fail.
#'
#' RA_20231005_downstream (2026-09-26, diag_ra1005.R): events$t_peak_s =
#' 6780 s comes from t_peak_s_probe in 02_integration.Rmd, the last grab the
#' field sheet labelled curve_limb == "rising". The NaCl grabs peak at 5220 s
#' (7.52 mg/L); 6420 s is already 6.13 and 6780 s is 4.75, on the falling
#' limb (there is a 20-min sampling gap 5220-6420 s, so the true peak lies in
#' 5220-6000 s). Every other event's labelled t_peak is within a few minutes
#' of its grab maximum. With 6780 s: v_fixed ends with alpha and As/A at
#' their bounds (RMSE 0.43) and v_fitted runs into v_input*2 (RMSE 0.11,
#' dAICc 120, rejected for being at the bound). With 5220 s: v_fitted is
#' selected, nothing at a bound, RMSE 0.064. The pre-arrival gap fill
#' (zeros 0-2880 s) is kept: the logger, useless for amplitude, still shows
#' the tracer front arriving at ~50 min, consistent with the first grab
#' (3120 s, 1.98 mg/L).
V_INPUT_TPEAK_OVERRIDE <- c(RA_20231005_downstream = 5220)

#' Stage-1 candidate forced for an event, overriding the AICc rule.
#'
#' RA_20231005_downstream (Kauan, 2026-09-26): with t_peak = 5220 s the AICc
#' rule picks v_fitted (RMSE 0.064 vs 0.181, dAICc 91), but that fit reaches
#' its RMSE by making the front a near-step exactly at the boundary between
#' the synthetic pre-arrival zeros (last at 2880 s) and the first grab
#' (3120 s, already 26% of peak), with D = 5.3e-4 m2/s: 10x below the
#' lowest D fitted on any logger curve here (0.005-0.11 m2/s). That shape is
#' set by where the gap fill ends, not by the data. The v_fixed candidate
#' (v = L/5220 s, D = 0.10, alpha = 1.0e-4, As/A = 0.10, nothing at a bound)
#' is smooth, matches the peak and the tail, and misses only the first
#' three grabs of the rising limb. Stage 2 barely changes: NH4-N 50 vs 54%
#' total uptake, SRP 97.5 vs 97%.
HYDRAULIC_MODEL_OVERRIDE <- c(RA_20231005_downstream = "v_fixed")

#' Hard-coded per-event/solute nutrient-grab exclusions, by their raw
#' `time_since_release_s` (BEFORE the time-shift correction below -- see
#' find_nutrient_time_shift()). All identified by Kauan, 2026-09-25, from
#' claude/tsm_uptake_pipeline_status.md's nutrient-fit diagnosis and the
#' resulting time-shift-corrected plots:
#'
#'  - CB_20230907_single SRP: 5 grabs (t=2280/2520/2820/3240/3840) jump to
#'    212/615/630/39/152 ug/L over ~26 min while the paired NaCl is already
#'    declining (15.5 -> 0 mg/L) -- a 20x spike-and-crash with no plausible
#'    transport/reaction explanation; a clock shift cannot produce this
#'    shape (it only translates the curve, it can't invert its slope), so
#'    these are treated as bad samples, not a timing artifact.
#'  - RA_20230906_N NH4-N: first 2 grabs (t=2460, 2700), taken while the
#'    paired NaCl was still exactly 0 -- i.e. before the co-injected tracer
#'    had arrived at all, the same pre-arrival-background issue the overall
#'    diagnosis was built on.
#'  - RA_20230906_P SRP: first grab (t=720), isolated ~38 min before the
#'    next one (t=3000) and, like the NH4-N pair above, far pre-arrival
#'    (paired NaCl = 0).
#'  - RA_20231005_downstream NH4-N: 2nd and 5th grabs by time order
#'    (t=3480, t=4260) -- this event's Stage-1 hydraulics are already the
#'    weakest/least identifiable in the set (see handoff.md), so these are
#'    a smaller, more surgical trim rather than a broad rule.
NUTRIENT_EXCLUDED_TIMES <- list(
  CB_20230907_single      = list(SRP     = c(2280, 2520, 2820, 3240, 3840)),
  RA_20230906_N           = list(`NH4-N` = c(2460, 2700)),
  RA_20230906_P           = list(SRP     = c(720)),
  RA_20231005_downstream  = list(`NH4-N` = c(3480, 4260))
)

#' Per-event/solute concentration-value thresholds: any grab (post
#' background-correction, i.e. on the fitted `conc_col`) above the given
#' value (ug/L) is dropped before Stage 2. Kauan, 2026-09-25: specified as
#' a value rule rather than specific points for these three -- in each case
#' the excluded grabs are isolated single-sample spikes sitting on an
#' otherwise low, flat background (the same spike pattern as
#' CB_20230907_single SRP above), not part of the main breakthrough curve;
#' the corresponding NH4-N series at RS/SR, which has a normal, smooth
#' arrival peak that legitimately exceeds these values, was NOT given a
#' threshold.
NUTRIENT_VALUE_THRESHOLD <- list(
  RA_20231005_downstream = list(SRP = 25),
  RS_20231011_single     = list(SRP = 50),
  SR_20231011_single     = list(SRP = 60)
)

#' Drop NUTRIENT_EXCLUDED_TIMES / NUTRIENT_VALUE_THRESHOLD rows from `nut`
#' (a master_tsm slice for one event x solute x conc_col), returning the
#' filtered tibble with attr(., "n_dropped_manual") set.
apply_nutrient_exclusions <- function(event_id, solute, conc_col, nut) {
  n_dropped <- 0L
  bad_t <- NUTRIENT_EXCLUDED_TIMES[[event_id]][[solute]]
  if (!is.null(bad_t)) {
    n_dropped <- n_dropped + sum(nut$time_since_release_s %in% bad_t)
    nut <- nut %>% filter(!time_since_release_s %in% bad_t)
  }
  thr <- NUTRIENT_VALUE_THRESHOLD[[event_id]][[solute]]
  if (!is.null(thr)) {
    over <- nut[[conc_col]] > thr
    n_dropped <- n_dropped + sum(over, na.rm = TRUE)
    nut <- nut %>% filter(!over)
  }
  attr(nut, "n_dropped_manual") <- n_dropped
  nut
}

#' Per-event/solute clock-offset correction for nutrient grabs (Kauan,
#' 2026-09-25): the nutrient-grab clock and the logger clock were not
#' synchronized, so `time_since_release_s` for the hand-taken nutrient
#' samples can be off by a few minutes. Found by grid search over
#' `TSM_METHOD$nutrient_shift_range_s`: for each candidate shift, evaluate
#' the SELECTED Stage-1 hydraulics' pure conservative shape (lambda =
#' lambda_s = 0, at the event's own `nacl_mass_g` -- reaction only reshapes
#' the tail a little and barely touches the rising edge, so this is a clean
#' shape reference for a pure timing offset) at `times = obs_t - shift`,
#' find the best-fit linear amplitude in closed form (this is a SHAPE
#' match, not a mass-balance one), and score by SSE against the observed
#' (already exclusion-filtered) grabs. Across nearly every event this
#' finds a consistent -180s to -345s shift (58-97% SSE reduction versus no
#' shift), i.e. the nutrient clock was running fast / under-counting
#' elapsed time -- see claude/tsm_uptake_pipeline_status.md for the full
#' validation (including why a single mis-timed cluster, CB_20230907_single
#' SRP, needed its own point exclusion above instead: a shift only
#' translates the curve, it can't explain that shape).
#'
#' NOTE on sign: this fits model(obs_t - shift) ~ obs_c, so a grab's TRUE
#' time is (obs_t - shift) -- the correction applied downstream is
#' `time_shifted <- time_since_release_s - shift_s`, not `+`.
find_nutrient_time_shift <- function(conc_col, hyd, mass_nacl_mg, nut,
                                     shift_range = TSM_METHOD$nutrient_shift_range_s) {
  obs_t <- nut$time_since_release_s
  obs_c <- pmax(nut[[conc_col]], 0)
  sse_for_shift <- function(shift) {
    model_c <- simulate_tsm(L = hyd$L, Q = hyd$Q, A = hyd$A, D = hyd$D, alpha = hyd$alpha,
                            As = hyd$As, lambda = 0, lambda_s = 0, mass = mass_nacl_mg,
                            times = obs_t - shift, n_cells = 40)$C
    if (any(!is.finite(model_c)) || sum(model_c^2) == 0) return(Inf)
    k <- sum(obs_c * model_c) / sum(model_c^2)   # closed-form best-fit amplitude
    if (!is.finite(k) || k <= 0) return(Inf)
    sum((obs_c - k * model_c)^2)
  }
  sses <- vapply(shift_range, sse_for_shift, numeric(1))
  if (all(!is.finite(sses))) return(list(shift = 0, sse_best = NA_real_, sse_zero = NA_real_))
  best <- shift_range[which.min(sses)]
  list(shift = best, sse_best = min(sses), sse_zero = sses[which(shift_range == 0)])
}

#' Profile of the Stage-2 cost over lambda_s -> RANGES for the channel vs
#' storage-zone split (2026-09-25, Kauan: the chapter needs zone-separated
#' uptake to compare with fish excretion).
#'
#' Why: As and alpha are well identified from the logger NaCl curve (Stage 1),
#' but lambda_s (reaction INSIDE the storage zone) can only come from the
#' sparse nutrient grabs. In most events the cost surface is flat from
#' lambda_s = 1e-7 up to ~1e-4 and rises after that, so the point estimate
#' lands on the lower bound: storage uptake is small, but "how small" is not
#' resolved. A synthetic test on the real hydraulics and grab times
#' (diag_storage_synthetic.R) recovered a true storage share of 45% to within
#' ~5-7 points, but a true share of 15-20% came back as 0-27%: with the
#' current grab design, storage shares below ~20-30% of total uptake are
#' indistinguishable from zero. So the split is reported as a min-max range
#' over every (lambda, lambda_s) pair that fits within `tol` of the best
#' sqrt-RMSE, alongside the point estimate.
#'
#' For each lambda_s on the grid, lambda is re-optimised (1-D), and for every
#' pair kept (plus the fitted pair itself) the Runkel (2007) partition and
#' the uptake metrics are recomputed.
#'
#' Interval ends are exact, not grid points (2026-09-26): on each side of the
#' accepted grid points, uniroot() finds the log10(lambda_s) where the
#' profile crosses tol x best, and the partition/metrics are evaluated
#' there. If the lowest grid point (1e-7) is still accepted, the lower end
#' is lambda_s = 0 itself (storage share exactly 0, "not detected"); if the
#' highest (1e-1) is still accepted, the upper end stays at the search limit
#' and `upper_at_limit` is TRUE.
profile_storage_uptake <- function(gf, L, Q, A, D, alpha, As, mass, fit_par, fit_rmse,
                                   depth_main, depth_storage, v, Camb, n_cells = 40,
                                   method = TSM_METHOD) {
  prof_at <- function(lls) {   # lambda_s = 10^lls (lls = -Inf -> 0); lambda re-optimised
    ls <- if (is.finite(lls)) 10^lls else 0
    f <- function(ll) {
      sim <- simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = alpha, As = As,
                          lambda = 10^ll, lambda_s = ls, mass = mass,
                          times = gf$time, n_cells = n_cells)$C
      sqrt_rmse(sim, gf$value)
    }
    o <- optimize(f, c(-7, -1))
    c(lambda = 10^o$minimum, lambda_s = ls, rmse = o$objective, lls = lls)
  }
  grid <- method$storage_profile_grid
  prof <- as.data.frame(do.call(rbind, lapply(grid, prof_at)))
  best <- min(c(prof$rmse, fit_rmse))
  thr  <- method$storage_profile_tol * best
  ok   <- which(prof$rmse <= thr)
  if (length(ok) == 0) ok <- which.min(prof$rmse)
  cross <- function(i_in, i_out) {   # exact crossing between an accepted and a rejected grid point
    r <- uniroot(function(x) prof_at(x)[["rmse"]] - thr, sort(c(grid[i_in], grid[i_out])), tol = 1e-3)
    prof_at(r$root)
  }
  lo <- if (min(ok) > 1) cross(min(ok), min(ok) - 1) else prof_at(-Inf)
  upper_at_limit <- max(ok) == length(grid)
  hi <- if (!upper_at_limit) cross(max(ok), max(ok) + 1) else unlist(prof[max(ok), ])
  keep <- rbind(as.data.frame(t(lo)), prof[ok, ], as.data.frame(t(hi)),
                data.frame(lambda = fit_par[["lambda"]], lambda_s = fit_par[["lambda_s"]],
                           rmse = fit_rmse, lls = log10(fit_par[["lambda_s"]])))

  res <- lapply(seq_len(nrow(keep)), function(k) {
    p <- tryCatch(partition_uptake(L, Q, A, D, alpha, As, keep$lambda[k], keep$lambda_s[k],
                                   mass, n_cells = n_cells), error = function(e) NULL)
    if (is.null(p)) return(NULL)
    m <- uptake_metrics(v, depth_main, depth_storage, keep$lambda[k], keep$lambda_s[k],
                        alpha, A, As, Camb)
    c(pct_total = p$pct_total_uptake, pct_main = p$pct_mainchannel, pct_storage = p$pct_storagezone,
      U_main = m$U_mainchannel * 3600, U_storage = m$U_storagezone * 3600, U_total = m$U_total * 3600)
  })
  res <- as.data.frame(do.call(rbind, res[!vapply(res, is.null, logical(1))]))
  rng <- function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  list(n_pairs = nrow(res), lambda_s_min = lo[["lambda_s"]], lambda_s_max = hi[["lambda_s"]],
       upper_at_limit = upper_at_limit,
       pct_total = rng(res$pct_total), pct_main = rng(res$pct_main),
       pct_storage = rng(res$pct_storage), U_main = rng(res$U_main),
       U_storage = rng(res$U_storage), U_total = rng(res$U_total))
}

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
           !(event_id %in% names(HYDRAULICS_BORROWED_FROM)),
           !(event_id %in% names(HYDRAULICS_SHARED_WITH))) %>%
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
    v_input <- if (eid %in% names(V_INPUT_TPEAK_OVERRIDE)) L / V_INPUT_TPEAK_OVERRIDE[[eid]] else e$water_velocity_ms
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
    if (eid %in% names(HYDRAULIC_MODEL_OVERRIDE)) {   # documented per-event choice, see the constant
      forced <- HYDRAULIC_MODEL_OVERRIDE[[eid]]
      f_forced <- if (forced == "v_fitted") f_vf else f_fix
      if (!is.null(f_forced)) {
        reason <- sprintf("override: %s (AICc rule: %s, %s)", forced,
                          if (choose_vf) "v_fitted" else "v_fixed", reason)
        choose_vf <- forced == "v_fitted"; sel <- f_forced
      }
    }

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
  # fitting -- always an explicit list passed in by the caller. Combined
  # with the hard-coded NUTRIENT_EXCLUDED_TIMES / NUTRIENT_VALUE_THRESHOLD
  # (see their comments above run_one_uptake in this file).
  if (!is.null(excluded_times) && length(excluded_times) > 0) {
    nut <- nut %>% filter(!time_since_release_s %in% excluded_times)
  }
  nut <- apply_nutrient_exclusions(event_id, solute, conc_col, nut)
  n_dropped_manual <- attr(nut, "n_dropped_manual") %||% 0L

  base_row <- tibble(event_id = event_id, solute = solute, conc_col = conc_col,
                     hydraulic_model = hyd$hydraulic_model %||% NA_character_)
  if (nrow(nut) < 6) return(base_row %>% mutate(fit_status = "too_few_grabs"))

  mass_mg <- lookup_injected_mass_mg(e, nitrogen_raw, phosphate_raw, solute)
  if (!is.finite(mass_mg)) return(base_row %>% mutate(fit_status = "no_addition_record"))

  L <- hyd$L; Q <- hyd$Q; A <- hyd$A; v <- hyd$v; width <- hyd$width
  hydraulics <- c(D = hyd$D, alpha = hyd$alpha, As = hyd$As)

  # Clock-offset correction (see find_nutrient_time_shift() above): search
  # once the manual exclusions are already applied, then shift this
  # event/solute's grab times before anything downstream sees them. Skipped
  # (shift forced to 0) when TSM_METHOD$nutrient_time_shift is off, or when
  # the event has no nacl_mass_g to build a conservative reference curve
  # from (e.g. a borrowed/shared hydraulics event with a missing addition
  # record -- rare, falls back to uncorrected timing rather than failing).
  # No shift when Stage 1 was fit on the hand-held probe grabs themselves
  # (conservative_source "probe_grab", currently RA_20231005_downstream):
  # those NaCl readings share the nutrient samples' own timestamps, so there
  # is no second clock to correct -- a nonzero "shift" there only absorbs
  # differences in curve shape (tested 2026-09-26: +330 s for NH4-N and the
  # -900 s bound for the near-flat SRP series).
  shift_s <- 0
  if (isTRUE(TSM_METHOD$nutrient_time_shift) && is.finite(e$nacl_mass_g) &&
      !identical(hyd$conservative_source, "probe_grab")) {
    ts <- find_nutrient_time_shift(conc_col, hyd, e$nacl_mass_g * 1000, nut)
    shift_s <- ts$shift
  }
  nut <- nut %>% mutate(time_since_release_s = time_since_release_s - shift_s) %>%
    filter(time_since_release_s >= 0) %>% arrange(time_since_release_s)
  if (nrow(nut) < 6) return(base_row %>% mutate(fit_status = "too_few_grabs"))

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

  # ranges for the channel/storage split (see profile_storage_uptake());
  # withheld (NA) when total uptake isn't significant, same as the point split
  na2 <- c(NA_real_, NA_real_)
  pr <- list(n_pairs = NA_integer_, lambda_s_min = NA_real_, lambda_s_max = NA_real_,
             upper_at_limit = NA, pct_total = na2, pct_main = na2,
             pct_storage = na2, U_main = na2, U_storage = na2, U_total = na2)
  if (isTRUE(TSM_METHOD$storage_profile) && uptake_significant) {
    pr <- tryCatch(
      profile_storage_uptake(gf, L, Q, A, hyd$D, hyd$alpha, hyd$As, mass_mg, fitU$par, fitU$rmse,
                             depth_main, depth_storage, v, Camb, n_cells = n_cells),
      error = function(err) pr)
    # with Stage-1 alpha/As at a bound (RA_20231005_downstream: alpha ~ 0,
    # As/A = 10) the storage zone barely exchanges, any lambda_s fits, and
    # the storage flux range is meaningless (U_storagezone ran to 1e4
    # mg/m2/h) -- keep the total, withhold the split
    pab <- hyd$params_at_bound %||% ""
    if (!is.na(pab) && nzchar(pab)) {
      pr$pct_main <- na2; pr$pct_storage <- na2; pr$U_main <- na2; pr$U_storage <- na2
    }
  }
  storage_detected <- isTRUE(pr$pct_storage[1] > TSM_METHOD$storage_detect_pct)

  tibble(
    event_id = event_id, solute = solute, conc_col = conc_col,
    fit_status = "ok", n_grabs = nrow(nut), mass_injected_mg = mass_mg,
    background_ugL = Camb,
    nutrient_time_shift_s = shift_s, n_dropped_manual = n_dropped_manual,
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
    uptake_significant = uptake_significant,
    # ---- profile ranges for the channel/storage split ----
    # every (lambda, lambda_s) pair within storage_profile_tol of the best
    # sqrt-RMSE; use these (not the point split) when the zone split matters
    n_profile_pairs = pr$n_pairs, lambda_s_min_profile_1s = pr$lambda_s_min,
    lambda_s_max_profile_1s = pr$lambda_s_max, storage_upper_at_limit = pr$upper_at_limit,
    pct_total_uptake_min = pr$pct_total[1], pct_total_uptake_max = pr$pct_total[2],
    pct_mainchannel_min = pr$pct_main[1], pct_mainchannel_max = pr$pct_main[2],
    pct_storagezone_min = pr$pct_storage[1], pct_storagezone_max = pr$pct_storage[2],
    U_mainchannel_mgm2h_min = pr$U_main[1], U_mainchannel_mgm2h_max = pr$U_main[2],
    U_storagezone_mgm2h_min = pr$U_storage[1], U_storagezone_mgm2h_max = pr$U_storage[2],
    U_total_mgm2h_min = pr$U_total[1], U_total_mgm2h_max = pr$U_total[2],
    storage_uptake_detected = storage_detected
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
  
  # Share hydraulics outright for events on HYDRAULICS_SHARED_WITH (same
  # point, same day -- see the constant's comment above). Unlike the
  # borrowing loop above, nothing is recomputed from the target event's own
  # Q/velocity: Q, A, v, D, alpha, As are copied unchanged from the source.
  # v_input is kept from the target's own events.csv row for reporting only
  # (so the gap between what was measured for N and what was actually used
  # stays visible in the output) -- it plays no role in the fit.
  for (eid in names(HYDRAULICS_SHARED_WITH)) {
    src <- HYDRAULICS_SHARED_WITH[[eid]]
    if (!(eid %in% names(hyd_fits)) && src %in% names(hyd_fits) && eid %in% events$event_id) {
      e <- events %>% filter(event_id == eid)
      s <- hyd_fits[[src]]
      hyd_fits[[eid]] <- list(
        event_id = eid, L = s$L, Q = s$Q, A = s$A, v = s$v,
        v_input = e$water_velocity_ms, width = s$width,
        D = s$D, alpha = s$alpha, As = s$As,
        rmse = NA_real_, lhs = NULL,
        hydraulic_model = s$hydraulic_model,
        selection_reason = paste("shared with", src, "(same point/day, unmodified)"),
        params_at_bound = s$params_at_bound, hydraulics_borrowed_from = src,
        conservative_source = "shared", model_version = method$model_version)
      message(sprintf("  %s: shares %s's hydraulics outright (same point/day) (%s)",
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