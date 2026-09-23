################################################################################
# tsm_calibrate.R
#
# Two-stage calibration of the transient-storage model (tsm_model.R) without
# needing OTIS-P, following the standard two-step logic of Runkel (2007):
#
#   Stage 1 (hydraulics): fit D, alpha, As against the CONSERVATIVE tracer
#            (NaCl) breakthrough curve, which is dense (logger, seconds
#            resolution) and well constrained. Q is taken as already known
#            from the independent dilution-gauging discharge estimate in
#            `events.csv`.
#
#            A (= Q/v) can EITHER be fixed from `events$water_velocity_ms`
#            (fit_velocity = FALSE, the original behaviour) OR fitted as a
#            fourth parameter (fit_velocity = TRUE), as OTIS-P does. Why the
#            option exists (2026-09): `water_velocity_ms` is reach length /
#            time-to-peak of the tracer curve itself. In a transient-storage
#            model the peak arrives LATER than L/v, because storage holds the
#            tracer back, so a peak-time velocity underestimates the
#            main-channel velocity and locks the model's timing. On
#            CB_20230907_single, freeing v (Q fixed) roughly halved the fit
#            error and moved v by +14% (diag_tsm_shape.R, variants V1 vs V2).
#            Q stays fixed in both modes: it is measured independently and
#            fixing it keeps mass balance anchored (Bonanno et al. 2022).
#
#   Stage 2 (uptake): with D, alpha, As (and A) FIXED at their Stage-1
#            values, fit lambda (main-channel) and lambda_s (storage-zone)
#            first-order uptake coefficients against the sparse nutrient
#            grab series.
#
# Both stages use the same two-step search: a global Latin Hypercube scan
# (cheap, many model runs, no gradient) to avoid local minima and to obtain a
# lightweight identifiability diagnostic (parameter vs RMSE, in the spirit of
# Bonanno et al. 2022's global identifiability analysis, though nowhere near
# as exhaustive as their 115,000-run iterative DYNIA procedure -- this is a
# fast, practical compromise for routinely processing many streams), followed
# by local Nelder-Mead refinement from the best LHS point.
#
# The cost function compares sqrt(simulated) to sqrt(observed) concentration
# (rather than raw concentration). This down-weights the sharp peak relative
# to the recession/tail, which is where most of the information about alpha,
# As (and lambda_s) actually lives (Runkel 2007; Bonanno et al. 2022) -- a
# plain linear RMSE is dominated by getting the peak height right and pays
# almost no attention to the tail.
#
# Reported parameters are CLAMPED to the search bounds (2026-09 fix).
# Previously the cost function clamped internally but the unclamped
# Nelder-Mead vector was returned, so a parameter pinned at a bound could be
# reported outside it (e.g. D = 8.7e-5 with a 1e-4 floor).
################################################################################

#' Derive a stable, portable per-job integer seed from a base seed and a
#' string key (e.g. an event_id, or "event_id|solute|conc_col").
#'
#' Why this exists: `fit_hydraulics()`/`fit_uptake()` draw their global LHS
#' scan from R's ordinary RNG stream. The batch driver used to call
#' `set.seed(seed)` exactly ONCE, before looping over all events (Stage 1)
#' or all event x solute x correction jobs (Stage 2) -- so every job's LHS
#' draw depended on wherever the RNG stream happened to be after every job
#' that ran before it. That is fragile in a way that bit us in practice: two
#' runs of the identical code on the identical data, on two different
#' machines (this was first run in the cloud sandbox on Linux; re-knit later
#' on Kauan's Mac), produced materially different fitted lambda/lambda_s for
#' several events (e.g. RS_20231011_single NH4-N: pct_total_uptake 0.01% ->
#' 26.57%, uptake_significant FALSE -> TRUE) even with the same `seed = 1`
#' argument -- almost certainly because `sample()` (used inside
#' `lhs_sample()` to break correlation across parameters) is only guaranteed
#' to reproduce identically for a given seed on the SAME R version/platform
#' (R's default `sample.kind` changed in R 3.6.0, and even the same version
#' can differ in exactly how many random draws upstream code consumed).
#' Giving every job its OWN seed, derived deterministically from a stable
#' string key rather than from call order, makes each job self-contained:
#' its result no longer depends on how many other jobs ran before it, what
#' order they ran in, or (mostly) which R/platform build ran them, only on
#' its own identity and the shared base `seed`. See `fit_all_hydraulics()`
#' and `run_tsm_uptake_all()`/`run_one_uptake()` for where this is used.
#'
#' Deliberately does not use any of R's own string-hashing (`rlang::hash()`,
#' environment addresses, etc.) since those are not guaranteed portable
#' across R versions/platforms either -- `utf8ToInt()` plus plain integer
#' arithmetic is base R and gives the same result everywhere.
job_seed <- function(seed, key) {
  # Deliberately done in double precision (not integer): R integers are
  # 32-bit, and `acc * 31L` overflows well before the modulo below can rein
  # it back in. Doubles are exact for integers up to 2^53, and `h` is kept
  # under 2147483647 after every step, so `h * 31` (~6.7e10) never leaves
  # that exact range.
  codes <- utf8ToInt(key)
  h <- 0
  for (code in codes) h <- (h * 31 + code) %% 2147483647
  as.integer((h + seed) %% 2147483647)
}

#' Simple Latin Hypercube sample
#'
#' @param bounds named list of c(min, max), one per parameter (natural units;
#'   the caller decides whether to work in log space by passing log10 bounds)
#' @param n number of draws
#' @return data.frame, n rows, one column per parameter
lhs_sample <- function(bounds, n) {
  p <- length(bounds)
  m <- matrix(NA_real_, nrow = n, ncol = p)
  colnames(m) <- names(bounds)
  for (j in seq_len(p)) {
    edges <- seq(0, 1, length.out = n + 1)
    u <- runif(n, edges[-(n + 1)], edges[-1])
    u <- sample(u)  # break correlation across parameters
    rng <- bounds[[j]]
    m[, j] <- rng[1] + u * (rng[2] - rng[1])
  }
  as.data.frame(m)
}

#' Root-mean-square error on sqrt-transformed concentrations
sqrt_rmse <- function(sim, obs) {
  sim <- pmax(sim, 0)
  obs <- pmax(obs, 0)
  sqrt(mean((sqrt(sim) - sqrt(obs))^2, na.rm = TRUE))
}

#' Clamp each element of a log10 parameter vector to its bounds
#' (bounds: list of c(min, max), same order as `x`)
clamp_to_bounds <- function(x, bounds) {
  unname(mapply(function(val, b) min(max(val, b[1]), b[2]), x, bounds))
}

#' Thin a time series to at most `n_max` points for faster fitting, always
#' keeping the first and last point. Used only inside the optimisation loop;
#' the final reported fit is re-evaluated / plotted against the full series.
thin_series <- function(time, value, n_max = 200) {
  n <- length(time)
  if (n <= n_max) return(data.frame(time = time, value = value))
  idx <- unique(round(seq(1, n, length.out = n_max)))
  data.frame(time = time[idx], value = value[idx])
}

#' Zero the pre-arrival baseline drift of a conservative-tracer curve.
#'
#' Why (2026-09): on logger curves the background-corrected signal often
#' creeps up BEFORE the tracer can have arrived (baseline drift -- e.g.
#' CB_20230907_single climbs 0 -> ~0.7 mg/L over the first ~1400 s, 47 of
#' 153 points). Physically that stretch is "nothing has arrived yet" (the
#' same reasoning as fill_pre_arrival_gap() for grab series), but left in,
#' those points pull the fit toward a wide, early, symmetric curve: with the
#' exact solver and v free, zeroing them moved CB_20230907_single from
#' D = 0.13, As/A = 0.16 (peak-fit error 0.85 mg/L) to D = 0.05,
#' As/A = 0.28 (0.63 mg/L). This is a STOPGAP inside the fit; the proper fix
#' is a drifting-baseline correction upstream, in the integration step.
#'
#' Onset = the last observation before the peak that is below
#' onset_frac * peak. Everything strictly before it is set to 0.
#'
#' @return list(value, n_zeroed, t_onset)
zero_pre_arrival <- function(time, value, onset_frac = 0.05) {
  if (length(value) < 3 || all(!is.finite(value))) {
    return(list(value = value, n_zeroed = 0L, t_onset = NA_real_))
  }
  i_pk  <- which.max(value)
  below <- which(value[seq_len(i_pk)] < onset_frac * value[i_pk])
  if (length(below) == 0) return(list(value = value, n_zeroed = 0L, t_onset = NA_real_))
  t_onset <- time[max(below)]
  pre <- time < t_onset
  value[pre] <- 0
  list(value = value, n_zeroed = sum(pre), t_onset = t_onset)
}

#' Fill the "pre-arrival" gap in a sparse grab/probe series with synthetic
#' background points, at the same sampling interval used later in that curve.
#'
#' Field rationale (Kauan): monitoring starts at the moment of tracer
#' release, but a hand-held-probe or nutrient grab is only WRITTEN DOWN once
#' the signal is seen to change -- if the reading hasn't moved, there is
#' nothing new to log. So every grab-based series (hand-held-probe NaCl
#' grabs, nutrient grabs -- NOT the continuously-auto-logging conductivity
#' logger, which records regardless of whether the value moved) has a real,
#' physically-justified stretch between t=0 (release) and the first logged
#' point during which the true background-corrected concentration was ~0
#' throughout. That stretch is not missing data -- it is unlogged constancy.
#' Leaving it empty starves the model of exactly the "nothing happened yet"
#' information the dense logger curves get for free, which was making sparse
#' grab-based fits needlessly sensitive to wherever the first grab happened
#' to land (this is what motivated the fix -- see the comment in
#' simulate_tsm(), tsm_model.R, for the separate but related time-alignment
#' bug found via the same review).
#'
#' @param time,value observed time (s, sorted ascending, deduplicated) and
#'   concentration (already background-corrected, so "no signal yet" == 0)
#' @param dt sampling interval (s) used for the synthetic fill points; if
#'   NULL, uses the median spacing of the observed series itself -- i.e.
#'   "this curve's own sampling frequency", per Kauan's instruction, rather
#'   than a single interval assumed across all events
#' @param min_gap_factor only fill when the first observed time is more than
#'   this many multiples of dt after t=0, so a series that already starts
#'   right at release doesn't get one redundant synthetic point
#' @return list(time, value, kind, n_added, dt_used). `kind` is a parallel
#'   character vector ("gap_fill" / "observed") for plotting/auditing which
#'   points are real vs. assumed; n_added/dt_used record what was assumed so
#'   callers can report it rather than silently changing the series.
fill_pre_arrival_gap <- function(time, value, dt = NULL, min_gap_factor = 1.5) {
  kind_obs <- rep("observed", length(time))
  if (length(time) < 2 || is.na(time[1]) || time[1] <= 0) {
    return(list(time = time, value = value, kind = kind_obs, n_added = 0L, dt_used = NA_real_))
  }
  if (is.null(dt)) dt <- stats::median(diff(time))
  if (!is.finite(dt) || dt <= 0 || time[1] < min_gap_factor * dt) {
    return(list(time = time, value = value, kind = kind_obs, n_added = 0L, dt_used = dt))
  }
  fill_times <- seq(0, time[1] - dt, by = dt)
  list(time = c(fill_times, time), value = c(rep(0, length(fill_times)), value),
       kind = c(rep("gap_fill", length(fill_times)), kind_obs),
       n_added = length(fill_times), dt_used = dt)
}

#' Stage 1 -- fit D, alpha, As (and optionally v, hence A) to a
#' conservative-tracer breakthrough curve
#'
#' @param obs_time,obs_conc observed time (s) and background-corrected
#'   conservative-tracer concentration (any consistent unit, e.g. mg/L)
#' @param L,Q reach length (m), discharge (m3/s). Q is always fixed.
#' @param A main-channel area (m2). With fit_velocity = FALSE it is used as
#'   given; with fit_velocity = TRUE it is only the CENTRE of the velocity
#'   search (v0 = Q/A) and the fitted A is returned in `par`.
#' @param mass  tracer mass injected, in units consistent with obs_conc (see
#'   tsm_model.R's `simulate_tsm` documentation)
#' @param n_lhs number of global-search draws
#' @param n_cells passed through to simulate_tsm() (ignored by the exact
#'   solver; kept for compatibility)
#' @param bounds_log10 list(D=c(min,max), alpha=c(min,max), As_ratio=c(min,max))
#'   in log10 units; As_ratio = As/A. May also contain v=c(min,max) (log10
#'   m/s); if absent and fit_velocity = TRUE, v is searched within
#'   v0 / v_factor .. v0 * v_factor.
#' @param seed if not NULL, `set.seed(seed)` right before the LHS draw, so
#'   this call's result depends only on `seed` (see `job_seed()`), not on
#'   the RNG state left over from whatever ran before it. NULL (default)
#'   preserves the old behaviour of using the ambient RNG state as-is.
#' @param fit_velocity FALSE (default, original behaviour): A fixed.
#'   TRUE: v = Q/A fitted as a fourth parameter, Q fixed.
#' @param v_factor width of the default velocity search (multiplicative)
#' @param zero_pre_arrival if TRUE, set the pre-arrival stretch of the curve
#'   to 0 before fitting (see zero_pre_arrival()). Default FALSE (original
#'   behaviour). Meant for logger curves with baseline drift.
#' @param onset_frac threshold passed to zero_pre_arrival()
#' @param n_starts number of local (Nelder-Mead) searches, started from the
#'   n_starts best LHS points; each is restarted once from its own end
#'   point, and the best result is kept. With 4 free parameters a single
#'   Nelder-Mead run from the single best LHS point was found to stop on a
#'   ridge worse than the 3-parameter optimum (CB_20230907_single: RMSE
#'   0.438 with v free vs 0.423 with v fixed, impossible at a true optimum
#'   since v fixed is a special case). n_starts = 1 without restart is the
#'   original behaviour.
#' @return list(par = c(D, alpha, As, A, v), rmse, A, v, lhs = data.frame of
#'   the global scan (log10 units, plus `rmse`; includes a `v` column when
#'   fit_velocity = TRUE) for an identifiability plot, fit_velocity, bounds,
#'   n_pre_zeroed, t_onset, n_fit = number of points the cost was computed
#'   on, for AICc model comparison)
fit_hydraulics <- function(obs_time, obs_conc, L, Q, A, mass,
                           n_lhs = 250, n_cells = 40,
                           bounds_log10 = list(D = c(-4, 1),
                                               alpha = c(-6, -1),
                                               As_ratio = c(-3, 1)),
                           seed = NULL,
                           fit_velocity = FALSE, v_factor = 2,
                           zero_pre_arrival = FALSE, onset_frac = 0.05,
                           n_starts = 3) {
  
  zp <- list(value = obs_conc, n_zeroed = 0L, t_onset = NA_real_)
  if (zero_pre_arrival) zp <- zero_pre_arrival(obs_time, obs_conc, onset_frac)
  thin <- thin_series(obs_time, zp$value, n_max = 200)
  v0 <- Q / A
  
  pnames <- c("D", "alpha", "As_ratio")
  if (fit_velocity) {
    if (is.null(bounds_log10$v)) bounds_log10$v <- log10(v0) + c(-1, 1) * log10(v_factor)
    pnames <- c(pnames, "v")
  }
  bounds_log10 <- bounds_log10[pnames]
  
  # log10 parameter vector -> natural-unit parameters, clamped to bounds
  unpack <- function(log_par) {
    lp <- clamp_to_bounds(log_par, bounds_log10)
    v  <- if (fit_velocity) 10^lp[4] else v0
    Af <- Q / v
    list(D = 10^lp[1], alpha = 10^lp[2], As = 10^lp[3] * Af, A = Af, v = v)
  }
  
  cost_fun <- function(log_par) {
    p <- unpack(log_par)
    sim <- tryCatch(
      simulate_tsm(L = L, Q = Q, A = p$A, D = p$D, alpha = p$alpha, As = p$As,
                   mass = mass, times = thin$time, n_cells = n_cells)$C,
      error = function(e) rep(NA_real_, nrow(thin))
    )
    if (any(!is.finite(sim))) return(1e6)
    sqrt_rmse(sim, thin$value)
  }
  
  if (!is.null(seed)) set.seed(seed)
  lhs <- lhs_sample(bounds_log10, n_lhs)
  lhs$rmse <- apply(lhs, 1, function(r) cost_fun(as.numeric(r[pnames])))
  
  ctrl <- list(maxit = if (fit_velocity) 800 else 400, reltol = 1e-9)
  starts <- lhs[order(lhs$rmse), pnames, drop = FALSE]
  starts <- starts[seq_len(min(max(1L, n_starts), nrow(starts))), , drop = FALSE]
  fit <- NULL
  for (k in seq_len(nrow(starts))) {
    f <- optim(as.numeric(starts[k, ]), cost_fun, method = "Nelder-Mead", control = ctrl)
    if (n_starts > 1) {   # restart from the end point: escapes NM's premature stops
      f <- optim(f$par, cost_fun, method = "Nelder-Mead", control = ctrl)
    }
    if (is.null(fit) || f$value < fit$value) fit <- f
  }
  
  p <- unpack(fit$par)
  par <- c(D = p$D, alpha = p$alpha, As = p$As, A = p$A, v = p$v)
  list(par = par, rmse = fit$value, A = p$A, v = p$v, lhs = lhs,
       fit_velocity = fit_velocity, bounds_log10 = bounds_log10,
       n_pre_zeroed = zp$n_zeroed, t_onset = zp$t_onset, n_fit = nrow(thin))
}

#' Stage 2 -- fit lambda, lambda_s to a nutrient breakthrough curve, with
#' hydraulics (D, alpha, As) fixed from Stage 1
#'
#' @param obs_time,obs_conc observed grab time (s) and background-corrected
#'   nutrient concentration (e.g. ug/L)
#' @param hydraulics named vector c(D, alpha, As) from fit_hydraulics()$par
#' @param mass nutrient mass injected, in units consistent with obs_conc
#' @param bounds_log10 list(lambda=c(min,max), lambda_s=c(min,max)) in log10
#'   units (1/s)
#' @param seed if not NULL, `set.seed(seed)` right before the LHS draw --
#'   see the `seed` argument of `fit_hydraulics()` and `job_seed()` above
#'   for why this matters (this is the Stage-2 half of the same fix).
#' @return list(par = c(lambda, lambda_s), rmse, lhs)
fit_uptake <- function(obs_time, obs_conc, L, Q, A, hydraulics, mass,
                       n_lhs = 250, n_cells = 40,
                       bounds_log10 = list(lambda = c(-7, -1),
                                           lambda_s = c(-7, -1)),
                       seed = NULL) {
  
  D  <- hydraulics[["D"]]
  al <- hydraulics[["alpha"]]
  As <- hydraulics[["As"]]
  
  clamp <- function(x, b) pmin(pmax(x, b[1]), b[2])
  cost_fun <- function(log_par) {
    log_par <- c(clamp(log_par[1], bounds_log10$lambda),
                 clamp(log_par[2], bounds_log10$lambda_s))
    lam  <- 10^log_par[1]
    lams <- 10^log_par[2]
    sim <- tryCatch(
      simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = al, As = As,
                   lambda = lam, lambda_s = lams,
                   mass = mass, times = obs_time, n_cells = n_cells)$C,
      error = function(e) rep(NA_real_, length(obs_time))
    )
    if (any(!is.finite(sim))) return(1e6)
    sqrt_rmse(sim, obs_conc)
  }
  
  if (!is.null(seed)) set.seed(seed)
  lhs <- lhs_sample(bounds_log10, n_lhs)
  lhs$rmse <- apply(lhs, 1, function(r) cost_fun(as.numeric(r[c("lambda", "lambda_s")])))
  
  best <- lhs[which.min(lhs$rmse), c("lambda", "lambda_s")]
  fit <- optim(as.numeric(best), cost_fun, method = "Nelder-Mead",
               control = list(maxit = 300, reltol = 1e-8))
  
  lp  <- clamp_to_bounds(fit$par, bounds_log10[c("lambda", "lambda_s")])
  par <- c(lambda = 10^lp[1], lambda_s = 10^lp[2])
  list(par = par, rmse = fit$value, lhs = lhs)
}