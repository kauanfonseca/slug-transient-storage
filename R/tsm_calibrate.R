################################################################################
# tsm_calibrate.R
#
# Two-stage calibration of the transient-storage model (tsm_model.R) without
# needing OTIS-P, following the standard two-step logic of Runkel (2007):
#
#   Stage 1 (hydraulics): fit D, alpha, As against the CONSERVATIVE tracer
#            (NaCl) breakthrough curve, which is dense (logger, seconds
#            resolution) and well constrained. Q, A (=Q/v) are taken as
#            already known from the independent dilution-gauging discharge
#            estimate in `events.csv` -- these are of good, field-checked
#            quality (see the project handoff) and fixing them removes two
#            of the parameters that most threaten identifiability (Bonanno
#            et al. 2022).
#
#   Stage 2 (uptake): with D, alpha, As FIXED at their Stage-1 values, fit
#            lambda (main-channel) and lambda_s (storage-zone) first-order
#            uptake coefficients against the sparse nutrient grab series.
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
################################################################################

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

#' Thin a time series to at most `n_max` points for faster fitting, always
#' keeping the first and last point. Used only inside the optimisation loop;
#' the final reported fit is re-evaluated / plotted against the full series.
thin_series <- function(time, value, n_max = 200) {
  n <- length(time)
  if (n <= n_max) return(data.frame(time = time, value = value))
  idx <- unique(round(seq(1, n, length.out = n_max)))
  data.frame(time = time[idx], value = value[idx])
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

#' Stage 1 -- fit D, alpha, As to a conservative-tracer breakthrough curve
#'
#' @param obs_time,obs_conc observed time (s) and background-corrected
#'   conservative-tracer concentration (any consistent unit, e.g. mg/L)
#' @param L,Q,A reach length (m), discharge (m3/s), main-channel area (m2)
#' @param mass  tracer mass injected, in units consistent with obs_conc (see
#'   tsm_model.R's `simulate_tsm` documentation)
#' @param n_lhs number of global-search draws
#' @param n_cells spatial resolution
#' @param bounds_log10 list(D=c(min,max), alpha=c(min,max), As_ratio=c(min,max))
#'   in log10 units; As_ratio = As/A
#' @return list(par = c(D, alpha, As), rmse, lhs = data.frame of the global
#'   scan with an added `rmse` column, for an identifiability plot)
fit_hydraulics <- function(obs_time, obs_conc, L, Q, A, mass,
                            n_lhs = 250, n_cells = 40,
                            bounds_log10 = list(D = c(-4, 1),
                                                 alpha = c(-6, -1),
                                                 As_ratio = c(-3, 1))) {

  thin <- thin_series(obs_time, obs_conc, n_max = 200)

  clamp <- function(x, b) pmin(pmax(x, b[1]), b[2])
  cost_fun <- function(log_par) {
    log_par <- c(clamp(log_par[1], bounds_log10$D),
                 clamp(log_par[2], bounds_log10$alpha),
                 clamp(log_par[3], bounds_log10$As_ratio))
    D  <- 10^log_par[1]
    al <- 10^log_par[2]
    As <- 10^log_par[3] * A
    sim <- tryCatch(
      simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = al, As = As,
                   mass = mass, times = thin$time, n_cells = n_cells)$C,
      error = function(e) rep(NA_real_, nrow(thin))
    )
    if (any(!is.finite(sim))) return(1e6)
    sqrt_rmse(sim, thin$value)
  }

  lhs <- lhs_sample(bounds_log10, n_lhs)
  lhs$rmse <- apply(lhs, 1, function(r) cost_fun(as.numeric(r[c("D", "alpha", "As_ratio")])))

  best <- lhs[which.min(lhs$rmse), c("D", "alpha", "As_ratio")]
  fit <- optim(as.numeric(best), cost_fun, method = "Nelder-Mead",
               control = list(maxit = 300, reltol = 1e-8))

  par <- c(D = 10^fit$par[1], alpha = 10^fit$par[2], As = 10^fit$par[3] * A)
  list(par = par, rmse = fit$value, A = A, lhs = lhs)
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
#' @return list(par = c(lambda, lambda_s), rmse, lhs)
fit_uptake <- function(obs_time, obs_conc, L, Q, A, hydraulics, mass,
                        n_lhs = 250, n_cells = 40,
                        bounds_log10 = list(lambda = c(-7, -1),
                                             lambda_s = c(-7, -1))) {

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

  lhs <- lhs_sample(bounds_log10, n_lhs)
  lhs$rmse <- apply(lhs, 1, function(r) cost_fun(as.numeric(r[c("lambda", "lambda_s")])))

  best <- lhs[which.min(lhs$rmse), c("lambda", "lambda_s")]
  fit <- optim(as.numeric(best), cost_fun, method = "Nelder-Mead",
               control = list(maxit = 300, reltol = 1e-8))

  par <- c(lambda = 10^fit$par[1], lambda_s = 10^fit$par[2])
  list(par = par, rmse = fit$value, lhs = lhs)
}
