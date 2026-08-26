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
