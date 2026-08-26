################################################################################
# plot_tsm_fit.R
#
# Diagnostic plots: (1) observed vs fitted breakthrough curve for a given
# event/solute, and (2) a lightweight identifiability diagnostic (parameter
# value vs RMSE across the global LHS scan; c.f. Bonanno et al. 2022's
# "global identifiability analysis" plots, Fig. 2) -- useful to eyeball
# whether a given stream's fit is trustworthy or sitting on a flat/
# non-identifiable ridge before trusting its numbers.
################################################################################

suppressPackageStartupMessages(library(ggplot2))

#' Plot observed vs simulated concentration for a fitted conservative or
#' nutrient breakthrough curve
#'
#' @param obs_time,obs_conc observed series
#' @param L,Q,A,D,alpha,As,lambda,lambda_s,mass as in simulate_tsm()
#' @param title plot title
#' @param unit y-axis unit label
plot_btc_fit <- function(obs_time, obs_conc, L, Q, A, D, alpha, As,
                          lambda = 0, lambda_s = 0, mass, title = "", unit = "mg/L",
                          n_cells = 40) {
  sim_times <- seq(0, max(obs_time), length.out = 400)
  sim <- simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = alpha, As = As,
                       lambda = lambda, lambda_s = lambda_s, mass = mass,
                       times = sim_times, n_cells = n_cells)
  ggplot() +
    geom_point(aes(x = obs_time, y = obs_conc), color = "grey30", size = 1.3, alpha = 0.7) +
    geom_line(data = sim, aes(x = time, y = C), color = "#1b6ca8", linewidth = 0.9) +
    labs(title = title, x = "time since release (s)", y = paste0("concentration (", unit, ")")) +
    theme_minimal(base_size = 12)
}

#' Bonanno-style "lite" identifiability plot: parameter value (natural units)
#' vs RMSE across the global LHS scan, for one parameter of a fit_hydraulics()
#' or fit_uptake() result. A clear, narrow minimum = identifiable; a flat
#' scatter across the whole range = poorly identified given this BTC.
#'
#' @param lhs data.frame returned as fit$lhs (has an `rmse` column plus one
#'   column per parameter, already in log10 units)
#' @param param column name to plot (e.g. "D", "alpha", "As_ratio", "lambda")
plot_identifiability <- function(lhs, param, natural_units = TRUE, xlab = param) {
  df <- lhs
  x <- if (natural_units) 10^df[[param]] else df[[param]]
  ggplot(df, aes(x = x, y = rmse)) +
    geom_point(alpha = 0.35, size = 1) +
    scale_x_log10() +
    labs(x = xlab, y = "RMSE (sqrt-concentration scale)",
         title = paste("Identifiability scan:", param)) +
    theme_minimal(base_size = 12)
}
