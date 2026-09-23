################################################################################
# diag_tsm_shape.R
#
# Troubleshooting why the Stage-1 transient-storage fit (fit_hydraulics(),
# tsm_calibrate.R) fails to reproduce the "shark fin" of the conservative-
# tracer breakthrough curve. Built for ONE event (EVENT_ID below), plus a
# cross-event screen at the end (item 9).
#
# Items (numbering follows the discussion):
#   2. Numerical-dispersion floor of the upwind grid vs the fitted D
#   3. Grid convergence: (a) pure ADE, numeric vs exact analytic solution;
#      (b) the fitted TSM re-simulated at finer grids; (c) refit on a fine grid;
#      (d) analytic ADE fitted with v and D free (no grid error at all)
#   4. Temporal moments (mass recovery, centroid, variance, skewness),
#      observed vs model
#   5. Velocity lock: is v = Q/A double-counting storage delay? Refits with
#      v free, and with v + mass recovery free
#   6. Baseline drift: raw SpC vs the background used, pre-arrival ramp,
#      clipped tail values, share of the sqrt-cost spent on the ramp
#   7. Cost-function / drift sensitivity: linear RMSE, pre-arrival zeroed,
#      pre-arrival excluded, and a combined variant
#   8. Residuals vs time, and mean residual by curve segment, for every variant
#   9. Cross-event screen: which events share the same symptoms
#
# Standalone. Does not modify any pipeline file or output. Writes to
#   data_derived/diagnostics/<EVENT_ID>/
#
# Run time: the refits in items 3c/5/7 run on an N_FINE-cell grid with a
# banded-Jacobian solver; expect a few minutes in total. Lower N_FINE to 120
# if it is too slow.
#
# NOTE on the pipeline fit being diagnosed: the reported D = 8.71e-05 m2/s is
# BELOW the lower bound of the Stage-1 search (log10 D in [-4, 1], i.e.
# 1e-4). fit_hydraulics() clamps parameters inside the cost function but
# returns the UNCLAMPED Nelder-Mead vector, so the value that was actually
# simulated is D = 1e-4. In other words, D sat on its floor: the optimiser
# wanted less dispersion than the model could give it. Items 2-3 test
# whether that floor is the search bound or the grid.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(deSolve)
  library(here)
})

source(here::here("R", "tsm_model.R"))      # simulate_tsm_fv() (legacy solver)
source(here::here("R", "tsm_calibrate.R"))  # thin_series(), sqrt_rmse()
source(here::here("R", "tsm_partition.R"))  # trapz()

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
EVENT_ID   <- "CB_20230907_single"
FIT_PAR    <- c(D = 8.7137e-05, alpha = 0.00078832, As = 0.7328)  # app_tsm_review output
N_BASE     <- 40                      # grid used by the pipeline fit
N_GRID     <- c(40, 80, 160, 320)     # grid-convergence test (item 3b)
N_FINE     <- 160                     # grid for the diagnostic refits
ONSET_FRAC <- 0.05                    # pre-arrival window ends at the last point
# before the peak below ONSET_FRAC * peak
MAXIT      <- 400

# search bounds (log10) for the diagnostic refits. D's lower bound is opened
# to 1e-5 on purpose: once grid dispersion is reduced, the true optimum may
# sit below the pipeline's 1e-4 floor.
PAR_BOUNDS <- list(D = c(-5, 1), alpha = c(-6, -1), As_ratio = c(-3, 1),
                   v = c(-3, 1), mfrac = c(log10(0.5), log10(2)))

OUT_DIR <- here::here("data_derived", "diagnostics", EVENT_ID)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 5.5) {
  print(p)
  ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 150, bg = "white")
}
save_tbl <- function(x, name) {
  write_csv(x, file.path(OUT_DIR, name))
  print(x, n = Inf, width = Inf)
}
rmse_lin <- function(sim, obs) sqrt(mean((pmax(sim, 0) - obs)^2, na.rm = TRUE))
cost_fn  <- function(sim, obs, type) {
  switch(type, sqrt = sqrt_rmse(sim, obs), linear = rmse_lin(sim, obs))
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
events           <- read_csv(here("data_derived", "events.csv"), show_col_types = FALSE)
btc_conservative <- read_csv(here("data_derived", "btc_conservative.csv"), show_col_types = FALSE)

e <- events %>% filter(event_id == EVENT_ID)
stopifnot(nrow(e) == 1)

L      <- e$reach_length_m
Q      <- e$discharge_Ls / 1000        # m3/s
v_in   <- e$water_velocity_ms          # m/s, derived from the tracer curve
A_in   <- Q / v_in                     # m2
mass_g <- e$nacl_mass_g                # g  -> model output in mg/L

obs <- btc_conservative %>%
  filter(event_id == EVENT_ID, time_since_release_s >= 0) %>%
  arrange(time_since_release_s) %>%
  distinct(time_since_release_s, .keep_all = TRUE)

t_obs <- obs$time_since_release_s
c_raw <- obs$nacl_mgL
c_obs <- pmax(c_raw, 0)                # what the pipeline fits

i_pk <- which.max(c_obs); t_pk <- t_obs[i_pk]; c_pk <- c_obs[i_pk]
t_onset <- t_obs[max(which(c_obs[seq_len(i_pk)] < ONSET_FRAC * c_pk))]
pre_idx <- t_obs < t_onset
T_PLOT  <- seq(0, max(t_obs), length.out = 600)

segment_of <- function(t, c) {
  case_when(t < t_onset              ~ "1_pre_arrival",
            t < t_pk                 ~ "2_rising",
            c >= ONSET_FRAC * c_pk   ~ "3_falling",
            TRUE                     ~ "4_tail")
}
seg_obs <- segment_of(t_obs, c_obs)

cat(sprintf(paste0("\n%s: L = %.1f m | Q = %.1f L/s | v_in = %.4f m/s | A_in = %.3f m2 | ",
                   "NaCl = %.0f g | n_obs = %d\n"),
            EVENT_ID, L, Q * 1000, v_in, A_in, mass_g, length(t_obs)))
cat(sprintf("As/A check: As/A_in = %.3f (app reported 0.156)\n", FIT_PAR[["As"]] / A_in))
cat(sprintf("peak %.2f mg/L at %.0f s | onset (%.0f%% of peak) at %.0f s | %d pre-arrival points\n",
            c_pk, t_pk, 100 * ONSET_FRAC, t_onset, sum(pre_idx)))

# ---------------------------------------------------------------------------
# Fast TSM solver used for the diagnostics
#
# Same physics and discretisation as simulate_tsm() (first-order upwind
# advection, central dispersion, slug in cell 1, sample in the last cell) but
#   * no 120-cell cap (tsm_grid() clamps n_cells to [10, 120]);
#   * channel and storage states INTERLEAVED (C1, Cs1, C2, Cs2, ...), which
#     makes the Jacobian banded (2 up / 2 down), so lsoda can use a banded
#     solver and fine grids stay affordable;
#   * alpha = 0 allowed (pure ADE) without the A/As division blowing up.
# Checked against simulate_tsm() at n = 40 right below.
# ---------------------------------------------------------------------------
tsm_derivs_il <- function(t, y, p) {
  C  <- y[c(TRUE, FALSE)]
  Cs <- y[c(FALSE, TRUE)]
  n  <- p$n
  C_up <- c(0, C[-n])
  C_dn <- c(C[-1], C[n])
  dC  <- -p$v * (C - C_up) / p$dx + p$D * (C_dn - 2 * C + C_up) / p$dx^2 +
    p$alpha * (Cs - C) - p$lambda * C
  dCs <- p$gamma * (C - Cs) - p$lambda_s * Cs
  out <- numeric(2 * n)
  out[c(TRUE, FALSE)] <- dC
  out[c(FALSE, TRUE)] <- dCs
  list(out)
}

sim_tsm2 <- function(L, Q, A, D, alpha, As, mass, times, n_cells,
                     lambda = 0, lambda_s = 0) {
  n  <- round(n_cells)
  dx <- L / n
  y0 <- numeric(2 * n)
  y0[1] <- mass / (A * dx)
  p <- list(n = n, dx = dx, v = Q / A, D = D, alpha = alpha,
            gamma = if (alpha > 0 && As > 0) alpha * A / As else 0,
            lambda = lambda, lambda_s = lambda_s)
  pre <- times[1] != 0
  tt  <- if (pre) c(0, times) else times
  out <- suppressWarnings(
    ode(y = y0, times = tt, func = tsm_derivs_il, parms = p, method = "lsoda",
        jactype = "bandint", bandup = 2, banddown = 2)
  )
  out <- as.data.frame(out)
  if (pre) out <- out[-1, , drop = FALSE]
  data.frame(time = out$time, C = out[[2 * n]])   # column 2n = C of the last cell
}

# exact ADE solution (no storage, no grid), instantaneous slug, resident
# concentration at x = L
ade_analytic <- function(t, L, A, v, D, mass) {
  out <- numeric(length(t))
  pos <- t > 0
  out[pos] <- mass / (A * sqrt(4 * pi * D * t[pos])) *
    exp(-(L - v * t[pos])^2 / (4 * D * t[pos]))
  out
}

moments <- function(t, c) {
  m0 <- trapz(t, c)
  if (!is.finite(m0) || m0 <= 0) {
    return(tibble(M0 = NA_real_, t_centroid = NA_real_, variance = NA_real_,
                  skewness = NA_real_, t_peak = NA_real_, c_peak = NA_real_))
  }
  t1  <- trapz(t, t * c) / m0
  mu2 <- trapz(t, (t - t1)^2 * c) / m0
  mu3 <- trapz(t, (t - t1)^3 * c) / m0
  tibble(M0 = m0, t_centroid = t1, variance = mu2, skewness = mu3 / mu2^1.5,
         t_peak = t[which.max(c)], c_peak = max(c))
}

# sanity check: fast solver == pipeline solver at n = 40
chk_pipe <- simulate_tsm_fv(L = L, Q = Q, A = A_in, D = FIT_PAR[["D"]], alpha = FIT_PAR[["alpha"]],
                            As = FIT_PAR[["As"]], mass = mass_g, times = T_PLOT, n_cells = N_BASE)$C
chk_fast <- sim_tsm2(L, Q, A_in, FIT_PAR[["D"]], FIT_PAR[["alpha"]], FIT_PAR[["As"]],
                     mass_g, T_PLOT, N_BASE)$C
rel_diff <- max(abs(chk_pipe - chk_fast)) / max(chk_pipe)
cat(sprintf("\nsolver check (n = %d): max |pipeline - fast| / peak = %.2e\n", N_BASE, rel_diff))
if (rel_diff > 1e-3) warning("fast solver disagrees with simulate_tsm(); treat item 3+ with care")

# ---------------------------------------------------------------------------
# Shared evaluation of any fitted curve against the observed series
# ---------------------------------------------------------------------------
evaluate_curve <- function(label, sim_fun, par_row) {
  s_obs  <- sim_fun(t_obs)
  s_plot <- sim_fun(T_PLOT)
  c_pk_m <- max(s_plot)
  t_pk_m <- T_PLOT[which.max(s_plot)]
  metrics <- tibble(
    label          = label,
    rmse_sqrt_all  = sqrt_rmse(s_obs, c_obs),
    rmse_lin_all   = rmse_lin(s_obs, c_obs),
    rmse_sqrt_post = sqrt_rmse(s_obs[!pre_idx], c_obs[!pre_idx]),
    rmse_lin_post  = rmse_lin(s_obs[!pre_idx], c_obs[!pre_idx]),
    c_peak_model   = c_pk_m,
    peak_ratio     = c_pk_m / c_pk,
    peak_shift_s   = t_pk_m - t_pk
  ) %>% bind_cols(par_row)
  list(
    metrics = metrics,
    curve   = tibble(label = label, time = T_PLOT, C = s_plot),
    resid   = tibble(label = label, time = t_obs, obs = c_obs, sim = s_obs,
                     resid_lin = c_obs - s_obs,
                     resid_sqrt = sqrt(c_obs) - sqrt(pmax(s_obs, 0)),
                     segment = seg_obs)
  )
}

# Refit on the fast solver. Parameters in log10; `free` picks which move.
# Bounds are enforced by clamping, and the CLAMPED values are what is
# reported (unlike fit_hydraulics(), see header note).
START <- list(D = log10(FIT_PAR[["D"]]), alpha = log10(FIT_PAR[["alpha"]]),
              As_ratio = log10(FIT_PAR[["As"]] / A_in), v = log10(v_in), mfrac = 0)

fit_diag <- function(label, t_fit, c_fit, free, cost = "sqrt", n_cells = N_FINE,
                     starts = list(START), maxit = MAXIT) {
  th <- thin_series(t_fit, c_fit, n_max = 200)
  unpack <- function(lp, base) {
    p <- base
    p[free] <- as.list(lp)
    for (nm in free) p[[nm]] <- min(max(p[[nm]], PAR_BOUNDS[[nm]][1]), PAR_BOUNDS[[nm]][2])
    p
  }
  sim_p <- function(p, times) {
    v <- 10^p$v
    A <- Q / v
    sim_tsm2(L, Q, A, D = 10^p$D, alpha = 10^p$alpha, As = 10^p$As_ratio * A,
             mass = mass_g * 10^p$mfrac, times = times, n_cells = n_cells)$C
  }
  best <- NULL
  for (s in starts) {
    obj <- function(lp) {
      sim <- tryCatch(sim_p(unpack(lp, s), th$time), error = function(err) NA_real_)
      if (any(!is.finite(sim))) return(1e6)
      cost_fn(sim, th$value, cost)
    }
    f <- optim(unlist(s[free]), obj, method = "Nelder-Mead",
               control = list(maxit = maxit, reltol = 1e-8))
    if (is.null(best) || f$value < best$value) best <- c(f, list(start = s))
  }
  p <- unpack(best$par, best$start)
  at_bound <- vapply(free, function(nm) {
    b <- PAR_BOUNDS[[nm]]
    (p[[nm]] - b[1]) < 0.02 * diff(b) || (b[2] - p[[nm]]) < 0.02 * diff(b)
  }, logical(1))
  v <- 10^p$v
  A <- Q / v
  par_row <- tibble(
    D = 10^p$D, alpha = 10^p$alpha, As = 10^p$As_ratio * A, As_ratio = 10^p$As_ratio,
    v = v, A = A, mass_recovery = 10^p$mfrac, n_cells = n_cells, cost = cost,
    free = paste(free, collapse = "+"),
    at_bound = if (any(at_bound)) paste(free[at_bound], collapse = "+") else "",
    converged = best$convergence == 0
  )
  message(sprintf("  %-28s done (%s cost = %.4f)", label, cost, best$value))
  evaluate_curve(label, function(tt) sim_p(p, tt), par_row)
}

pipeline_row <- tibble(D = FIT_PAR[["D"]], alpha = FIT_PAR[["alpha"]], As = FIT_PAR[["As"]],
                       As_ratio = FIT_PAR[["As"]] / A_in, v = v_in, A = A_in,
                       mass_recovery = 1, n_cells = N_BASE, cost = "sqrt",
                       free = "(pipeline fit)", at_bound = "D (see header)", converged = NA)
V0 <- evaluate_curve(
  "V0_pipeline_n40",
  function(tt) simulate_tsm_fv(L = L, Q = Q, A = A_in, D = FIT_PAR[["D"]], alpha = FIT_PAR[["alpha"]],
                               As = FIT_PAR[["As"]], mass = mass_g, times = tt, n_cells = N_BASE)$C,
  pipeline_row
)

obs_df <- tibble(time = t_obs, conc = c_obs, pre = pre_idx)
base_obs_layer <- list(
  geom_point(data = obs_df, aes(time, conc), colour = "grey45", size = 1.1, alpha = 0.8),
  geom_vline(xintercept = t_onset, linetype = "dotted", colour = "grey40"),
  labs(x = "time since release (s)", y = "NaCl (mg/L)", colour = NULL),
  theme_minimal(base_size = 12),
  theme(legend.position = "bottom")
)

################################################################################
# 2. Numerical-dispersion floor
################################################################################
cat("\n=== 2. Numerical-dispersion floor (first-order upwind: D_num ~ v*dx/2) ===\n")
disp_tbl <- tibble(n_cells = sort(unique(c(N_GRID, 120)))) %>%
  mutate(dx_m           = L / n_cells,
         D_num_m2s      = v_in * dx_m / 2,
         D_fit_m2s      = FIT_PAR[["D"]],
         D_num_over_fit = D_num_m2s / D_fit_m2s,
         cell_Peclet    = v_in * dx_m / D_fit_m2s,
         D_eff_m2s      = D_fit_m2s + D_num_m2s)
save_tbl(disp_tbl, "02_numerical_dispersion.csv")
cat("D_num_over_fit >> 1 means the grid, not D, controls the spread of the curve.\n",
    "tsm_grid() caps n_cells at 120 in the pipeline, so the 120 row is the best the pipeline can do.\n")

################################################################################
# 3. Grid convergence
################################################################################
cat("\n=== 3a. Pure ADE (alpha = 0): numeric grid vs exact solution ===\n")
t_ade <- seq(0, 2.5 * L / v_in, length.out = 2000)
ade_exact <- ade_analytic(t_ade, L, A_in, v_in, FIT_PAR[["D"]], mass_g)
ade_curves <- map_dfr(N_GRID, function(n) {
  tibble(label = sprintf("numeric n=%d", n), time = t_ade,
         C = sim_tsm2(L, Q, A_in, FIT_PAR[["D"]], 0, 1, mass_g, t_ade, n)$C)
}) %>% bind_rows(tibble(label = "exact", time = t_ade, C = ade_exact))

ade_tbl <- ade_curves %>%
  group_by(label) %>%
  group_modify(~ moments(.x$time, .x$C)) %>%
  ungroup() %>%
  mutate(D_eff_from_variance = variance * v_in^3 / (2 * L))  # approx, for D << vL
save_tbl(ade_tbl, "03a_ade_numeric_vs_exact.csv")
cat("D_eff_from_variance for the numeric rows should track D + v*dx/2 (item 2);",
    "for 'exact' it should return D.\n")

save_plot(
  ggplot(ade_curves, aes(time, C, colour = label)) + geom_line(linewidth = 0.8) +
    labs(title = sprintf("3a. Pure ADE at D = %.2g m2/s: grid vs exact", FIT_PAR[["D"]]),
         x = "time since release (s)", y = "NaCl (mg/L)", colour = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom"),
  "03a_ade_numeric_vs_exact.png"
)

cat("\n=== 3b. Fitted TSM parameters re-simulated on finer grids (no refit) ===\n")
grid_curves <- map_dfr(N_GRID, function(n) {
  tibble(label = sprintf("n=%d", n), time = T_PLOT,
         C = sim_tsm2(L, Q, A_in, FIT_PAR[["D"]], FIT_PAR[["alpha"]], FIT_PAR[["As"]],
                      mass_g, T_PLOT, n)$C)
})
save_plot(
  ggplot() + base_obs_layer +
    geom_line(data = grid_curves, aes(time, C, colour = label), linewidth = 0.8) +
    labs(title = "3b. Same D, alpha, As; only the grid changes"),
  "03b_grid_convergence_fitted_params.png"
)

cat("\n=== 3c/5/7. Refits (this is the slow part) ===\n")

# item 5 needs a second starting v: the velocity that would make the model
# centroid L(1+beta)/v land on the observed centroid
mom_obs <- moments(t_obs, c_obs)
beta0   <- FIT_PAR[["As"]] / A_in
v_centroid_consistent <- L * (1 + beta0) / mom_obs$t_centroid
START_V2 <- modifyList(START, list(v = log10(v_centroid_consistent)))

c_prezero <- replace(c_obs, pre_idx, 0)

fits <- list(
  V0 = V0,
  # 3c: same free parameters and cost as the pipeline, finer grid
  V1 = fit_diag("V1_fine_sqrt", t_obs, c_obs, c("D", "alpha", "As_ratio")),
  # 5: velocity free
  V2 = fit_diag("V2_fine_sqrt_vfree", t_obs, c_obs, c("D", "alpha", "As_ratio", "v"),
                starts = list(START, START_V2)),
  # 4/5: velocity and mass recovery free
  V3 = fit_diag("V3_fine_sqrt_vfree_mfree", t_obs, c_obs,
                c("D", "alpha", "As_ratio", "v", "mfrac"), starts = list(START, START_V2)),
  # 7: cost function and pre-arrival treatment (v fixed, as in the pipeline)
  V4 = fit_diag("V4_fine_linear", t_obs, c_obs, c("D", "alpha", "As_ratio"), cost = "linear"),
  V5 = fit_diag("V5_fine_sqrt_prezero", t_obs, c_prezero, c("D", "alpha", "As_ratio")),
  V6 = fit_diag("V6_fine_sqrt_preexcl", t_obs[!pre_idx], c_obs[!pre_idx],
                c("D", "alpha", "As_ratio")),
  # combined: everything suspected fixed at once
  V7 = fit_diag("V7_fine_sqrt_prezero_vfree", t_obs, c_prezero,
                c("D", "alpha", "As_ratio", "v"), starts = list(START, START_V2))
)

# 3d: exact ADE (no storage, no grid), v and D free, mass recovery free
fit_ade <- local({
  lp0 <- c(v = log10(L / t_pk),
           D = log10(max(mom_obs$variance * (L / t_pk)^3 / (2 * L), 1e-5)),
           mfrac = 0)
  sim_p <- function(lp, tt) {
    v <- 10^lp[["v"]]
    ade_analytic(tt, L, Q / v, v, 10^lp[["D"]], mass_g * 10^lp[["mfrac"]])
  }
  f <- optim(lp0, function(lp) {
    lp[["mfrac"]] <- min(max(lp[["mfrac"]], PAR_BOUNDS$mfrac[1]), PAR_BOUNDS$mfrac[2])
    s <- sim_p(lp, t_obs)
    if (any(!is.finite(s))) return(1e6)
    sqrt_rmse(s, c_obs)
  }, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-10))
  lp <- f$par
  lp[["mfrac"]] <- min(max(lp[["mfrac"]], PAR_BOUNDS$mfrac[1]), PAR_BOUNDS$mfrac[2])
  v <- 10^lp[["v"]]
  evaluate_curve("V8_exact_ADE_vfree_mfree", function(tt) sim_p(lp, tt),
                 tibble(D = 10^lp[["D"]], alpha = 0, As = 0, As_ratio = 0, v = v, A = Q / v,
                        mass_recovery = 10^lp[["mfrac"]], n_cells = NA_real_, cost = "sqrt",
                        free = "v+D+mfrac (analytic)", at_bound = "",
                        converged = f$convergence == 0))
})
fits$V8 <- fit_ade

variant_tbl <- map_dfr(fits, "metrics")
all_curves  <- map_dfr(fits, "curve")
all_resid   <- map_dfr(fits, "resid")

cat("\n=== Variant comparison (all metrics on the SAME observed series) ===\n")
save_tbl(variant_tbl %>% arrange(rmse_lin_post), "03-07_variant_comparison.csv")

save_plot(
  ggplot() + base_obs_layer +
    geom_line(data = all_curves, aes(time, C, colour = label), linewidth = 0.7) +
    facet_wrap(~ label, ncol = 3) + theme(legend.position = "none") +
    labs(title = sprintf("%s: all variants (dotted line = onset)", EVENT_ID)),
  "03-07_variant_fits.png", w = 13, h = 10
)

################################################################################
# 4. Temporal moments
################################################################################
cat("\n=== 4. Temporal moments ===\n")
t_long <- seq(0, max(choose_horizon(L, Q, A_in, FIT_PAR[["alpha"]], FIT_PAR[["As"]])$T_max,
                     2 * max(t_obs)), length.out = 4000)
mom_tbl <- bind_rows(
  moments(t_obs, c_obs) %>% mutate(series = "observed (clamped >= 0)", window = "observed"),
  moments(t_obs, c_raw) %>% mutate(series = "observed (raw)",          window = "observed"),
  moments(t_obs, c_obs * !pre_idx) %>% mutate(series = "observed, pre-arrival zeroed", window = "observed"),
  moments(T_PLOT, fits$V0$curve$C) %>% mutate(series = "pipeline fit, n=40", window = "observed"),
  moments(t_long, sim_tsm2(L, Q, A_in, FIT_PAR[["D"]], FIT_PAR[["alpha"]], FIT_PAR[["As"]],
                           mass_g, t_long, N_BASE)$C) %>%
    mutate(series = "pipeline fit, n=40", window = "full horizon"),
  moments(t_long, sim_tsm2(L, Q, A_in, FIT_PAR[["D"]], FIT_PAR[["alpha"]], FIT_PAR[["As"]],
                           mass_g, t_long, N_FINE)$C) %>%
    mutate(series = sprintf("pipeline params, n=%d", N_FINE), window = "full horizon")
) %>%
  mutate(mass_recovery = Q * M0 / mass_g) %>%   # (m3/s)(mg.s/L)/g = dimensionless
  relocate(series, window)
save_tbl(mom_tbl, "04_temporal_moments.csv")
cat("mass_recovery: Q * M0 / mass. The model's full-horizon row is 1 by construction;",
    "an observed value far from 1 means the fixed Q and the logger curve disagree,",
    "so the model cannot match peak height and area at the same time.\n",
    "skewness: the shark fin should show as clearly positive; compare the model rows.\n")

################################################################################
# 5. Velocity lock
################################################################################
cat("\n=== 5. Velocity: is storage delay counted twice? ===\n")
t_adv <- L / v_in
vel_tbl <- tibble(
  L_m = L, v_in_ms = v_in, t_adv_L_over_v_s = t_adv,
  t_peak_obs_s = t_pk, t_centroid_obs_s = mom_obs$t_centroid,
  ratio_tadv_to_tpeak     = t_adv / t_pk,
  ratio_tadv_to_tcentroid = t_adv / mom_obs$t_centroid,
  beta_fit = beta0,
  model_centroid_expected_s = t_adv * (1 + beta0),   # TSM first moment, L(1+beta)/v
  model_centroid_bias_pct   = 100 * (t_adv * (1 + beta0) / mom_obs$t_centroid - 1),
  v_centroid_consistent_ms  = v_centroid_consistent,
  v_refit_V2_ms = variant_tbl$v[variant_tbl$label == "V2_fine_sqrt_vfree"],
  v_refit_V7_ms = variant_tbl$v[variant_tbl$label == "V7_fine_sqrt_prezero_vfree"]
)
save_tbl(vel_tbl, "05_velocity_lock.csv")
cat("A ratio_tadv_to_tpeak or ratio_tadv_to_tcentroid of ~1 identifies which one v came from.\n",
    "If v came from the centroid, the model centroid is late by ~beta (model_centroid_bias_pct),\n",
    "and the optimiser can only compensate by shrinking As/A, which removes the tail.\n",
    "Compare V1 vs V2/V3 in the variant table: a large drop in RMSE and a larger As_ratio\n",
    "once v is free points to the velocity lock.\n")

################################################################################
# 6. Baseline drift
################################################################################
cat("\n=== 6. Baseline drift and clipping ===\n")
pre_df   <- obs_df %>% filter(pre)
drift_lm <- if (nrow(pre_df) >= 5) lm(conc ~ time, data = pre_df) else NULL
tail_idx <- seg_obs == "4_tail"

sq_res_v0 <- fits$V0$resid$resid_sqrt^2
drift_tbl <- tibble(
  t_onset_s             = t_onset,
  n_pre_arrival         = sum(pre_idx),
  pre_mean_mgL          = mean(c_obs[pre_idx]),
  pre_slope_mgL_per_h   = if (!is.null(drift_lm)) coef(drift_lm)[[2]] * 3600 else NA_real_,
  pre_level_at_onset    = if (!is.null(drift_lm)) predict(drift_lm, tibble(time = t_onset)) else NA_real_,
  pre_level_pct_of_peak = 100 * pre_level_at_onset / c_pk,
  n_negative_raw        = sum(c_raw < 0, na.rm = TRUE),
  n_negative_raw_tail   = sum(c_raw[tail_idx] < 0, na.rm = TRUE),
  tail_mean_mgL         = mean(c_obs[tail_idx]),
  last_value_mgL        = tail(c_raw, 1),
  # share of the pipeline fit's squared sqrt-residuals coming from the ramp
  pct_sqrt_cost_from_pre_arrival = 100 * sum(sq_res_v0[pre_idx]) / sum(sq_res_v0),
  pct_points_pre_arrival         = 100 * mean(pre_idx)
)
save_tbl(drift_tbl, "06_baseline_drift.csv")
cat("If pct_sqrt_cost_from_pre_arrival is much larger than pct_points_pre_arrival,",
    "the ramp is steering the fit.\n")

# raw SpC, if the conservative table carries it
spc_cols <- intersect(c("spc_uscm", "spc_corr_uscm", "background_spc", "background_uscm"),
                      names(obs))
if (length(spc_cols) > 0) {
  save_plot(
    obs %>% select(time_since_release_s, all_of(spc_cols)) %>%
      pivot_longer(-time_since_release_s) %>%
      ggplot(aes(time_since_release_s, value, colour = name)) + geom_line() +
      geom_vline(xintercept = t_onset, linetype = "dotted") +
      labs(title = "6a. SpC columns in btc_conservative", x = "time since release (s)",
           y = "uS/cm", colour = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom"),
    "06a_spc_in_btc_conservative.png"
  )
}

# the whole logger day with the background actually used, and the slug window
hobo_path <- here("data", "hobo_with_meta.rds")
if (file.exists(hobo_path)) {
  hday <- readRDS(hobo_path) %>%
    filter(stream == e$stream, as.Date(date) == as.Date(e$date))
  if (nrow(hday) > 0) {
    p6 <- ggplot(hday, aes(datetime, spc_uscm)) +
      geom_line(linewidth = 0.4) +
      geom_line(aes(y = background_spc), colour = "red", linetype = "dashed") +
      facet_wrap(~ source_file, ncol = 1, scales = "free") +
      labs(title = sprintf("6b. Raw logger SpC, %s %s (red = background used)",
                           e$stream, as.Date(e$date)),
           x = NULL, y = "SpC (uS/cm)") +
      theme_minimal(base_size = 11)
    log_path <- here("data", "slug_trimmed_log.csv")
    if (file.exists(log_path)) {
      win <- read_csv(log_path, show_col_types = FALSE) %>%
        filter(stream == e$stream, as.Date(date) == as.Date(e$date))
      if (nrow(win) > 0) {
        p6 <- p6 + geom_rect(data = win, inherit.aes = FALSE,
                             aes(xmin = start_datetime, xmax = end_datetime,
                                 ymin = -Inf, ymax = Inf),
                             fill = "steelblue", alpha = 0.12)
      }
    }
    save_plot(p6, "06b_raw_logger_day.png", h = 3 + 2.5 * n_distinct(hday$source_file))
  } else {
    cat("(hobo_with_meta.rds has no rows for this stream/date)\n")
  }
}

################################################################################
# 8. Residuals
################################################################################
cat("\n=== 8. Residuals by curve segment ===\n")
seg_tbl <- all_resid %>%
  group_by(label, segment) %>%
  summarise(n = n(), mean_resid_lin = mean(resid_lin), mean_resid_sqrt = mean(resid_sqrt),
            .groups = "drop") %>%
  pivot_wider(id_cols = label, names_from = segment,
              values_from = c(mean_resid_lin), names_prefix = "mean_obs_minus_model_")
save_tbl(seg_tbl, "08_residuals_by_segment.csv")
cat("Positive = model too low. The shark-fin signature for a too-dispersive model is:\n",
    "negative early on the rising limb, positive late on the rising limb and at the peak, ",
    "negative on the early falling limb.\n")

save_plot(
  ggplot(all_resid, aes(time, resid_lin, colour = segment)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_point(size = 0.8) +
    facet_wrap(~ label, ncol = 3) +
    labs(title = "8. Observed - model (mg/L)", x = "time since release (s)",
         y = "residual (mg/L)", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom"),
  "08_residuals_linear.png", w = 13, h = 10
)
save_plot(
  ggplot(all_resid, aes(time, resid_sqrt, colour = segment)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_point(size = 0.8) +
    facet_wrap(~ label, ncol = 3) +
    labs(title = "8. sqrt(obs) - sqrt(model): what the pipeline cost actually sees",
         x = "time since release (s)", y = "residual (sqrt mg/L)", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom"),
  "08_residuals_sqrt.png", w = 13, h = 10
)

################################################################################
# 9. Cross-event screen (no refitting)
################################################################################
cat("\n=== 9. Cross-event screen ===\n")

# fitted hydraulics, from a run of the batch pipeline if it is in memory,
# else from the app's accepted overrides
hyd_tbl <- NULL
if (exists("summary_tbl") && !is.null(attr(summary_tbl, "hydraulic_fits"))) {
  hyd_tbl <- map_dfr(attr(summary_tbl, "hydraulic_fits"), function(x) {
    tibble(event_id = x$event_id, D_fit = x$D, alpha_fit = x$alpha, As_fit = x$As,
           rmse_fit = x$rmse, fit_source = x$conservative_source)
  })
} else if (file.exists(here("data_derived", "tsm_manual_overrides.csv"))) {
  hyd_tbl <- read_csv(here("data_derived", "tsm_manual_overrides.csv"),
                      col_types = cols(.default = "c")) %>%
    filter(stage == "hydraulics") %>%
    group_by(event_id) %>% slice_tail(n = 1) %>% ungroup() %>%
    transmute(event_id, D_fit = as.numeric(D_m2s), alpha_fit = as.numeric(alpha_1s),
              As_fit = as.numeric(As_m2), rmse_fit = as.numeric(rmse),
              fit_source = conservative_source)
}

obs_mom_all <- btc_conservative %>%
  filter(time_since_release_s >= 0) %>%
  arrange(event_id, time_since_release_s) %>%
  distinct(event_id, time_since_release_s, .keep_all = TRUE) %>%
  group_by(event_id) %>%
  filter(n() >= 8) %>%
  group_modify(~ moments(.x$time_since_release_s, pmax(.x$nacl_mgL, 0))) %>%
  ungroup()

cross <- events %>%
  select(event_id, stream, reach_length_m, discharge_Ls, water_velocity_ms, nacl_mass_g,
         any_of(c("discharge_source", "flag_mixing_incomplete", "flag_no_baseline",
                  "flag_weak_signal"))) %>%
  inner_join(obs_mom_all, by = "event_id") %>%
  mutate(
    Q_m3s             = discharge_Ls / 1000,
    t_adv_s           = reach_length_m / water_velocity_ms,
    ratio_tadv_tpeak  = t_adv_s / t_peak,
    ratio_tadv_tcent  = t_adv_s / t_centroid,
    mass_recovery     = Q_m3s * M0 / nacl_mass_g,
    D_num_n40         = water_velocity_ms * (reach_length_m / 40) / 2,
    D_num_n120        = water_velocity_ms * (reach_length_m / 120) / 2
  )
if (!is.null(hyd_tbl)) {
  cross <- cross %>%
    left_join(hyd_tbl, by = "event_id") %>%
    mutate(D_at_floor         = D_fit <= 1.05e-4,
           Dnum40_over_Dfit   = D_num_n40 / D_fit,
           beta_fit           = As_fit / (Q_m3s / water_velocity_ms),
           centroid_bias_pct  = 100 * (t_adv_s * (1 + beta_fit) / t_centroid - 1))
} else {
  cat("(no fitted hydraulics found: run run_tsm_uptake_all() into `summary_tbl` first,",
      "or accept fits in the app, to get the D/RMSE columns)\n")
}
save_tbl(cross %>% arrange(desc(D_num_n40)), "09_cross_event_screen.csv")

if (!is.null(hyd_tbl) && any(!is.na(cross$D_fit))) {
  save_plot(
    ggplot(cross, aes(D_num_n40, D_fit, colour = rmse_fit)) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = 1e-4, linetype = "dotted", colour = "red") +
      geom_point(size = 2.5) +
      geom_text(aes(label = event_id), size = 2.6, vjust = -0.9, show.legend = FALSE) +
      scale_x_log10() + scale_y_log10() +
      labs(title = "9. Fitted D vs grid dispersion at n = 40",
           subtitle = "above the dashed line: fitted D exceeds the grid floor | red dotted: search floor",
           x = "D_num = v*dx/2 at n=40 (m2/s)", y = "fitted D (m2/s)", colour = "RMSE") +
      theme_minimal(base_size = 12),
    "09_Dfit_vs_Dnum.png", h = 6
  )
}

################################################################################
# Summary
################################################################################
cat("\n==================== SUMMARY:", EVENT_ID, "====================\n")
v <- variant_tbl
show <- function(lbl) {
  r <- v[v$label == lbl, ]
  cat(sprintf("  %-30s RMSE lin(post) %.3f | peak %.0f%% of obs, shift %+.0f s | D %.2g, As/A %.3f, v %.4f, mrec %.2f %s\n",
              lbl, r$rmse_lin_post, 100 * r$peak_ratio, r$peak_shift_s, r$D, r$As_ratio,
              r$v, r$mass_recovery, if (nzchar(r$at_bound)) paste0("[bound: ", r$at_bound, "]") else ""))
}
walk(v$label, show)
cat("\nReading guide:\n",
    " V0 -> V1 : effect of the grid alone (item 3)\n",
    " V1 -> V2 : effect of freeing v (item 5)\n",
    " V2 -> V3 : effect of also freeing mass recovery (item 4)\n",
    " V1 -> V4/V5/V6 : effect of the cost function and pre-arrival ramp (item 7)\n",
    " V7 : everything suspected fixed at once\n",
    " V8 : exact ADE, no storage, no grid. If V8 beats V0, the pipeline's TSM is being\n",
    "      held back by its setup, not by its physics.\n",
    "Outputs written to: ", OUT_DIR, "\n")