################################################################################
# tsm_model.R
#
# One-dimensional transient-storage transport model (Bencala & Walters 1983;
# Runkel 1998, OTIS) for a SLUG (instantaneous) addition, solved EXACTLY in R
# -- no spatial grid, no ODE solver, no OTIS/OTIS-P executables.
#
# Governing equations (main channel + storage zone, constant Q, A along the
# reach -- appropriate for a short slug-injection reach with no significant
# lateral inflow):
#
#   dC/dt  = -v dC/dx + D d2C/dx2 + alpha*(Cs - C) - lambda*C
#   dCs/dt = alpha*(A/As)*(C - Cs)                 - lambda_s*Cs
#
# C  : main-channel solute concentration (mass/volume), background-corrected
# Cs : storage-zone concentration
# v  : mean velocity = Q/A
# D  : longitudinal dispersion coefficient
# alpha   : main channel <-> storage zone exchange coefficient
# As      : storage zone cross-sectional area
# lambda, lambda_s : first-order removal (uptake) coefficients, channel and
#                     storage zone. Zero for a conservative tracer (NaCl).
#
# ---------------------------------------------------------------------------
# WHY THIS FILE CHANGED (2026-09)
# ---------------------------------------------------------------------------
# The previous solver discretised the reach into n_cells finite volumes with
# first-order upwind advection. That scheme adds an artificial ("numerical")
# dispersion of about v*dx/2 on top of D. For CB_20230907_single (L = 162 m,
# v = 0.0817 m/s, 40 cells) that is ~0.165 m2/s -- about 1900x the fitted D,
# and on its own already wider than the observed breakthrough curve. The
# optimiser responded by driving D to its lower bound and bending alpha/As,
# so the fitted D, alpha, As described the grid's smearing rather than the
# stream (see R/diag_tsm_shape.R, items 2-3: the 40-cell model reproduced
# the correct mass and centroid but a peak 40x too low for a pure-ADE test).
# Refining the grid only helps slowly (still ~0.02 m2/s at 320 cells).
#
# The old solver is kept below as simulate_tsm_fv() for comparison only.
#
# ---------------------------------------------------------------------------
# THE EXACT SOLUTION
# ---------------------------------------------------------------------------
# Slug of mass M released at x = 0, t = 0, in an unbounded channel. Taking
# the Laplace transform in t (s = Laplace variable) and eliminating Cs gives
#
#   D C'' - v C' - kappa(s) C = -(M/A) delta(x),
#   kappa(s) = s + lambda + alpha - alpha*gamma / (s + gamma + lambda_s),
#   gamma    = alpha * A / As   (storage-side exchange rate, 1/s)
#
# i.e. the plain advection-dispersion (ADE) problem with s replaced by
# kappa(s). Inverting that substitution term by term gives an exact
# time-domain convolution ("time subordination": a molecule spends time tau
# moving in the channel, and while it moves it makes a Poisson number of
# trips into storage, each of exponential duration):
#
#   C(x,t) = exp(-(alpha+lambda) t) * C_ade(x,t)
#          + int_0^t C_ade(x,tau) exp(-(alpha+lambda) tau)
#                    * exp(-(gamma+lambda_s) u) * sqrt(alpha*gamma*tau/u)
#                    * I1(2 sqrt(alpha*gamma*tau*u)) dtau ,     u = t - tau
#
#   Cs(x,t) = int_0^t C_ade(x,tau) exp(-(alpha+lambda) tau)
#                    * gamma * exp(-(gamma+lambda_s) u)
#                    * I0(2 sqrt(alpha*gamma*tau*u)) dtau
#
#   C_ade(x,tau) = M / (A sqrt(4 pi D tau)) * exp(-(x - v tau)^2 / (4 D tau))
#
# I0, I1: modified Bessel functions of the first kind. The first term of
# C is the fraction of the slug that has never entered storage; the
# integral is everything that has. This is the classical mobile-immobile
# solution (cf. De Smedt & Wierenga 1979; Goltz & Roberts 1986; De Smedt
# 2005 for the stream TSM) extended with first-order losses in each zone.
#
# The only numerical step is a one-dimensional integral in tau, over the
# short window where C_ade is non-negligible, done by the trapezoid rule on
# a node spacing set by the width of C_ade itself -- so accuracy does not
# depend on reach length, velocity or D the way a spatial grid does.
#
# Checks (tsm_selftest() at the bottom of this file):
#   * alpha = 0 reproduces the closed-form ADE solution exactly;
#   * zeroth moment:  Q * int C dt = M  (lambda = lambda_s = 0);
#   * first moment:   t_mean = (x/v + 2D/v^2) (1 + As/A), the exact
#                     centroid of this model;
#   * second moment:  variance equals the exact second cumulant of the
#                     Laplace solution (i.e. the SHAPE is right, not just
#                     mass and timing);
#   * with uptake, zeroth moment matches the exact Laplace value at s = 0.
#
# Boundary conditions: unbounded channel (no upstream/downstream walls),
# the standard form for slug tests. NOTE: the old solver cannot be used as
# a reference here. Besides its numerical dispersion, it placed the slug in
# cell 1 next to a C = 0 upstream boundary, so dispersion carried part of
# the injected mass OUT through that boundary: ~13% of the mass lost at
# D = 0.05 m2/s with 40 cells, and more as the grid is refined (the loss
# rate scales with D/dx). It was invisible in the pipeline fits only
# because D had been pushed down to 1e-4, where the loss is negligible.
#
# Concentration reported is the RESIDENT (in-channel) concentration at x,
# same quantity as before.
################################################################################

suppressPackageStartupMessages({
  library(deSolve)   # only needed by simulate_tsm_fv() (legacy/comparison)
})

#' Closed-form ADE solution for an instantaneous slug in an unbounded
#' channel, on the log scale (avoids underflow far from the peak)
log_c_ade <- function(x, tau, v, D, A, mass) {
  log(mass / (A * sqrt(4 * pi * D * tau))) - (x - v * tau)^2 / (4 * D * tau)
}

#' Window of tau (s) outside which C_ade(x, tau) is below exp(-drop) of its
#' maximum, i.e. contributes nothing to the solution
ade_support <- function(x, v, D, drop = 30) {
  t_adv <- x / v
  grid  <- exp(seq(log(t_adv * 1e-4), log(t_adv * 1e3), length.out = 4000))
  lc    <- log_c_ade(x, grid, v, D, A = 1, mass = 1)
  keep  <- which(lc >= max(lc) - drop)
  lo <- grid[max(1, min(keep) - 1)]
  hi <- grid[min(length(grid), max(keep) + 1)]
  # width of the ADE pulse: sets the node spacing
  sigma <- sqrt(2 * D * x / v^3)
  list(lo = lo, hi = hi, sigma = sigma)
}

#' Simulate a slug (instantaneous) addition through a transient-storage reach
#' -- EXACT solution, see file header.
#'
#' @param L        reach length (m)
#' @param Q        discharge (m3/s)
#' @param A        main-channel cross-sectional area (m2). If NULL, computed
#'                 as Q/v.
#' @param v        mean velocity (m/s). If NULL, computed as Q/A.
#' @param D        dispersion coefficient (m2/s), > 0
#' @param alpha    storage-zone exchange coefficient (1/s); 0 = no storage
#' @param As       storage-zone cross-sectional area (m2); 0 = no storage
#' @param lambda   main-channel first-order uptake coefficient (1/s); 0 for a
#'                 conservative tracer
#' @param lambda_s storage-zone first-order uptake coefficient (1/s)
#' @param mass     mass injected, in units consistent with the desired output
#'                 concentration: grams for mg/L, mg for ug/L
#'                 (conc = mass / volume, 1 m3 = 1000 L, mg/L = g/m3).
#' @param times    output times (s). Any order; need not include 0.
#' @param n_cells  IGNORED. Kept so every existing caller (fit_hydraulics,
#'                 fit_uptake, partition_uptake, plot_btc_fit, the review
#'                 app) works unchanged. There is no spatial grid any more.
#' @param sample_x distance from the injection point to the sampling station
#'                 (default = L, i.e. the outlet)
#' @param storage_conc if TRUE, also compute Cs (costs ~2x). Default FALSE:
#'                 no caller in the pipeline uses Cs, so it is returned as NA.
#' @param nodes_per_sigma node density of the tau integral, per smallest
#'                 time scale of the integrand (ADE pulse width or storage-
#'                 time spread). 10 gives relative errors < 5e-4 of the
#'                 peak (tested); raise it for reference runs.
#'
#' @return data.frame(time, C, Cs) at the sampling station, one row per
#'   element of `times`, in the same order
simulate_tsm <- function(L, Q, A = NULL, v = NULL, D, alpha, As,
                         lambda = 0, lambda_s = 0, mass, times,
                         n_cells = 40, sample_x = L,
                         storage_conc = FALSE, nodes_per_sigma = 10) {
  
  if (is.null(A) && is.null(v)) stop("supply either A or v")
  if (is.null(A)) A <- Q / v
  if (is.null(v)) v <- Q / A
  
  D <- unname(D); alpha <- unname(alpha); As <- unname(As)
  lambda <- unname(lambda); lambda_s <- unname(lambda_s)
  x <- sample_x
  times <- as.numeric(times)
  
  if (!is.finite(D) || D <= 0) stop("D must be > 0")
  if (!is.finite(v) || v <= 0) stop("v must be > 0")
  if (!is.finite(x) || x <= 0) stop("sample_x must be > 0")
  if (any(c(alpha, As, lambda, lambda_s) < 0)) stop("alpha, As, lambda, lambda_s must be >= 0")
  
  n  <- length(times)
  C  <- numeric(n)
  Cs <- if (storage_conc) numeric(n) else rep(NA_real_, n)
  pos <- which(is.finite(times) & times > 0)
  if (length(pos) == 0) return(data.frame(time = times, C = C, Cs = Cs))
  
  storage <- alpha > 0 && As > 0
  
  # ---- no storage zone: closed-form ADE with first-order loss ---------------
  if (!storage) {
    tt <- times[pos]
    C[pos] <- exp(log_c_ade(x, tt, v, D, A, mass) - lambda * tt)
    if (storage_conc) Cs[pos] <- 0
    return(data.frame(time = times, C = C, Cs = Cs))
  }
  
  # ---- with storage: exact convolution (file header) ------------------------
  gam   <- alpha * A / As          # 1/s
  gam_s <- gam + lambda_s
  ag    <- alpha * gam
  
  sup <- ade_support(x, v, D)
  if (max(times[pos]) <= sup$lo) {   # every output time precedes the pulse
    return(data.frame(time = times, C = C, Cs = Cs))
  }
  
  # Node spacing must resolve BOTH scales in the integrand:
  #   * the ADE pulse, width sigma = sqrt(2 D x / v^3);
  #   * the storage kernel, whose finest scale is one exchange trip,
  #     1 / (gamma + lambda_s). With a small storage zone (As/A -> 0) gamma
  #     is large and this is seconds, far below sigma.
  # (2026-09 fix: the first version sized nodes on sigma only; for
  # RA_20231005_downstream -- As/A = 0.0017, 1/gamma ~ 4 s, node spacing
  # ~100 s -- the kernel was aliased and the curve oscillated.)
  for (i in pos) {
    t <- times[i]
    base_t <- log_c_ade(x, t, v, D, A, mass) - (alpha + lambda) * t
    direct <- exp(base_t)            # fraction never stored
    
    # window of tau that contributes: inside the ADE support AND within the
    # storage-time range reachable by time t (Poisson number of trips,
    # each exponential; generous upper quantile)
    m     <- alpha * t                       # expected number of trips
    n_hi  <- m + 10 * sqrt(m) + 10
    n_lo  <- max(0, m - 10 * sqrt(m) - 10)
    u_max <- (n_hi + 10 * sqrt(n_hi) + 10) / gam_s
    u_min <- if (n_lo > 0) max(0, (n_lo - 10 * sqrt(n_lo)) / gam_s) else 0
    # kernel scale: one trip (1/gamma) when few trips, the spread of the
    # total storage time, sqrt(2m)/gamma, when many (2026-09 speed fix: using
    # 1/gamma always made fast-exchange cases need ~20,000 nodes per time)
    w     <- max(1, sqrt(2 * m)) / gam_s
    h_max <- min(sup$sigma, w) / nodes_per_sigma
    a <- max(sup$lo, t - u_max)
    b <- min(sup$hi, t - u_min)
    if (b <= a) {
      C[i] <- direct
      if (storage_conc) Cs[i] <- 0
      next
    }
    n_tau <- max(50L, min(as.integer(ceiling((b - a) / h_max)) + 1L, 20000L))
    tau <- seq(a, b, length.out = n_tau)
    h   <- tau[2] - tau[1]
    u   <- t - tau
    z   <- 2 * sqrt(ag * tau * u)
    # z <= alpha*tau + gamma*u (AM-GM), so this exponent never overflows
    e   <- exp(log_c_ade(x, tau, v, D, A, mass) - (alpha + lambda) * tau - gam_s * u + z)
    
    at0 <- u <= 0                     # node at tau = t: use the analytic limit
    g1 <- numeric(n_tau)
    g1[!at0] <- e[!at0] * sqrt(ag * tau[!at0] / u[!at0]) *
      besselI(z[!at0], 1, expon.scaled = TRUE)
    g1[at0]  <- exp(base_t) * ag * t
    C[i] <- direct + h * (sum(g1) - 0.5 * (g1[1] + g1[n_tau]))
    
    if (storage_conc) {
      g0 <- gam * e * besselI(z, 0, expon.scaled = TRUE)
      Cs[i] <- h * (sum(g0) - 0.5 * (g0[1] + g0[n_tau]))
    }
  }
  
  data.frame(time = times, C = C, Cs = Cs)
}

################################################################################
# Self-test: run tsm_selftest() after any change to simulate_tsm().
# Returns a data.frame of checks; stops if any fails.
################################################################################
tsm_selftest <- function(L = 162, Q = 0.3833, v = 0.0817, mass = 7632,
                         alpha = 7.8832e-4, As_ratio = 0.156, tol = 1e-3) {
  A  <- Q / v
  As <- As_ratio * A
  trap <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
  out <- list()
  for (D in c(1e-3, 0.05, 0.5)) {
    # 1. alpha = 0 vs closed form
    tt  <- seq(1, 3 * L / v, length.out = 3000)
    num <- simulate_tsm(L, Q, A, D = D, alpha = 0, As = As, mass = mass, times = tt)$C
    ref <- mass / (A * sqrt(4 * pi * D * tt)) * exp(-(L - v * tt)^2 / (4 * D * tt))
    out[[length(out) + 1]] <- data.frame(D = D, check = "alpha=0 equals closed-form ADE",
                                         value = max(abs(num - ref)) / max(ref), target = 0)
    
    # 2-3. mass and centroid with storage (long, fine output grid)
    sig <- sqrt(2 * D * L / v^3)
    t_mean <- (L / v + 2 * D / v^2) * (1 + As_ratio)
    tt <- sort(unique(c(seq(1, 12 * t_mean, length.out = 20000),
                        seq(max(1, L / v - 15 * sig), L / v + 15 * sig, length.out = 4000))))
    cc <- simulate_tsm(L, Q, A, D = D, alpha = alpha, As = As, mass = mass, times = tt)$C
    m0 <- trap(tt, cc)
    out[[length(out) + 1]] <- data.frame(D = D, check = "mass recovery Q*int(C)/M",
                                         value = Q * m0 / mass, target = 1)
    m1 <- trap(tt, tt * cc) / m0
    out[[length(out) + 1]] <- data.frame(D = D, check = "centroid / exact (x/v+2D/v^2)(1+As/A)",
                                         value = m1 / t_mean, target = 1)
    
    # variance vs exact second cumulant of the Laplace solution
    gam0 <- alpha * A / As
    lnC  <- function(s) {
      k <- s + alpha - alpha * gam0 / (s + gam0)
      w <- sqrt(v^2 + 4 * D * k)
      L * (v - w) / (2 * D) - log(w)
    }
    hs <- 0.01 / t_mean
    var_exact <- (lnC(hs) - 2 * lnC(0) + lnC(-hs)) / hs^2
    out[[length(out) + 1]] <- data.frame(D = D, check = "variance / exact Laplace cumulant",
                                         value = trap(tt, (tt - m1)^2 * cc) / m0 / var_exact,
                                         target = 1)
    
    # 4. with uptake: zeroth moment vs exact Laplace value at s = 0
    lam <- 2e-4; lam_s <- 1e-3
    gam <- alpha * A / As
    kap0 <- lam + alpha - alpha * gam / (gam + lam_s)
    w0   <- sqrt(v^2 + 4 * D * kap0)
    m0_exact <- (mass / A) * exp(L * (v - w0) / (2 * D)) / w0
    cu <- simulate_tsm(L, Q, A, D = D, alpha = alpha, As = As, lambda = lam,
                       lambda_s = lam_s, mass = mass, times = tt)$C
    out[[length(out) + 1]] <- data.frame(D = D, check = "uptake: int(C) / exact Laplace M0",
                                         value = trap(tt, cu) / m0_exact, target = 1)
  }
  # 5. fast exchange (tiny storage zone, 1/gamma of seconds): the curve
  #    must be smooth and converged -- refining the nodes 4x changes nothing
  L2 <- 135; Q2 <- 0.210; v2 <- 0.0199; A2 <- Q2 / v2
  tt <- seq(0, 10200, by = 10)
  for (asr in c(0.0017, 0.02, 0.5)) {
    c1 <- simulate_tsm(L2, Q2, A2, D = 0.2, alpha = 4e-4, As = asr * A2, mass = 5947, times = tt)$C
    c4 <- simulate_tsm(L2, Q2, A2, D = 0.2, alpha = 4e-4, As = asr * A2, mass = 5947, times = tt,
                       nodes_per_sigma = 100)$C
    out[[length(out) + 1]] <- data.frame(D = 0.2, check = sprintf("converged, As/A = %g (4x nodes)", asr),
                                         value = max(abs(c1 - c4)) / max(c4), target = 0)
  }
  res <- do.call(rbind, out)
  res$ok <- abs(res$value - res$target) < tol
  print(res, row.names = FALSE)
  if (!all(res$ok)) stop("simulate_tsm() self-test FAILED")
  invisible(res)
}

################################################################################
# LEGACY finite-volume solver -- kept for comparison only (diag_tsm_shape.R).
# Do NOT use for fitting: first-order upwind advection adds numerical
# dispersion ~ v*dx/2, and the slug-in-cell-1 + C = 0 upstream boundary
# leaks mass upstream at a rate ~ D/dx (see header).
################################################################################

#' Build a discretisation grid for a reach
tsm_grid <- function(L, n_cells = 40) {
  n_cells <- max(10, min(round(n_cells), 120))
  dx <- L / n_cells
  list(n = n_cells, dx = dx, x_centers = (seq_len(n_cells) - 0.5) * dx)
}

#' Right-hand side of the TSM ODE system (legacy)
tsm_derivs <- function(t, y, parms) {
  n  <- parms$n
  dx <- parms$dx
  v  <- parms$v
  D  <- parms$D
  alpha <- parms$alpha
  A  <- parms$A
  As <- parms$As
  lambda   <- parms$lambda
  lambda_s <- parms$lambda_s
  
  C  <- y[1:n]
  Cs <- y[(n + 1):(2 * n)]
  
  C_up   <- c(0, C[1:(n - 1)])
  C_down <- c(C[2:n], C[n])
  
  advection  <- -v * (C - C_up) / dx
  dispersion <- D * (C_down - 2 * C + C_up) / dx^2
  exchange   <- alpha * (Cs - C)
  
  dC  <- advection + dispersion + exchange - lambda * C
  dCs <- alpha * (A / As) * (C - Cs) - lambda_s * Cs
  
  list(c(dC, dCs))
}

#' Legacy finite-volume slug simulation (the pre-2026-09 simulate_tsm()).
#' Same arguments and return value as simulate_tsm().
simulate_tsm_fv <- function(L, Q, A = NULL, v = NULL, D, alpha, As,
                            lambda = 0, lambda_s = 0, mass, times,
                            n_cells = 40, sample_x = L) {
  
  if (is.null(A) && is.null(v)) stop("supply either A or v")
  if (is.null(A)) A <- Q / v
  if (is.null(v)) v <- Q / A
  
  grid <- tsm_grid(L, n_cells)
  n <- grid$n
  dx <- grid$dx
  
  sample_cell <- max(1, min(n, round(sample_x / dx + 0.5)))
  
  y0 <- rep(0, 2 * n)
  y0[1] <- mass / (A * dx)
  
  parms <- list(n = n, dx = dx, v = v, D = D, alpha = alpha,
                A = A, As = As, lambda = lambda, lambda_s = lambda_s)
  
  prepended_zero <- times[1] != 0
  if (prepended_zero) times <- c(0, times)
  
  out <- suppressWarnings(
    ode(y = y0, times = times, func = tsm_derivs, parms = parms, method = "lsoda")
  )
  out <- as.data.frame(out)
  if (prepended_zero) out <- out[-1, , drop = FALSE]
  
  data.frame(
    time = out$time,
    C    = out[[1 + sample_cell]],
    Cs   = out[[1 + n + sample_cell]]
  )
}