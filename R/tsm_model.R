################################################################################
# tsm_model.R
#
# One-dimensional transient-storage transport model (Bencala & Walters 1983;
# Runkel 1998, OTIS), solved directly in R by finite volumes + deSolve,
# without requiring the OTIS/OTIS-P executables.
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
# The reach is discretised into `n_cells` finite volumes of width dx = L/n_cells.
# A slug (instantaneous) addition is represented as an initial condition: all
# of the injected mass is placed, fully mixed, into the first cell. This is
# the standard treatment for tracer/nutrient SLUG tests (as opposed to a
# constant-rate injection, which would instead need a boundary condition) and
# avoids having to know/guess an injection duration.
#
# Boundary conditions:
#   upstream (x=0-):  C = 0  (ambient/background water flowing in from above
#                      the injection point, in background-corrected units)
#   downstream (x=L+): zero-gradient (transmissive; mass leaves only by
#                      advection, consistent with OTIS's own open outlet)
#
# This reproduces the physics of OTIS's transient-storage module (Runkel
# 1998) without needing OTIS itself, following the same rationale used by
# Bonanno et al. (2022) for building a fast, R-native alternative to OTIS-P
# when doing many parameter evaluations (there: for identifiability
# analysis; here: for batch fitting many streams).
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

#' Build a discretisation grid for a reach
#'
#' @param L reach length (m)
#' @param n_cells number of finite-volume cells (a single dx is used)
#' @return list(n, dx, x_centers)
tsm_grid <- function(L, n_cells = 40) {
  n_cells <- max(10, min(round(n_cells), 120))
  dx <- L / n_cells
  list(n = n_cells, dx = dx, x_centers = (seq_len(n_cells) - 0.5) * dx)
}

#' Right-hand side of the TSM ODE system
#'
#' State vector y = c(C[1..n], Cs[1..n])  (channel block, then storage block)
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

  # ghost cells: upstream Dirichlet C=0, downstream zero-gradient
  C_up   <- c(0, C[1:(n - 1)])
  C_down <- c(C[2:n], C[n])

  advection  <- -v * (C - C_up) / dx
  dispersion <- D * (C_down - 2 * C + C_up) / dx^2
  exchange   <- alpha * (Cs - C)

  dC  <- advection + dispersion + exchange - lambda * C
  dCs <- alpha * (A / As) * (C - Cs) - lambda_s * Cs

  list(c(dC, dCs))
}

#' Simulate a slug (instantaneous) addition through a transient-storage reach
#'
#' @param L        reach length (m)
#' @param Q        discharge (m3/s)
#' @param A        main-channel cross-sectional area (m2). If NULL, Q/velocity
#'                 must be supplied via `v`.
#' @param v        mean velocity (m/s). If NULL, computed as Q/A.
#' @param D        dispersion coefficient (m2/s)
#' @param alpha    storage-zone exchange coefficient (1/s)
#' @param As       storage-zone cross-sectional area (m2)
#' @param lambda   main-channel first-order uptake coefficient (1/s); 0 for a
#'                 conservative tracer
#' @param lambda_s storage-zone first-order uptake coefficient (1/s)
#' @param mass     mass injected, in units consistent with the desired output
#'                 concentration (e.g. supply grams for mg/L output, or mg for
#'                 ug/L output: conc = mass / volume with volume in m3 == L
#'                 of water, since 1 m3 = 1000 L and mg/L = g/m3).
#'                 Concretely: if you want C in mg/L, mass must be in grams;
#'                 if you want C in ug/L, mass must be in mg.
#' @param times    output times (s), must include 0
#' @param n_cells  spatial resolution (see tsm_grid)
#' @param sample_x distance from the injection point to the sampling station
#'                 (default = L, i.e. the outlet). Must be <= L.
#'
#' @return data.frame(time, C, Cs) at the sampling station
simulate_tsm <- function(L, Q, A = NULL, v = NULL, D, alpha, As,
                          lambda = 0, lambda_s = 0, mass, times,
                          n_cells = 40, sample_x = L) {

  if (is.null(A) && is.null(v)) stop("supply either A or v")
  if (is.null(A)) A <- Q / v
  if (is.null(v)) v <- Q / A

  grid <- tsm_grid(L, n_cells)
  n <- grid$n
  dx <- grid$dx

  # sampling cell = the cell whose centre is nearest sample_x
  sample_cell <- max(1, min(n, round(sample_x / dx + 0.5)))

  y0 <- rep(0, 2 * n)
  y0[1] <- mass / (A * dx)   # instantaneous, fully-mixed slug in cell 1

  parms <- list(n = n, dx = dx, v = v, D = D, alpha = alpha,
                A = A, As = As, lambda = lambda, lambda_s = lambda_s)

  if (times[1] != 0) times <- c(0, times)

  out <- suppressWarnings(
    ode(y = y0, times = times, func = tsm_derivs, parms = parms, method = "lsoda")
  )
  out <- as.data.frame(out)

  data.frame(
    time = out$time,
    C    = out[[1 + sample_cell]],
    Cs   = out[[1 + n + sample_cell]]
  )
}
