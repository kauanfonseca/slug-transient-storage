################################################################################
# tsm_partition.R
#
# Runkel (2007) mass-balance partitioning of nutrient uptake between the main
# channel and the transient-storage (transitional/hyporheic) zone, plus the
# standard nutrient-spiralling uptake metrics (Sw, vf, U), computed from the
# calibrated TSM (tsm_model.R + tsm_calibrate.R) instead of OTIS-P.
#
# Partitioning logic (Runkel 2007, section "Time-series data sets", eqs 22-24,
# applied here at the single reach/outlet level rather than reach-by-reach
# along a multi-reach network -- see that paper's Green Creek example, Table 3,
# for the network version this simplifies):
#
#   Simulation A (conservative):        lambda = 0,       lambda_s = 0
#   Simulation B (fully reactive):      lambda = lambda,  lambda_s = lambda_s
#   Simulation C (main-channel uptake): lambda = lambda,  lambda_s = 0
#   Simulation D (storage-zone uptake): lambda = 0,       lambda_s = lambda_s
#
#   mass_i = Q * integral_0^T C_i(t) dt        (mass passing the outlet)
#   total_uptake      = mass_A - mass_B
#   pct_mainchannel   = (mass_A - mass_C) / total_uptake   (raw, C alone overuptakes)
#   pct_storagezone   = (mass_A - mass_D) / total_uptake   (raw, D alone overuptakes)
#   -> renormalised so pct_mainchannel + pct_storagezone == 100%, because
#      simulations C and D each attribute the *entire* channel/storage
#      interaction to a single removal pathway and so each individually
#      overestimates its own share (Runkel 2007, p.58, "simulations C and D
#      overestimate uptake").
#
# Biologically, for these experiments:
#   - the main-channel share is attributed to biotic assimilation in the
#     water column / epilithon (periphyton uptake), the dominant near-surface
#     process (Runkel 2007; Tromboni et al. 2017);
#   - the storage-zone share is the transient-storage / hyporheic-type
#     removal (microbial uptake, sorption, or -- for SRP specifically --
#     CaCO3 co-precipitation, see the circularity caveat in the project
#     handoff Section 4.6/7). For SRP, running this on `conc_corr_coprec_ugL`
#     (background AND co-precipitation corrected) instead of `conc_corr_ugL`
#     is how the abiotic co-precipitation fraction gets discounted so that
#     what remains approximates biotic assimilation only.
################################################################################

#' Simple trapezoidal integral
trapz <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)

#' Choose a simulation horizon long enough for the BTC to fully return to
#' background, so the mass-balance integrals close (avoids the "NA%"
#' truncation problem noted in the project's own earlier prototype).
#'
#' @param L,Q,A reach length (m), discharge (m3/s), channel area (m2)
#' @param alpha,As storage exchange coefficient (1/s), storage area (m2)
#' @param n_points number of output times used for the integration
choose_horizon <- function(L, Q, A, alpha, As, n_points = 600, cap_s = 86400) {
  v <- Q / A
  travel_time <- L / v
  t_sto <- As / (alpha * A)          # mean storage-zone residence time (Thackston & Schnelle 1970)
  T_max <- min(cap_s, max(8 * travel_time, 8 * t_sto, 3600))
  list(times = seq(0, T_max, length.out = n_points), T_max = T_max, t_sto = t_sto,
       travel_time = travel_time)
}

#' Runkel (2007) four-simulation mass-balance partition of uptake between the
#' main channel and the transient-storage zone
#'
#' @param L,Q,A,D,alpha,As TSM parameters (see tsm_model.R)
#' @param lambda,lambda_s fitted first-order uptake coefficients (1/s)
#' @param mass nutrient mass injected (consistent units, see tsm_model.R)
#' @param n_cells spatial resolution
#' @return list with mass_A..D, pct_mainchannel, pct_storagezone (normalised
#'   to 100%), pct_total_uptake (mass_A-mass_B as % of mass_A), and the
#'   simulation horizon used
partition_uptake <- function(L, Q, A, D, alpha, As, lambda, lambda_s, mass,
                              n_cells = 40) {

  D <- unname(D); alpha <- unname(alpha); As <- unname(As)
  lambda <- unname(lambda); lambda_s <- unname(lambda_s)

  hz <- choose_horizon(L, Q, A, alpha, As)
  times <- hz$times

  run_mass <- function(lam, lams) {
    sim <- simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = alpha, As = As,
                         lambda = lam, lambda_s = lams, mass = mass,
                         times = times, n_cells = n_cells, sample_x = L)
    Q * trapz(sim$time, pmax(sim$C, 0))
  }

  mass_A <- run_mass(0, 0)
  mass_B <- run_mass(lambda, lambda_s)
  mass_C <- run_mass(lambda, 0)
  mass_D <- run_mass(0, lambda_s)

  total_uptake <- mass_A - mass_B
  pct_mc_raw <- (mass_A - mass_C) / total_uptake * 100
  pct_sz_raw <- (mass_A - mass_D) / total_uptake * 100

  if (!is.finite(total_uptake) || total_uptake <= 0 ||
      !is.finite(pct_mc_raw) || !is.finite(pct_sz_raw) ||
      (pct_mc_raw + pct_sz_raw) == 0) {
    pct_mc <- NA_real_; pct_sz <- NA_real_
  } else {
    pct_mc <- pct_mc_raw / (pct_mc_raw + pct_sz_raw) * 100
    pct_sz <- pct_sz_raw / (pct_mc_raw + pct_sz_raw) * 100
  }

  list(mass_A = mass_A, mass_B = mass_B, mass_C = mass_C, mass_D = mass_D,
       pct_total_uptake = total_uptake / mass_A * 100,
       pct_mainchannel = pct_mc, pct_storagezone = pct_sz,
       horizon_s = hz$T_max, t_sto_s = hz$t_sto, travel_time_s = hz$travel_time)
}

#' Standard nutrient-spiralling uptake metrics (uptake length Sw, uptake
#' velocity vf, areal uptake rate U), reported separately for the main
#' channel and the storage zone (Runkel 2007 eqs 8, 25, 26), plus a
#' whole-reach aggregate Sw comparable to the conventional (OTIS-P / Sw
#' linear-regression) literature value (Runkel 2007 eq 18 effective storage
#' coefficient, combined via eq 12's Case-II form).
#'
#' @param v mean channel velocity (m/s)
#' @param depth_main main-channel mean depth (m) = A / width
#' @param depth_storage storage-zone "depth" (m) = As / width, a bookkeeping
#'   convention (Runkel 2007 eq 26) letting vf be computed the same way as
#'   for the channel; not a physical depth of the storage zone itself
#' @param lambda,lambda_s fitted coefficients (1/s)
#' @param alpha,A,As storage exchange coefficient and areas
#' @param Camb ambient (background) nutrient concentration, for U = vf * C
uptake_metrics <- function(v, depth_main, depth_storage, lambda, lambda_s,
                            alpha, A, As, Camb) {

  v <- unname(v); depth_main <- unname(depth_main); depth_storage <- unname(depth_storage)
  lambda <- unname(lambda); lambda_s <- unname(lambda_s)
  alpha <- unname(alpha); A <- unname(A); As <- unname(As); Camb <- unname(Camb)

  Sw_mc <- if (lambda > 0) v / lambda else NA_real_
  vf_mc <- lambda * depth_main
  U_mc  <- vf_mc * Camb

  vf_sz <- lambda_s * depth_storage
  U_sz  <- vf_sz * Camb

  # effective storage-zone coefficient projected onto the main channel
  # (Runkel 2007, eq. 18), and the resulting whole-reach aggregate Sw
  # (eq. 12, Case II form) -- the number directly comparable to a
  # conventional Sw obtained by linear regression on steady-state data
  lambda_eff_s <- (alpha * lambda_s * As) / (alpha * A + lambda_s * As)
  Sw_total <- if ((lambda + lambda_eff_s) > 0) v / (lambda + lambda_eff_s) else NA_real_
  vf_total <- v * depth_main / Sw_total
  U_total  <- vf_total * Camb

  list(Sw_mainchannel_m = Sw_mc, vf_mainchannel_mps = vf_mc, U_mainchannel = U_mc,
       vf_storagezone_mps = vf_sz, U_storagezone = U_sz,
       lambda_eff_storagezone = lambda_eff_s,
       Sw_total_m = Sw_total, vf_total_mps = vf_total, U_total = U_total)
}
