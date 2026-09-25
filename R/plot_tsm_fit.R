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
#' @param obs_kind optional parallel character vector ("observed" /
#'   "gap_fill") from fill_pre_arrival_gap() -- when supplied, gap-fill
#'   points (synthetic pre-arrival background, added because sparse grab
#'   series only get logged once the signal changes -- see
#'   fill_pre_arrival_gap() in tsm_calibrate.R) are drawn as hollow markers
#'   distinct from the real observed points, so the fill is always visible,
#'   never hidden in the fit plot
#' @param L,Q,A,D,alpha,As,lambda,lambda_s,mass as in simulate_tsm()
#' @param title plot title
#' @param unit y-axis unit label
plot_btc_fit <- function(obs_time, obs_conc, L, Q, A, D, alpha, As,
                         lambda = 0, lambda_s = 0, mass, title = "", unit = "mg/L",
                         n_cells = 40, obs_kind = NULL) {
  sim_times <- seq(0, max(obs_time), length.out = 400)
  sim <- simulate_tsm(L = L, Q = Q, A = A, D = D, alpha = alpha, As = As,
                      lambda = lambda, lambda_s = lambda_s, mass = mass,
                      times = sim_times, n_cells = n_cells)
  if (is.null(obs_kind)) obs_kind <- rep("observed", length(obs_time))
  obs_df <- data.frame(time = obs_time, conc = obs_conc, kind = obs_kind)
  ggplot() +
    geom_point(data = obs_df, aes(x = time, y = conc, shape = kind, color = kind), size = 1.5, alpha = 0.75) +
    scale_shape_manual(values = c(observed = 16, gap_fill = 1), guide = if (all(obs_kind == "observed")) "none" else "legend") +
    scale_color_manual(values = c(observed = "grey30", gap_fill = "#c0392b"), guide = if (all(obs_kind == "observed")) "none" else "legend") +
    geom_line(data = sim, aes(x = time, y = C), color = "#1b6ca8", linewidth = 0.9) +
    labs(title = title, x = "time since release (s)", y = paste0("concentration (", unit, ")"),
         shape = "", color = "") +
    theme_minimal(base_size = 12)
}

#' Print the conservative-tracer (NaCl) fit for EVERY event with fitted or
#' borrowed hydraulics -- one figure per event, printed in a loop so each
#' shows up as its own plot when run in an R Markdown chunk (or in the
#' console/RStudio plot pane if run interactively).
#'
#' @param hyd_fits the `attr(summary_tbl, "hydraulic_fits")` list from
#'   run_tsm_uptake_all()
#' @param events,btc_conservative,master_tsm the same tables passed into the
#'   pipeline -- master_tsm supplies the hand-held-probe NaCl series
#'   (`nacl_mgL_grab`) for events whose `conservative_source` is
#'   `"probe_grab"` (see fit_all_hydraulics() for why some events use that
#'   instead of the logger)
plot_all_hydraulic_fits <- function(hyd_fits, events, btc_conservative, master_tsm, n_cells = 40) {
  for (eid in names(hyd_fits)) {
    hf <- hyd_fits[[eid]]
    e <- events[events$event_id == eid, ]
    
    if (identical(hf$conservative_source, "probe_grab")) {
      sub <- master_tsm[master_tsm$event_id == eid & !is.na(master_tsm$nacl_mgL_grab) &
                          master_tsm$time_since_release_s >= 0, c("time_since_release_s", "nacl_mgL_grab")]
      sub <- unique(sub)
      names(sub)[names(sub) == "nacl_mgL_grab"] <- "nacl_mgL"
    } else {
      sub <- btc_conservative[btc_conservative$event_id == eid & btc_conservative$time_since_release_s >= 0, ]
    }
    sub <- sub[order(sub$time_since_release_s), ]
    if (nrow(sub) == 0) {
      message(eid, ": no own conservative-tracer series (hydraulics borrowed) -- skipped")
      next
    }
    borrowed <- if (!is.null(hf$hydraulics_borrowed_from)) {
      paste0(" [hydraulics borrowed from ", hf$hydraulics_borrowed_from, "]")
    } else ""
    src_tag <- switch(hf$conservative_source %||% "",
                      probe_grab = " [probe grabs]",
                      logger_rescaled = sprintf(" [logger shape, mass x %.2f]", hf$logger_recovery),
                      "")
    # re-derive the same pre-arrival gap fill used at fit time (probe grabs
    # only -- see fit_all_hydraulics()) purely for plotting, so the fill is
    # always visible (hollow red points) rather than hidden inside the fit
    gf <- if (identical(hf$conservative_source, "probe_grab") && isTRUE(hf$n_gap_filled > 0)) {
      fill_pre_arrival_gap(sub$time_since_release_s, pmax(sub$nacl_mgL, 0), dt = hf$gap_fill_dt_s)
    } else {
      list(time = sub$time_since_release_s, value = pmax(sub$nacl_mgL, 0),
           kind = rep("observed", nrow(sub)))
    }
    # the mass the fit actually used (rescaled for "logger_rescaled" events)
    mass_plot <- if (!is.null(hf$nacl_mass_fit_g)) hf$nacl_mass_fit_g else e$nacl_mass_g
    p <- plot_btc_fit(gf$time, gf$value,
                      hf$L, hf$Q, hf$A, hf$D, hf$alpha, hf$As, mass = mass_plot,
                      n_cells = n_cells, unit = "mg/L", obs_kind = gf$kind,
                      title = sprintf("%s -- NaCl fit (RMSE=%.3f)%s%s", eid,
                                      ifelse(is.na(hf$rmse), NA, hf$rmse), src_tag, borrowed))
    print(p)
  }
}

#' Print the nutrient-uptake fit for EVERY event x solute x correction
#' combination with fit_status == "ok" -- one figure each. Rebuilds the plot
#' straight from the fitted lambda/lambda_s and hydraulics already stored in
#' `summary_tbl`/`hyd_fits` (no re-fitting, so this is fast).
#'
#' Reapplies the SAME per-event/solute manual exclusions
#' (NUTRIENT_EXCLUDED_TIMES / NUTRIENT_VALUE_THRESHOLD, run_tsm_uptake.R)
#' and clock-offset shift (already found and stored as
#' `nutrient_time_shift_s` by run_one_uptake() -- not re-searched here) that
#' the fit itself used, so the plotted points always match what the
#' lambda/lambda_s/RMSE in `summary_tbl` were actually fit to. Before this,
#' the plot re-pulled the untouched raw grabs from master_tsm, so a
#' correctly-fit curve could still be shown against the old, uncorrected
#' points (2026-09-25).
#'
#' @param summary_tbl the tibble returned by run_tsm_uptake_all()
#' @param hyd_fits attr(summary_tbl, "hydraulic_fits")
#' @param master_tsm the master_tsm table (for the observed grab series)
plot_all_uptake_fits <- function(summary_tbl, hyd_fits, master_tsm, n_cells = 40) {
  ok <- summary_tbl[summary_tbl$fit_status == "ok", ]
  for (i in seq_len(nrow(ok))) {
    r <- ok[i, ]
    hf <- hyd_fits[[r$event_id]]
    nut <- master_tsm[master_tsm$event_id == r$event_id & master_tsm$solute == r$solute &
                        !is.na(master_tsm[[r$conc_col]]), ]
    nut <- nut[order(nut$time_since_release_s), ]
    # same hard-coded exclusions run_one_uptake() applied before fitting
    if (exists("apply_nutrient_exclusions")) {
      nut <- apply_nutrient_exclusions(r$event_id, r$solute, r$conc_col, nut)
    }
    # same clock-offset correction, using the shift ALREADY found and
    # stored by run_one_uptake() (nutrient_time_shift_s) -- not re-searched
    shift_s <- if (!is.null(r$nutrient_time_shift_s) && is.finite(r$nutrient_time_shift_s)) r$nutrient_time_shift_s else 0
    nut$time_since_release_s <- nut$time_since_release_s - shift_s
    nut <- nut[nut$time_since_release_s >= 0, ]
    nut <- nut[order(nut$time_since_release_s), ]
    flag <- if (isTRUE(r$lambda_at_bound) || isTRUE(r$lambda_s_at_bound)) "  [check identifiability]" else ""
    shift_tag <- if (shift_s != 0) sprintf(" [t-shift %+ds]", as.integer(shift_s)) else ""
    # re-derive the pre-arrival gap fill applied at fit time, purely for
    # plotting (nutrient grabs are always sparse -- see run_one_uptake())
    gf <- if (isTRUE(r$n_gap_filled_uptake > 0)) {
      fill_pre_arrival_gap(nut$time_since_release_s, pmax(nut[[r$conc_col]], 0), dt = r$gap_fill_dt_uptake_s)
    } else {
      list(time = nut$time_since_release_s, value = pmax(nut[[r$conc_col]], 0),
           kind = rep("observed", nrow(nut)))
    }
    p <- plot_btc_fit(gf$time, gf$value,
                      hf$L, hf$Q, hf$A, hf$D, hf$alpha, hf$As,
                      lambda = r$lambda_1s, lambda_s = r$lambda_s_1s,
                      mass = r$mass_injected_mg, n_cells = n_cells, unit = "ug/L", obs_kind = gf$kind,
                      title = sprintf("%s -- %s (%s)%s%s", r$event_id, r$solute, r$correction, shift_tag, flag))
    print(p)
  }
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

`%||%` <- function(a, b) if (is.null(a)) b else a