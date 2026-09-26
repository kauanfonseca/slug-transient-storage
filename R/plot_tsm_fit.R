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

#' Overview of every current fit, one page per stream: small multiples of
#' NaCl (Stage 1), NH4-N and SRP (Stage 2) for each event, each panel on
#' its own scales. Built only from saved outputs (no refitting), so it can
#' be run from the CSVs alone:
#'   plot_fits_grid(read_csv("data_derived/tsm_uptake_summary.csv"),
#'                  read_csv("data_derived/tsm_hydraulics_summary.csv"),
#'                  events, btc_conservative, master_tsm)
#' Nutrient panels show exactly what Stage 2 saw: grabs after
#' apply_nutrient_exclusions() and the stored nutrient_time_shift_s, gap-fill
#' zeros as hollow points, and the excluded grabs as red crosses at their
#' shifted times. NaCl panels use the mass the Stage-1 fit used
#' (nacl_mass_fit_g for logger_rescaled; for the "shared" RA_20230906_N, its
#' own logger recovery, since its hydraulics are P's). Raw correction only.
#'
#' @return a named list of ggplots (one per stream); print() them
plot_fits_grid <- function(summary_tbl, hyd_tbl, events, btc_conservative, master_tsm,
                           n_cells = 40, onset_frac = 0.05) {
  tracer_lv <- c("NaCl (mg/L)", "NH4-N (ug/L)", "SRP (ug/L)")
  obs <- list(); mod <- list(); lab <- list()

  # ---- NaCl, Stage 1 ----
  for (i in seq_len(nrow(hyd_tbl))) {
    h <- hyd_tbl[i, ]; eid <- h$event_id
    e <- events[events$event_id == eid, ]
    src <- h$conservative_source
    if (src %in% c("probe_grab", "borrowed")) {
      sub <- unique(master_tsm[master_tsm$event_id == eid & !is.na(master_tsm$nacl_mgL_grab) &
                                 master_tsm$time_since_release_s >= 0,
                               c("time_since_release_s", "nacl_mgL_grab")])
      names(sub) <- c("t", "c")
    } else {
      sub <- btc_conservative[btc_conservative$event_id == eid & btc_conservative$time_since_release_s >= 0,
                              c("time_since_release_s", "nacl_mgL")]
      names(sub) <- c("t", "c")
    }
    sub <- sub[order(sub$t), ]; sub$c <- pmax(sub$c, 0)
    if (nrow(sub) == 0) next
    mass <- if (is.finite(h$nacl_mass_fit_g %||% NA_real_)) h$nacl_mass_fit_g else e$nacl_mass_g
    if (identical(src, "shared")) {   # P's hydraulics on N's own logger: use N's own recovered mass
      zp <- zero_pre_arrival(sub$t, sub$c, onset_frac)
      mass <- h$Q_m3s * sum(diff(sub$t) * (head(zp$value, -1) + tail(zp$value, -1)) / 2)
    }
    tt <- seq(0, max(sub$t), length.out = 300)
    C <- simulate_tsm(L = h$L_m, Q = h$Q_m3s, A = h$A_m2, D = h$D_m2s, alpha = h$alpha_1s, As = h$As_m2,
                      lambda = 0, lambda_s = 0, mass = mass, times = tt, n_cells = n_cells)$C
    tag <- switch(src, logger_rescaled = "logger, rescaled", shared = paste("shared:", h$hydraulics_borrowed_from),
                  borrowed = paste("borrowed:", h$hydraulics_borrowed_from), src)
    obs[[length(obs) + 1]] <- data.frame(event_id = eid, tracer = tracer_lv[1], t = sub$t / 60, c = sub$c, kind = "kept")
    mod[[length(mod) + 1]] <- data.frame(event_id = eid, tracer = tracer_lv[1], t = tt / 60, c = C)
    lab[[length(lab) + 1]] <- data.frame(event_id = eid, tracer = tracer_lv[1],
      txt = sprintf("%s | %s\nRMSE %s", h$hydraulic_model, tag,
                    ifelse(is.na(h$hydraulic_rmse), "-", sprintf("%.2f", h$hydraulic_rmse))))
  }

  # ---- nutrients, Stage 2 ----
  ok <- summary_tbl[summary_tbl$fit_status == "ok" & summary_tbl$conc_col == "conc_corr_ugL", ]
  for (i in seq_len(nrow(ok))) {
    r <- ok[i, ]; eid <- r$event_id
    h <- hyd_tbl[hyd_tbl$event_id == eid, ]
    trc <- if (r$solute == "NH4-N") tracer_lv[2] else tracer_lv[3]
    nut <- master_tsm[master_tsm$event_id == eid & master_tsm$solute == r$solute & !is.na(master_tsm$conc_corr_ugL), ]
    kept <- apply_nutrient_exclusions(eid, r$solute, "conc_corr_ugL", nut)
    gone <- nut[!nut$time_since_release_s %in% kept$time_since_release_s, ]
    sh <- r$nutrient_time_shift_s %||% 0
    kt <- kept$time_since_release_s - sh; kc <- pmax(kept$conc_corr_ugL, 0)
    kc <- kc[kt >= 0]; kt <- kt[kt >= 0]
    gf <- fill_pre_arrival_gap(kt, kc)
    kind <- ifelse(gf$kind == "gap_fill", "gap fill", "kept")
    gt <- gone$time_since_release_s - sh
    tt <- seq(0, max(c(gf$time, gt[gt >= 0])), length.out = 300)
    C <- simulate_tsm(L = h$L_m, Q = h$Q_m3s, A = r$A_m2, D = r$D_m2s, alpha = r$alpha_1s, As = r$As_m2,
                      lambda = r$lambda_1s, lambda_s = r$lambda_s_1s, mass = r$mass_injected_mg,
                      times = tt, n_cells = n_cells)$C
    obs[[length(obs) + 1]] <- data.frame(event_id = eid, tracer = trc, t = gf$time / 60, c = gf$value, kind = kind)
    if (nrow(gone)) obs[[length(obs) + 1]] <- data.frame(event_id = eid, tracer = trc, t = gt[gt >= 0] / 60,
                                                         c = pmax(gone$conc_corr_ugL[gt >= 0], 0), kind = "excluded")
    mod[[length(mod) + 1]] <- data.frame(event_id = eid, tracer = trc, t = tt / 60, c = C)
    st <- if (isTRUE(r$uptake_significant) && is.finite(r$pct_storagezone_min %||% NA_real_))
      sprintf("\nstorage %.0f-%.0f%% of uptake", r$pct_storagezone_min, r$pct_storagezone_max) else ""
    lab[[length(lab) + 1]] <- data.frame(event_id = eid, tracer = trc,
      txt = sprintf("shift %+ds | RMSE %.2f\nuptake %.0f%%%s", as.integer(sh), r$uptake_rmse, r$pct_total_uptake, st))
  }

  obs <- do.call(rbind, obs); mod <- do.call(rbind, mod); lab <- do.call(rbind, lab)
  for (d in c("obs", "mod", "lab")) {
    x <- get(d); x$stream <- sub("_.*", "", x$event_id)
    x$tracer <- factor(x$tracer, levels = tracer_lv); assign(d, x)
  }
  # excluded spikes (up to 630 ug/L) would flatten the fitted curve; cap them
  # at the panel top so the fit stays readable, and mark them as off-scale
  ymax <- aggregate(c ~ event_id + tracer, data = rbind(obs[obs$kind != "excluded", c("event_id", "tracer", "c")],
                                                         mod[, c("event_id", "tracer", "c")]), FUN = max)
  obs <- merge(obs, setNames(ymax, c("event_id", "tracer", "ymax")), by = c("event_id", "tracer"), all.x = TRUE)
  off <- obs$kind == "excluded" & obs$c > 1.25 * obs$ymax
  obs$c[off] <- 1.25 * obs$ymax[off]; obs$kind[off] <- "excluded (off-scale)"
  obs$kind <- factor(obs$kind, levels = c("kept", "gap fill", "excluded", "excluded (off-scale)"))

  plots <- lapply(split(unique(obs$stream), unique(obs$stream)), function(st) {
    o <- obs[obs$stream == st, ]; m <- mod[mod$stream == st, ]; l <- lab[lab$stream == st, ]
    ggplot() +
      geom_line(data = m, aes(t, c), color = "#1b6ca8", linewidth = 0.7) +
      geom_point(data = o, aes(t, c, shape = kind, color = kind), size = 1.3, stroke = 0.6) +
      geom_text(data = l, aes(x = -Inf, y = Inf, label = txt), hjust = -0.04, vjust = 1.15,
                size = 2.4, color = "grey25", lineheight = 0.9) +
      scale_shape_manual(values = c(kept = 16, `gap fill` = 1, excluded = 4, `excluded (off-scale)` = 2), drop = FALSE) +
      scale_color_manual(values = c(kept = "grey25", `gap fill` = "grey55", excluded = "#c0392b",
                                    `excluded (off-scale)` = "#c0392b"), drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.35))) +
      facet_wrap(~ event_id + tracer, scales = "free", ncol = 3,
                 labeller = labeller(.multi_line = FALSE)) +
      labs(title = paste("Stream", st, "- current TSM fits (line = model)"),
           x = "time since release (min; nutrients clock-shifted)", y = NULL, shape = NULL, color = NULL) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "top", panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey92"), strip.text = element_text(face = "bold", size = 7.5))
  })
  plots
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