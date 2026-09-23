################################################################################
# diag_ra_logger.R
#
# Quick test: run the NEW Stage-1 method (exact solver, v_fixed vs v_fitted
# with AICc selection, pre-arrival zeroing -- exactly fit_all_hydraulics()
# from run_tsm_uptake.R) on the LOGGER curves of the RA events, which the
# pipeline currently fits from the hand-held-probe grabs instead
# (events$discharge_source == "probe").
#
# Three fits per event, side by side:
#   probe_grab       : what the pipeline does now (grab series)
#   logger           : logger curve, Q fixed from events (same as any logger event)
#   logger_rescaled  : logger curve, injected mass scaled to the logger's own
#                      recovered mass (Q * integral C dt). Diagnostic only:
#                      it lets the fit use the SHAPE of the logger curve even
#                      where its area disagrees with Q (earlier screen: logger
#                      recovery 36% for RA_20230906_N, 61% for _P, 0.7% for
#                      RA_20231005_downstream). If this variant fits well and
#                      "logger" does not, the logger shape is usable but its
#                      area is not (incomplete mixing at the logger, or a
#                      conversion issue) -- that is a data question, not a
#                      model question.
#
# Standalone: does not modify any pipeline file or output.
# Writes to data_derived/diagnostics/ra_logger/
#
# Run time: 3 variants x 2 candidate models x N events, ~1-2 min per fit.
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr)
  library(ggplot2); library(here)
})

source(here("R", "tsm_model.R"))
source(here("R", "tsm_calibrate.R"))
source(here("R", "tsm_partition.R"))
source(here("R", "run_tsm_uptake.R"))   # fit_all_hydraulics(), tsm_hydraulics_table()

# ---- settings ---------------------------------------------------------------
IDS   <- c("RA_20230906_N", "RA_20230906_P")
N_LHS <- 250      # same as the batch; lower (e.g. 80) for a quicker first look
OUT   <- here("data_derived", "diagnostics", "ra_logger")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- data -------------------------------------------------------------------
events <- read_csv(here("data_derived", "events.csv"), show_col_types = FALSE)
btc    <- read_csv(here("data_derived", "btc_conservative.csv"), show_col_types = FALSE)
mt     <- read_csv(here("data_derived", "master_tsm.csv"), show_col_types = FALSE)

ev <- events %>% filter(event_id %in% IDS)
missing_logger <- setdiff(IDS, unique(btc$event_id))
if (length(missing_logger)) warning("no logger series in btc_conservative for: ",
                                    paste(missing_logger, collapse = ", "))

# logger recovery on the same (zeroed, clamped) series the fit will see
logger_recovery <- btc %>%
  filter(event_id %in% IDS, time_since_release_s >= 0) %>%
  arrange(event_id, time_since_release_s) %>%
  distinct(event_id, time_since_release_s, .keep_all = TRUE) %>%
  group_by(event_id) %>%
  summarise(n_logger = n(),
            M0_mgLs  = {
              zp <- zero_pre_arrival(time_since_release_s, pmax(nacl_mgL, 0),
                                     TSM_METHOD$onset_frac)
              trapz(time_since_release_s, zp$value)
            },
            .groups = "drop") %>%
  left_join(ev %>% select(event_id, discharge_Ls, nacl_mass_g), by = "event_id") %>%
  mutate(logger_recovery = (discharge_Ls / 1000) * M0_mgLs / nacl_mass_g)

cat("\n=== logger recovery (Q * integral C dt / mass) ===\n")
print(logger_recovery, width = Inf)

# ---- the three variants -----------------------------------------------------
ev_probe  <- ev                                            # pipeline as-is
ev_logger <- ev %>% mutate(discharge_source = "logger", has_logger = TRUE)
ev_resc   <- ev_logger %>%
  left_join(logger_recovery %>% select(event_id, logger_recovery), by = "event_id") %>%
  mutate(nacl_mass_g = nacl_mass_g * logger_recovery) %>%
  select(-logger_recovery)

run_variant <- function(ev_v, label) {
  message("\n--- ", label, " ---")
  fits <- fit_all_hydraulics(ev_v, btc, mt, n_lhs = N_LHS)
  list(fits = fits,
       table = tsm_hydraulics_table(fits) %>% mutate(series = label, .before = 1))
}
res <- list(
  probe_grab      = run_variant(ev_probe,  "probe_grab"),
  logger          = run_variant(ev_logger, "logger"),
  logger_rescaled = run_variant(ev_resc,   "logger_rescaled")
)

cmp <- map_dfr(res, "table") %>%
  left_join(logger_recovery %>% select(event_id, logger_recovery), by = "event_id") %>%
  select(event_id, series, conservative_source, hydraulic_model, selection_reason,
         v_input_ms, v_ms, D_m2s, alpha_1s, As_over_A, hydraulic_rmse,
         params_at_bound, n_fit, n_pre_zeroed, logger_recovery) %>%
  arrange(event_id, series)

cat("\n=== hydraulics: probe grabs vs logger ===\n")
print(cmp, n = Inf, width = Inf)
write_csv(cmp, file.path(OUT, "ra_logger_vs_probe_hydraulics.csv"))

# ---- plots: both series + the three fitted curves, one figure per event ----
mass_for <- function(eid, label) {
  m <- ev$nacl_mass_g[ev$event_id == eid]
  if (label == "logger_rescaled") m <- m * logger_recovery$logger_recovery[logger_recovery$event_id == eid]
  m
}

for (eid in IDS) {
  lg <- btc %>% filter(event_id == eid, time_since_release_s >= 0) %>%
    transmute(time = time_since_release_s, conc = pmax(nacl_mgL, 0), series = "logger")
  pg <- mt %>% filter(event_id == eid, !is.na(nacl_mgL_grab), time_since_release_s >= 0) %>%
    distinct(time_since_release_s, nacl_mgL_grab) %>%
    transmute(time = time_since_release_s, conc = pmax(nacl_mgL_grab, 0), series = "probe_grab")
  obs <- bind_rows(lg, pg)
  if (nrow(obs) == 0) next
  t_plot <- seq(0, max(obs$time), length.out = 600)
  
  curves <- imap_dfr(res, function(r, label) {
    h <- r$fits[[eid]]
    if (is.null(h)) return(NULL)
    tibble(time = t_plot, series = label,
           C = simulate_tsm(L = h$L, Q = h$Q, A = h$A, D = h$D, alpha = h$alpha, As = h$As,
                            mass = mass_for(eid, label), times = t_plot)$C,
           label = sprintf("%s: %s, v %.4f, D %.3g, As/A %.2f", label, h$hydraulic_model,
                           h$v, h$D, h$As / h$A))
  })
  
  p <- ggplot() +
    geom_point(data = obs, aes(time, conc, shape = series), colour = "grey35", size = 1.4, alpha = 0.8) +
    scale_shape_manual(values = c(logger = 16, probe_grab = 2), name = "observed") +
    geom_line(data = curves, aes(time, C, colour = label), linewidth = 0.9) +
    labs(title = sprintf("%s -- probe grabs vs logger (logger recovery %.0f%%)", eid,
                         100 * logger_recovery$logger_recovery[logger_recovery$event_id == eid]),
         x = "time since release (s)", y = "NaCl (mg/L)", colour = "fit") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical")
  print(p)
  ggsave(file.path(OUT, paste0(eid, "_logger_vs_probe.png")), p, width = 10, height = 6,
         dpi = 150, bg = "white")
}

cat("\nWritten to:", OUT, "\n")