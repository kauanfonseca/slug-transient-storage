# =============================================================================
# PHREEQC speciation for the SRP breakthrough curves
#
# Run ONCE, outside the vignettes:
#
#   source(here::here("R", "phreeqc_speciation.R"))
#
# Depends on a local PHREEQC install and a machine-specific database path, so
# it is kept out of 02_integration.Rmd, which reads only the CSV this writes.
#
# Output: data_derived/phreeqc_activities.csv
#
# EVERY input is the DOWNSTREAM solution: the water passing the sampling
# station is what the model must speciate. pH, temperature and major ions all
# come from the downstream record, and Ca uses ca_mgL_down, never the
# upstream-downstream mean.
#
# Returns phosphate activities only. Calcite saturation is NOT computed: the
# input carries no carbonate species, and the campaign-mean alkalinity
# available is not an event-specific measurement. Evidence for active calcite
# precipitation comes instead from the upstream-downstream response ratios in
# Section 9.4 of 02_integration.Rmd, which use only quantities measured as a
# synoptic pair.
# =============================================================================

library(phreeqc)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(here)

PHREEQC_DB <- "/Users/kauanfonseca/Applications/phreeqc/phreeqc-3.5.0-14000/database/phreeqc.dat"

stopifnot(file.exists(PHREEQC_DB))
phrLoadDatabaseString(paste(readLines(PHREEQC_DB), collapse = "\n"))

M_P  <- 30.973762   # g/mol
M_CA <- 40.078      # g/mol

p_raw <- read_csv(here("data", "nutrients", "nutrient_addition_phosphate_data.csv"),
                  show_col_types = FALSE) %>%
  mutate(srp_molL = measured_SRP_ug_L / M_P / 1e6)

# guard: ph and temp_c must be the downstream values. Field sheets originally
# carried the upstream reading for at least one event; if any event matches the
# synoptic upstream pair exactly, that correction was not applied.
synoptic_path <- here("data", "nutrients", "synoptic_upstream.csv")
if (file.exists(synoptic_path)) {
  chk <- p_raw %>%
    group_by(stream, date) %>%
    summarise(ph = first(na.omit(ph)), temp_c = first(na.omit(temp_c)),
              .groups = "drop") %>%
    left_join(read_csv(synoptic_path, show_col_types = FALSE) %>%
                select(stream, date, ph_up, temp_c_up),
              by = c("stream", "date")) %>%
    filter(abs(ph - ph_up) < 0.02 & abs(temp_c - temp_c_up) < 0.2)
  if (nrow(chk) > 0) {
    print(chk)
    warning("Eventos acima tem ph/temp_c identicos ao montante: verifique se sao valores de jusante.")
  }
}

make_solution <- function(sample_id, srp_molL, ph, temp_c,
                          ca_down, mg, na, k, so4, cl) {
  lines <- c(sprintf("SOLUTION %s", sample_id),
             sprintf("    temp %s", temp_c),
             sprintf("    pH %s", ph),
             "    units mol/L")
  add <- function(label, value) {
    if (!is.na(value)) lines <<- c(lines, sprintf("    %s %g", label, value))
  }
  # ca_down, mg, na, k, so4, cl arrive in mg/L (see the _mgL column names in
  # the nutrient spreadsheet). Converting mg/L to mol/L is a TWO-step scaling
  # -- mg -> g (/1000), then g/L -> mol/L (/molar mass) -- and the previous
  # version of this function only did the second step, so every major ion
  # below was passed to PHREEQC 1000x too concentrated (e.g. ~1 "mol/L" Ca
  # for a stream that is actually ~1 mmol/L, i.e. brine-strength instead of
  # freshwater). That inflates ionic strength and ion-pairing enormously and
  # is the most likely reason a_PO4/a_HPO4 came out implausibly small
  # (~1e-14 to 1e-9, flagged but not explained in 03_btc_review.Rmd Sec.
  # 1.1) and the coprecipitation correction looked negligible. P was NOT
  # affected -- srp_molL above already applies both steps (/M_P/1e6). The
  # /1000 step here matches the identical mg/L -> mol/L conversion already
  # done correctly for calcium in 02_integration.Rmd's coprec-apply chunk
  # (delta_ca_molL = delta_ca_mgL / (1000 * M_CA)).
  add("Ca",   ca_down / 1000 / M_CA)
  add("Mg",   mg      / 1000 / 24.305)
  add("Na",   na      / 1000 / 22.99)
  add("K",    k       / 1000 / 39.098)
  add("S(6)", so4     / 1000 / 96.06)
  add("Cl",   cl      / 1000 / 35.45)
  add("P",    srp_molL)
  lines <- c(lines, "END", "",
             "SELECTED_OUTPUT",
             "    -molalities HPO4-2 PO4-3 H2PO4-",
             "    -activities HPO4-2 PO4-3 H2PO4-",
             "END")
  paste(lines, collapse = "\n")
}

# dummy run: the first call returns NA activities otherwise
phrRunString(paste(
  "SOLUTION 0", "  temp 25", "  pH 7", "END", "",
  "SELECTED_OUTPUT",
  "  -molalities HPO4-2 PO4-3 H2PO4-",
  "  -activities HPO4-2 PO4-3 H2PO4-",
  "END", sep = "\n"))

run_row <- function(row, sample_id) {
  phrRunString(make_solution(sample_id, row$srp_molL, row$ph, row$temp_c,
                             row$ca_mgL_down, row$mg_mgL, row$na_mgL,
                             row$k_mgL, row$so42_mgL, row$cl_mgL))
  out <- phrGetSelectedOutput()$n1
  tibble(a_PO4   = 10^(out$la_PO4.3),
         a_HPO4  = 10^(out$la_HPO4.2),
         a_H2PO4 = 10^(out$la_H2PO4.))
}

ids <- seq_len(nrow(p_raw))

activities <- p_raw %>%
  mutate(phreeqc = map2(split(p_raw, ids), ids, ~ run_row(.x, .y))) %>%
  unnest_wider(phreeqc) %>%
  transmute(stream, date, sample_btc, a_PO4, a_HPO4, a_H2PO4)

# guard: PHREEQC returns NA silently when a solution fails to converge
cat("\nlinhas sem atividade:\n")
activities %>% filter(if_any(c(a_PO4, a_HPO4), is.na)) %>% print(n = Inf)

write_csv(activities, here("data_derived", "phreeqc_activities.csv"))
nrow(activities)