# =============================================================================
# QAQC — campaign discharge records
#
# Reads every handheld slug record of the campaign, repairs the known
# transcription problems, computes discharge by integral dilution gauging, and
# assigns a discharge to each calcium sampling so the calcium transfer of
# Section 9 can be normalised to an areal rate.
#
# Sources
#   data/field/discharge_spc_records.csv     20 event-dates, all sites
#   data/field/discharge_RA_20230804.csv     RA first attempt, separate file
#   data/field/calcium_reach_change.csv      the calcium pairs to be matched
#
# Outputs
#   data_derived/discharge_campaign_qaqc.csv
#   data_derived/discharge_for_calcium.csv
#
# Standalone: run in the console, does not touch the vignettes.
#
# READING NOTE. These files use a POINT as decimal separator. Passing
# locale(decimal_mark = ",") makes readr treat the point as a thousands
# separator and concatenate the digits — 1.3833 becomes 13833 and 423.55
# becomes 42355. That looks like corrupt data but is purely a reading
# artefact. Do not set a locale.
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(pracma)
library(here)

NACL_SLOPE <- 0.461159   # mg/L NaCl per uS/cm

# =============================================================================
# 0. Known corrections, each with its reason
# =============================================================================

# Dates were converted automatically when discharge_spc_records.csv was built
# and some months shifted from 08 to 09. Corrected against the field sheets.
DATE_FIX <- tribble(
  ~Site, ~Date_errada,           ~Date_certa,
  "CC",  as.Date("2023-02-02"),  as.Date("2023-08-02"),
  "CB",  as.Date("2023-09-30"),  as.Date("2023-08-31"),
  "RA",  as.Date("2023-09-30"),  as.Date("2023-08-30")
)

# RS 2023-08-29: the sheet carries 1258 g, giving 153 L/s against 1145 L/s for
# the same reach six weeks later from a near-identical curve (integrals 8214
# and 9126). Equal integrals mean the discrepancy is entirely in the mass.
MASS_OVERRIDE <- tribble(
  ~Site, ~Date,                  ~nacl_g_corrigido, ~motivo_massa,
  "RS",  as.Date("2023-08-29"),  10291,             "ficha de campo; planilha trazia 1258 g"
)

# RC 2023-08-10 was logged by two instruments. The HOBO series is raw
# conductivity with no temperature compensation, so it is not comparable with
# the specific conductance of every other event, and mixing the two produced a
# curve ending far below background. Only the Hanna series is used.
INSTRUMENT_DROP <- tribble(
  ~Site, ~Date,                  ~Probe_meter,      ~motivo_instrumento,
  "RC",  as.Date("2023-08-10"),  "HOBO_CondLogger", "condutividade bruta, sem compensacao de temperatura"
)

# CC 2023-08-02: two slugs the same day. The first came out with an anomalous
# curve, so the reach was lengthened from 126 m to 286 m and the addition
# repeated. Reach length separates them in the key. The first attempt is kept
# as a comparison anchor but excluded from the analysis; the calcium pair
# corresponds to the 286 m reach.
ATTEMPT_EXCLUDE <- tribble(
  ~Site, ~Date,                  ~reach_m, ~motivo_ancora,
  "CC",  as.Date("2023-08-02"),  126,      "primeira tentativa; curva anomala, trecho depois alongado"
)

# RA 2023-08-04: first attempt at this site, on a longer reach (206 m) whose
# downstream station later became the upstream station of the sampling reach,
# and containing a pool over 2 m deep. Two problems make the curve unreliable:
# the recession drops straight to zero instead of decaying asymptotically, so
# the record was stopped rather than the tracer having finished passing; and
# the field notebooks disagree on the mass, 6 kg in one and 12 kg in the other.
# The integral with 6 kg gives 241 L/s, which agrees with the mean of the four
# other RA gaugings (240 L/s) — two independent routes to the same value. The
# 12 kg reading would give 481 L/s, supported by nothing else. 240 is adopted
# and the alternative is carried in the sensitivity table.
RA_20230804_Q      <- 240
RA_20230804_Q_ALT  <- 481

# Calcium was sampled on the day of gauging except where noted.
CA_DATE_MAP <- tribble(
  ~Site, ~date_ca,               ~date_q,                ~motivo_data,
  "RS",  as.Date("2023-08-28"),  as.Date("2023-08-29"),  "vazao medida no dia seguinte"
)

# =============================================================================
# 1. Read and repair the main file
# =============================================================================

d_raw <- read_delim(here("data", "field", "discharge_spc_records.csv"),
                    delim = ";", show_col_types = FALSE) %>%
  mutate(Date = dmy(Date))

stopifnot(!any(is.na(d_raw$Date)))

d <- d_raw %>%
  left_join(DATE_FIX, by = c("Site", "Date" = "Date_errada")) %>%
  mutate(Date = coalesce(Date_certa, Date)) %>%
  select(-Date_certa) %>%
  anti_join(INSTRUMENT_DROP %>% select(Site, Date, Probe_meter),
            by = c("Site", "Date", "Probe_meter")) %>%
  mutate(key = paste(Site, format(Date, "%Y%m%d"), Reach_Lenght_m))

cat("\n=== estrutura ===\n")
cat("linhas lidas:", nrow(d_raw), " apos reparos:", nrow(d),
    " eventos:", n_distinct(d$key), "\n")

cat("\n--- amplitude das variaveis ---\n")
d %>%
  summarise(across(c(SpC_uScm, Backg_SpC_uScm, NaCl_conc_mgL,
                     Time_s, NaCl_added_g, Reach_Lenght_m),
                   list(min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything()) %>% print(n = Inf)

# SpC outside 1-2000 uS/cm means the decimal was lost on reading. SA is the
# genuine exception: its background really is ~5 uS/cm.
cat("\n--- SpC fora de escala (deve ser vazio) ---\n")
d %>% filter(SpC_uScm > 2000 | (SpC_uScm < 1 & Site != "SA")) %>%
  count(Site, Date) %>% print(n = Inf)

cat("\n--- eventos com mais de um instrumento (deve ser vazio) ---\n")
d %>% distinct(Site, Date, Reach_Lenght_m, Probe_meter) %>%
  count(Site, Date, Reach_Lenght_m, name = "n_instr") %>%
  filter(n_instr > 1) %>% print(n = Inf)

# Mass must be constant within a slug. RA 2023-09-06 legitimately carries two
# values because the P and N slugs were released separately that day.
cat("\n--- massa inconsistente dentro do evento ---\n")
d %>%
  group_by(Site, Date, Reach_Lenght_m, Slug_type) %>%
  summarise(n_massas = n_distinct(NaCl_added_g),
            massas = paste(unique(NaCl_added_g), collapse = " / "),
            .groups = "drop") %>%
  filter(n_massas > 1) %>% print(n = Inf)

# =============================================================================
# 2. Per-event QAQC and discharge
# =============================================================================

qaqc <- d %>%
  group_by(Site, Date, key, Slug_type, Reach_Lenght_m) %>%
  summarise(
    n        = n(),
    dt_med_s = median(diff(Time_s)),
    dt_max_s = max(diff(Time_s)),
    t_max_s  = max(Time_s),

    bg         = first(Backg_SpC_uScm),
    spc_max    = max(SpC_uScm),
    peak_ratio = spc_max / bg,
    # fraction of the peak still present at the last reading; ~0 = full return
    retorno    = (last(SpC_uScm) - bg) / (spc_max - bg),

    conc_neg   = sum(NaCl_conc_mgL < 0, na.rm = TRUE),
    nacl_g_raw = first(NaCl_added_g),
    n_g        = first(N_added_g),
    p_g        = first(P_added_g),
    probe      = first(Probe_meter),
    travertino = first(Travertine_formation),

    integral_Ct = trapz(Time_s, pmax(NaCl_conc_mgL, 0)),
    .groups = "drop"
  ) %>%
  rename(reach_m = Reach_Lenght_m) %>%
  left_join(MASS_OVERRIDE,   by = c("Site", "Date")) %>%
  left_join(ATTEMPT_EXCLUDE, by = c("Site", "Date", "reach_m")) %>%
  mutate(
    nacl_g = coalesce(nacl_g_corrigido, nacl_g_raw),
    Q_Ls   = nacl_g * 1000 / integral_Ct,

    tem_nutrientes = n_g > 0 | p_g > 0,

    flag_massa_corrigida = !is.na(nacl_g_corrigido),
    flag_ancora          = !is.na(motivo_ancora),
    # abs(): a curve ending below background is as wrong as one ending above
    flag_truncado        = abs(retorno) > 0.05,
    flag_mistura         = peak_ratio > 2,
    flag_sinal_fraco     = peak_ratio < 1.05
  )

# ---- RA 2023-08-04, from its own file --------------------------------------

ra_path <- here("data", "field", "discharge_RA_20230804.csv")

if (file.exists(ra_path)) {
  ra <- read_csv(ra_path, show_col_types = FALSE)
  # NOTE: the column named Backg_SpC_uScm in this file is in fact the
  # background-CORRECTED conductance (SpC minus background), not the
  # background itself. The background is the plateau before arrival.
  bg_ra <- ra$SpC_uScm[1]

  ra_row <- ra %>%
    summarise(
      Site = first(Site), Date = first(Date),
      Slug_type = first(Slug_type), reach_m = first(Reach_Lenght_m),
      n = n(), dt_med_s = median(diff(Times_s)),
      dt_max_s = max(diff(Times_s)), t_max_s = max(Times_s),
      bg = bg_ra, spc_max = max(SpC_uScm),
      peak_ratio = spc_max / bg,
      retorno = (last(SpC_uScm) - bg) / (spc_max - bg),
      conc_neg = sum(NaCl_conc_mgL < 0, na.rm = TRUE),
      nacl_g_raw = first(NaCl_added_g),
      n_g = 0, p_g = 0,
      probe = first(Probe_meter), travertino = NA_character_,
      integral_Ct = trapz(Times_s, pmax(NaCl_conc_mgL, 0))
    ) %>%
    mutate(
      key = paste(Site, format(Date, "%Y%m%d"), reach_m),
      nacl_g_corrigido = NA_real_, motivo_massa = NA_character_,
      motivo_ancora = NA_character_,
      nacl_g = nacl_g_raw,
      Q_Ls = RA_20230804_Q,          # adopted value, see header
      tem_nutrientes = FALSE,
      flag_massa_corrigida = FALSE, flag_ancora = FALSE,
      flag_truncado = TRUE,          # recession stopped, does not decay to zero
      flag_mistura = peak_ratio > 2,
      flag_sinal_fraco = peak_ratio < 1.05
    )

  cat("\n--- RA 2023-08-04, arquivo proprio ---\n")
  ra_row %>% select(key, n, t_max_s, bg, spc_max, integral_Ct,
                    nacl_g_raw, Q_Ls, flag_truncado) %>% print()
  cat("integral da curva com 6 kg daria",
      round(6000 * 1000 / ra_row$integral_Ct, 1), "L/s;",
      "com 12 kg,", round(12000 * 1000 / ra_row$integral_Ct, 1), "L/s\n")

  qaqc <- bind_rows(qaqc, ra_row)
}

qaqc <- qaqc %>% arrange(Site, Date, reach_m)

cat("\n=== QAQC por evento ===\n")
qaqc %>%
  select(key, Slug_type, probe, tem_nutrientes, n, dt_med_s, dt_max_s,
         t_max_s, bg, spc_max, peak_ratio, retorno, conc_neg) %>%
  print(n = Inf)

cat("\n=== vazao e flags ===\n")
qaqc %>%
  select(key, nacl_g_raw, nacl_g, reach_m, integral_Ct, Q_Ls,
         flag_massa_corrigida, flag_ancora, flag_truncado, flag_mistura,
         flag_sinal_fraco, travertino) %>%
  print(n = Inf)

# =============================================================================
# 3. Sampling gaps
# =============================================================================
# A long unsampled stretch is bridged by a straight line in the trapezoidal
# integral. Where one interval carries much of the area, the estimate depends
# more on that interpolation than on the data.

gaps <- d %>%
  filter(!is.na(NaCl_conc_mgL)) %>%
  arrange(key, Time_s) %>%
  group_by(key) %>%
  mutate(dt = lead(Time_s) - Time_s,
         area = (NaCl_conc_mgL + lead(NaCl_conc_mgL)) / 2 * dt) %>%
  summarise(max_gap_s = max(dt, na.rm = TRUE),
            frac_area_max_gap = max(area, na.rm = TRUE) / sum(area, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(flag_recessao_esparsa = frac_area_max_gap > 0.3)

cat("\n=== lacunas de amostragem ===\n")
gaps %>% arrange(desc(frac_area_max_gap)) %>% print(n = Inf)

qaqc <- qaqc %>% left_join(gaps, by = "key")

# =============================================================================
# 4. Consistency between dates at the same site
# =============================================================================
# Discharge falls through the dry season, so differences between dates are
# expected. An extreme ratio points to a mass or integral error: this is how
# the RS problem was found, where the ratio was 8 before the correction.

cat("\n=== variacao entre datas no mesmo ponto ===\n")
qaqc %>%
  filter(!flag_ancora) %>%
  group_by(Site) %>% filter(n() > 1) %>%
  summarise(n_datas = n(),
            datas = paste(format(Date, "%d/%m"), collapse = ", "),
            Q_min = min(Q_Ls), Q_max = max(Q_Ls),
            razao = Q_max / Q_min,
            cv = sd(Q_Ls) / mean(Q_Ls),
            .groups = "drop") %>%
  arrange(desc(razao)) %>% print(n = Inf)

cat("\n=== CC 2023-08-02: tentativa descartada vs adotada ===\n")
qaqc %>% filter(Site == "CC", Date == as.Date("2023-08-02")) %>%
  select(key, reach_m, n, integral_Ct, nacl_g, Q_Ls, flag_ancora) %>% print()

# =============================================================================
# 5. Independent check against the pipeline
# =============================================================================
# The nutrient events exist in both sources and are computed from the same
# grabs by independent code paths, so they must agree. The key carries the slug
# label, otherwise RA 2023-09-06 (P and N) joins as a cartesian product.

if (exists("events")) {
  cat("\n=== conferencia contra o pipeline ===\n")
  qaqc %>%
    filter(tem_nutrientes) %>%
    mutate(event_id = case_when(
      Slug_type == "nutrient_spiral_N" ~ paste0(Site, "_", format(Date, "%Y%m%d"), "_N"),
      Slug_type == "nutrient_spiral_P" ~ paste0(Site, "_", format(Date, "%Y%m%d"), "_P"),
      Site == "RA" & Date == as.Date("2023-10-05") ~ "RA_20231005_downstream",
      TRUE ~ paste0(Site, "_", format(Date, "%Y%m%d"), "_single")
    )) %>%
    inner_join(events %>% select(event_id, Q_pipeline = discharge_Ls,
                                 Q_probe = discharge_Ls_probe,
                                 fonte = discharge_source,
                                 reach_pipeline_m = reach_length_m),
               by = "event_id") %>%
    transmute(event_id, Q_Ls, Q_probe, Q_pipeline, fonte,
              razao_probe = Q_Ls / Q_probe,
              reach_m, reach_pipeline_m,
              reach_difere = reach_m != reach_pipeline_m) %>%
    arrange(desc(abs(razao_probe - 1))) %>% print(n = Inf)
} else {
  cat("\n(objeto `events` ausente: rode 01 e 02 para a conferencia)\n")
}

# =============================================================================
# 6. Discharge for each calcium sampling
# =============================================================================
# Driven by the calcium table, so every calcium date appears whether or not a
# gauging matches. Same day where it exists, the mapped exception otherwise,
# and NA where neither applies rather than a silent substitution.

ca_path <- here("data", "field", "calcium_reach_change.csv")

if (file.exists(ca_path)) {

  ca_dates <- read_csv(ca_path, show_col_types = FALSE) %>%
    mutate(date_ca = ymd(as.character(date))) %>%
    select(Site = stream, date_ca, reach_ca_m = reach_length_m)

  discharge_for_ca <- ca_dates %>%
    left_join(CA_DATE_MAP, by = c("Site", "date_ca")) %>%
    mutate(date_q = coalesce(date_q, date_ca),
           motivo_data = coalesce(motivo_data, "mesmo dia")) %>%
    left_join(
      qaqc %>% filter(!flag_ancora) %>%
        select(Site, date_q = Date, Q_Ls, reach_q_m = reach_m,
               flag_massa_corrigida, flag_truncado, flag_mistura,
               flag_sinal_fraco, flag_recessao_esparsa),
      by = c("Site", "date_q")
    ) %>%
    mutate(dist_dias = as.numeric(abs(date_ca - date_q)),
           reach_difere = !is.na(reach_q_m) & reach_ca_m != reach_q_m,
           flag_sem_vazao = is.na(Q_Ls)) %>%
    arrange(Site, date_ca)

  cat("\n=== vazao atribuida a cada coleta de calcio ===\n")
  discharge_for_ca %>%
    select(Site, date_ca, date_q, dist_dias, Q_Ls, reach_ca_m, reach_q_m,
           reach_difere, motivo_data, flag_sem_vazao, flag_truncado,
           flag_sinal_fraco) %>%
    print(n = Inf)

  cat("\n--- coletas de calcio sem vazao ---\n")
  discharge_for_ca %>% filter(flag_sem_vazao) %>%
    select(Site, date_ca) %>% print(n = Inf)

  # Sensitivity: where the assigned discharge is uncertain, the areal rate
  # (Section 9.7 of 02_integration.Rmd) is linear in Q, so the ratio between
  # candidate discharges IS the spread of the rate. Carried into the exported
  # file — not just printed here — so 02_integration.Rmd can report the
  # spread instead of hiding it inside a single choice.
  discharge_for_ca <- discharge_for_ca %>%
    left_join(qaqc %>% filter(!flag_ancora) %>%
                group_by(Site) %>%
                summarise(Q_media_ponto = mean(Q_Ls),
                          Q_min_ponto = min(Q_Ls),
                          Q_max_ponto = max(Q_Ls),
                          n_datas_ponto = n(), .groups = "drop"),
              by = "Site") %>%
    mutate(Q_alt = if_else(Site == "RA" & date_ca == as.Date("2023-08-04"),
                           RA_20230804_Q_ALT, NA_real_),
           razao_media = Q_Ls / Q_media_ponto)

  cat("\n=== sensibilidade da vazao atribuida ===\n")
  discharge_for_ca %>%
    select(Site, date_ca, Q_Ls, Q_alt, Q_media_ponto, Q_min_ponto,
           Q_max_ponto, n_datas_ponto, razao_media) %>%
    print(n = Inf)

  write_csv(discharge_for_ca, here("data_derived", "discharge_for_calcium.csv"))
} else {
  cat("\n(calcium_reach_change.csv ausente)\n")
}

# =============================================================================
# 7. Curves
# =============================================================================

d %>%
  ggplot(aes(Time_s, NaCl_conc_mgL)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.5) +
  facet_wrap(~ key, scales = "free", ncol = 4) +
  labs(x = "tempo desde a liberacao (s)", y = "NaCl (mg/L)") +
  theme_minimal(base_size = 9)

# =============================================================================
# 8. Export
# =============================================================================

write_csv(qaqc, here("data_derived", "discharge_campaign_qaqc.csv"))
cat("\nescrito: discharge_campaign_qaqc.csv, discharge_for_calcium.csv\n")
