################################################################################
# app_tsm_review.R
#
# Interactive Shiny app for event-by-event review and manual curation of the
# transient-storage model fits (R/tsm_model.R, tsm_calibrate.R,
# tsm_partition.R, run_tsm_uptake.R). Companion to the batch pipeline
# (vignettes/04_tsm_uptake.Rmd): use the batch pass to get a first cut across
# all events, then use this app to fix the event-specific problems a
# one-size-fits-all script cannot catch -- which conservative-tracer series
# to trust (logger vs. hand-held probe), which grab points to exclude as
# field/instrument errors, and eyeballing every fit before accepting it.
#
# Nothing here is silently changed: every accepted decision is appended to
# data_derived/tsm_manual_overrides.csv, recording exactly what was
# overridden and why (chosen conservative-tracer source, excluded
# timestamps, resulting fitted parameters, RMSE, a timestamp). That mirrors
# this project's own QAQC convention (see handoff.md) of never dropping data
# without an explicit, auditable flag/reason.
#
# Workflow:
#   1. Pick an event. Hydraulics tab: check the source (logger vs. probe
#      grab -- defaults to whatever events$discharge_source already says is
#      trustworthy), exclude any bad points, Refit, look at the plot +
#      identifiability scan, and Accept when it looks right.
#   2. Uptake tab: pick the solute (and, for SRP, raw vs. coprec-corrected),
#      exclude any outlier grabs, Refit, look at the plot + identifiability
#      scan + partition summary, and Accept.
#   3. Move to the next event. The "Overrides salvos" tab shows everything
#      accepted so far in this session (and everything ever accepted, since
#      it reads straight from the CSV).
#
# This app does its own from-scratch Stage-1/Stage-2 fits per event (it does
# NOT reuse the batch pipeline's run_tsm_uptake_all() results) so that
# excluding a point or switching source takes effect immediately. A single
# event's fit takes a few seconds to a couple of minutes depending on
# n_lhs/n_cells -- same solver as the batch script, just run one event at a
# time.
#
# Run with:
#   shiny::runApp("R/app_tsm_review.R")
# (install shiny + DT once, if not already installed:
#   install.packages(c("shiny", "DT")))
################################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(here)
})

source(here::here("R", "tsm_model.R"))
source(here::here("R", "tsm_calibrate.R"))
source(here::here("R", "tsm_partition.R"))
source(here::here("R", "plot_tsm_fit.R"))
source(here::here("R", "run_tsm_uptake.R"))

# ---------------------------------------------------------------------------
# Data (loaded once at app start; edits to the underlying CSVs require a
# restart of the app to pick up)
# ---------------------------------------------------------------------------
DATA_DIR       <- here::here("data_derived")
NUTRIENTS_DIR  <- here::here("data", "nutrients")
OVERRIDES_PATH <- file.path(DATA_DIR, "tsm_manual_overrides.csv")

events           <- read_csv(file.path(DATA_DIR, "events.csv"), show_col_types = FALSE)
btc_conservative <- read_csv(file.path(DATA_DIR, "btc_conservative.csv"), show_col_types = FALSE)
master_tsm       <- read_csv(file.path(DATA_DIR, "master_tsm.csv"), show_col_types = FALSE)
nitrogen_raw <- read_csv(file.path(NUTRIENTS_DIR, "nutrient_addition_nitrogen_data.csv"), show_col_types = FALSE) %>%
  mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_NH4Cl_g, molar_mass_NH4Cl_nutrient)
phosphate_raw <- read_csv(file.path(NUTRIENTS_DIR, "nutrient_addition_phosphate_data.csv"), show_col_types = FALSE) %>%
  mutate(date = as.character(date)) %>% distinct(stream, date, added_mass_PO4_g, molar_mass_P_nutrient)

EVENT_IDS <- sort(unique(events$event_id))

OVERRIDE_COLS <- c("timestamp", "event_id", "stage", "solute", "conc_col",
                    "conservative_source", "excluded_times", "n_cells", "n_lhs",
                    "D_m2s", "alpha_1s", "As_m2", "lambda_1s", "lambda_s_1s",
                    "rmse", "n_gap_filled", "gap_fill_dt_s", "notes")

read_overrides <- function() {
  if (file.exists(OVERRIDES_PATH)) {
    read_csv(OVERRIDES_PATH, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else {
    as_tibble(setNames(replicate(length(OVERRIDE_COLS), character(0), simplify = FALSE), OVERRIDE_COLS))
  }
}

# Appends one row, coercing everything to character first -- avoids type
# clashes when binding onto a CSV re-read as all-character columns. This is
# an audit log, not an input table for the pipeline, so character-typed
# numbers (parsed back with as.numeric() by anything that reads it later)
# are an acceptable trade-off for never failing to log a decision.
append_override <- function(row) {
  existing <- read_overrides()
  row_chr <- lapply(row, function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x))
  for (nm in OVERRIDE_COLS) if (is.null(row_chr[[nm]])) row_chr[[nm]] <- NA_character_
  row_tbl <- as_tibble(row_chr[OVERRIDE_COLS])
  write_csv(bind_rows(existing, row_tbl), OVERRIDES_PATH)
}

# Observed conservative-tracer series for one event, from the chosen source.
# "probe_grab" pulls the hand-held-probe NaCl grabs from master_tsm (the same
# series fit_all_hydraulics() falls back to for events whose discharge_source
# is "probe" -- see run_tsm_uptake.R for why some events need this).
get_conservative_series <- function(eid, source) {
  if (identical(source, "probe_grab")) {
    master_tsm %>%
      filter(event_id == eid, !is.na(nacl_mgL_grab), time_since_release_s >= 0) %>%
      distinct(time_since_release_s, nacl_mgL_grab) %>%
      arrange(time_since_release_s) %>%
      rename(nacl_mgL = nacl_mgL_grab)
  } else {
    btc_conservative %>%
      filter(event_id == eid, time_since_release_s >= 0) %>%
      arrange(time_since_release_s)
  }
}

near_bound <- function(val, bounds, tol = 0.05) {
  lv <- log10(val)
  (lv - bounds[1]) < tol * diff(bounds) || (bounds[2] - lv) < tol * diff(bounds)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("TSM review — curadoria manual dos ajustes (hidráulica + uptake)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("event_id", "Event", choices = EVENT_IDS),
      hr(),
      helpText("1) Ajuste a hidráulica (aba Hydraulics) e clique Accept.",
               "2) Só depois ajuste o uptake do soluto desejado (aba Uptake) -- ",
               "precisa da hidráulica aceita para o evento atual.",
               "Cada Accept grava uma linha em data_derived/tsm_manual_overrides.csv."),
      hr(),
      h5("Status desta sessão"),
      verbatimTextOutput("session_status")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Hydraulics",
          br(),
          uiOutput("hyd_borrow_note"),
          fluidRow(
            column(4, radioButtons("cons_source", "Conservative-tracer source",
                                    choices = c("logger" = "logger",
                                                "probe_grab (sonda manual)" = "probe_grab"),
                                    selected = "logger")),
            column(4, numericInput("hyd_n_cells", "n_cells", value = 40, min = 10, max = 120, step = 5)),
            column(4, numericInput("hyd_n_lhs", "n_lhs", value = 250, min = 40, max = 1000, step = 10))
          ),
          uiOutput("hyd_gap_fill_ui"),
          helpText("Selecione linhas na tabela abaixo para EXCLUIR pontos antes de reajustar",
                    "(ex.: erro de leitura, pico de conductividade espúrio)."),
          DTOutput("hyd_table"),
          br(),
          actionButton("hyd_refit", "Refit hydraulics", class = "btn-primary"),
          actionButton("hyd_accept", "Accept hydraulics", class = "btn-success"),
          br(), br(),
          plotOutput("hyd_fit_plot", height = "320px"),
          fluidRow(
            column(6, selectInput("hyd_ident_param", "Identifiability plot: parâmetro",
                                   choices = c("D", "alpha", "As_ratio")))
          ),
          plotOutput("hyd_ident_plot", height = "300px"),
          verbatimTextOutput("hyd_fit_summary")
        ),
        tabPanel("Uptake",
          br(),
          fluidRow(
            column(4, selectInput("solute", "Solute", choices = c("NH4-N", "SRP"))),
            column(4, uiOutput("conc_col_ui")),
            column(4, numericInput("up_n_cells", "n_cells", value = 40, min = 10, max = 120, step = 5))
          ),
          numericInput("up_n_lhs", "n_lhs", value = 250, min = 40, max = 1000, step = 10),
          checkboxInput("up_fill_gap", "Preencher lacuna pré-chegada com background (recomendado)", value = TRUE),
          helpText("Requer hidráulica aceita para este evento (aba Hydraulics).",
                    "Selecione linhas na tabela para EXCLUIR outliers antes de reajustar."),
          DTOutput("up_table"),
          br(),
          actionButton("up_refit", "Refit uptake", class = "btn-primary"),
          actionButton("up_accept", "Accept uptake", class = "btn-success"),
          br(), br(),
          plotOutput("up_fit_plot", height = "320px"),
          fluidRow(
            column(6, selectInput("up_ident_param", "Identifiability plot: parâmetro",
                                   choices = c("lambda", "lambda_s")))
          ),
          plotOutput("up_ident_plot", height = "300px"),
          verbatimTextOutput("up_fit_summary")
        ),
        tabPanel("Overrides salvos",
          br(),
          helpText("Todo Accept fica registrado aqui e em data_derived/tsm_manual_overrides.csv."),
          DTOutput("overrides_table")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    hyd_fit = NULL,          # current (not-yet-accepted) hydraulics fit
    accepted_hyd = list(),   # event_id -> accepted hydraulics list
    up_fit = NULL,           # current (not-yet-accepted) uptake fit
    accepted_up = list()     # "event_id|solute|conc_col" -> accepted uptake list
  )

  # default the source radio to whatever events$discharge_source says is
  # trustworthy for this event, and clear any in-progress (unaccepted) fit
  observeEvent(input$event_id, {
    e <- events %>% filter(event_id == input$event_id)
    default_src <- if (nrow(e) > 0 && isTRUE(e$discharge_source == "probe")) "probe_grab" else "logger"
    updateRadioButtons(session, "cons_source", selected = default_src)
    rv$hyd_fit <- NULL
    rv$up_fit  <- NULL
  }, ignoreInit = FALSE)

  output$hyd_borrow_note <- renderUI({
    eid <- input$event_id
    if (eid %in% names(HYDRAULICS_BORROWED_FROM)) {
      src <- HYDRAULICS_BORROWED_FROM[[eid]]
      div(class = "alert alert-warning",
          sprintf(paste0("Este evento não tem curva conservativa própria de boa qualidade. ",
                          "No pipeline em lote, a hidráulica é emprestada de %s (mesmo trecho). ",
                          "Você ainda pode tentar um ajuste manual aqui (ex. usando probe_grab, se houver), ",
                          "mas considere se faz mais sentido aceitar diretamente os valores do evento fonte."),
                  src))
    }
  })

  # ---- Hydraulics: observed series ------------------------------------------
  hyd_series <- reactive({
    req(input$event_id)
    get_conservative_series(input$event_id, input$cons_source)
  })

  output$hyd_table <- renderDT({
    datatable(hyd_series(), selection = "multiple", rownames = FALSE,
              options = list(pageLength = 10, order = list(list(0, "asc"))))
  })

  # the pre-arrival gap fill only makes sense for the hand-held-probe grabs
  # (the logger already records continuously, so has no such gap -- see
  # fill_pre_arrival_gap(), tsm_calibrate.R) -- hide the toggle otherwise
  output$hyd_gap_fill_ui <- renderUI({
    if (identical(input$cons_source, "probe_grab")) {
      checkboxInput("hyd_fill_gap", "Preencher lacuna pré-chegada com background (recomendado)", value = TRUE)
    }
  })

  observeEvent(input$hyd_refit, {
    df_full <- hyd_series()
    sel <- input$hyd_table_rows_selected
    excluded_times <- if (!is.null(sel) && length(sel) > 0) df_full$time_since_release_s[sel] else numeric(0)
    df <- if (length(excluded_times) > 0) df_full[-sel, , drop = FALSE] else df_full
    validate(need(nrow(df) >= 8, "Poucos pontos após exclusão (mínimo 8) para ajustar a hidráulica."))

    e <- events %>% filter(event_id == input$event_id)
    L <- e$reach_length_m
    Q <- e$discharge_Ls / 1000
    A <- Q / e$water_velocity_ms

    do_fill <- identical(input$cons_source, "probe_grab") && isTRUE(input$hyd_fill_gap)
    obs_conc <- pmax(df$nacl_mgL, 0)
    gf <- if (do_fill) {
      fill_pre_arrival_gap(df$time_since_release_s, obs_conc)
    } else {
      list(time = df$time_since_release_s, value = obs_conc,
           kind = rep("observed", nrow(df)), n_added = 0L, dt_used = NA_real_)
    }

    fit <- withProgress(message = "Ajustando hidráulica...", value = 0.3, {
      out <- tryCatch(
        fit_hydraulics(gf$time, gf$value, L, Q, A,
                        e$nacl_mass_g, n_lhs = input$hyd_n_lhs, n_cells = input$hyd_n_cells),
        error = function(err) NULL
      )
      incProgress(0.7)
      out
    })
    validate(need(!is.null(fit),
                  "Falha no ajuste (erro no solver) -- tente outro n_cells ou revise a exclusão de pontos."))

    rv$hyd_fit <- list(event_id = input$event_id, L = L, Q = Q, A = A, v = Q / A,
                        width = e$reach_mean_width_m, D = fit$par[["D"]],
                        alpha = fit$par[["alpha"]], As = fit$par[["As"]], rmse = fit$rmse,
                        lhs = fit$lhs, conservative_source = input$cons_source,
                        excluded_times = excluded_times, obs = df, mass_g = e$nacl_mass_g,
                        plot_time = gf$time, plot_conc = gf$value, plot_kind = gf$kind,
                        n_gap_filled = gf$n_added, gap_fill_dt_s = gf$dt_used)
  })

  output$hyd_fit_plot <- renderPlot({
    req(rv$hyd_fit)
    h <- rv$hyd_fit
    plot_btc_fit(h$plot_time, h$plot_conc,
                 h$L, h$Q, h$A, h$D, h$alpha, h$As, mass = h$mass_g,
                 n_cells = input$hyd_n_cells, unit = "mg/L", obs_kind = h$plot_kind,
                 title = sprintf("%s -- NaCl fit (RMSE=%.3f, source=%s)", h$event_id, h$rmse, h$conservative_source))
  })

  output$hyd_ident_plot <- renderPlot({
    req(rv$hyd_fit, rv$hyd_fit$lhs)
    plot_identifiability(rv$hyd_fit$lhs, input$hyd_ident_param, xlab = input$hyd_ident_param)
  })

  output$hyd_fit_summary <- renderPrint({
    req(rv$hyd_fit)
    h <- rv$hyd_fit
    cat(sprintf("D = %.5g m2/s | alpha = %.5g 1/s | As = %.4g m2 (As/A = %.3f) | RMSE = %.4f\n",
                h$D, h$alpha, h$As, h$As / h$A, h$rmse))
    cat(sprintf("source = %s | pontos excluídos = %d\n", h$conservative_source, length(h$excluded_times)))
    if (h$n_gap_filled > 0) {
      cat(sprintf("lacuna pré-chegada preenchida: %d pontos sintéticos (background=0) a cada %.0f s, de t=0 até a 1ª amostra real\n",
                  h$n_gap_filled, h$gap_fill_dt_s))
    }
    if (input$event_id %in% names(rv$accepted_hyd)) {
      cat("\n[ACEITA para esta sessão -- a aba Uptake vai usar este ajuste]\n")
    } else {
      cat("\n[AINDA NÃO ACEITA -- clique 'Accept hydraulics' antes de ajustar o uptake]\n")
    }
  })

  observeEvent(input$hyd_accept, {
    req(rv$hyd_fit)
    h <- rv$hyd_fit
    rv$accepted_hyd[[input$event_id]] <- h
    append_override(list(
      timestamp = as.character(Sys.time()), event_id = input$event_id, stage = "hydraulics",
      solute = NA, conc_col = NA, conservative_source = h$conservative_source,
      excluded_times = paste(h$excluded_times, collapse = ";"),
      n_cells = input$hyd_n_cells, n_lhs = input$hyd_n_lhs,
      n_gap_filled = h$n_gap_filled, gap_fill_dt_s = h$gap_fill_dt_s,
      D_m2s = h$D, alpha_1s = h$alpha, As_m2 = h$As, lambda_1s = NA, lambda_s_1s = NA,
      rmse = h$rmse, notes = "accepted via app_tsm_review"
    ))
    showNotification(sprintf("Hidráulica de %s aceita e registrada.", input$event_id), type = "message")
  })

  # ---- Uptake: observed series -----------------------------------------------
  output$conc_col_ui <- renderUI({
    if (identical(input$solute, "SRP")) {
      selectInput("conc_col", "Concentration column",
                  choices = c("raw (conc_corr_ugL)" = "conc_corr_ugL",
                              "coprec-corrected, biotic-only (conc_corr_coprec_ugL)" = "conc_corr_coprec_ugL"))
    } else {
      selectInput("conc_col", "Concentration column",
                  choices = c("raw (conc_corr_ugL)" = "conc_corr_ugL"))
    }
  })

  up_series <- reactive({
    req(input$event_id, input$solute, input$conc_col)
    master_tsm %>%
      filter(event_id == input$event_id, solute == input$solute, !is.na(.data[[input$conc_col]])) %>%
      arrange(time_since_release_s) %>%
      select(time_since_release_s, conc_ugL = all_of(input$conc_col), background_ugL)
  })

  output$up_table <- renderDT({
    datatable(up_series(), selection = "multiple", rownames = FALSE,
              options = list(pageLength = 10, order = list(list(0, "asc"))))
  })

  observeEvent(input$up_refit, {
    hyd <- rv$accepted_hyd[[input$event_id]]
    validate(need(!is.null(hyd), "Aceite a hidráulica deste evento primeiro (aba Hydraulics)."))

    df_full <- up_series()
    sel <- input$up_table_rows_selected
    excluded_times <- if (!is.null(sel) && length(sel) > 0) df_full$time_since_release_s[sel] else numeric(0)
    df <- if (length(excluded_times) > 0) df_full[-sel, , drop = FALSE] else df_full
    validate(need(nrow(df) >= 6, "Poucos pontos após exclusão (mínimo 6) para ajustar o uptake."))

    e <- events %>% filter(event_id == input$event_id)
    mass_mg <- lookup_injected_mass_mg(e, nitrogen_raw, phosphate_raw, input$solute)
    validate(need(is.finite(mass_mg), "Não há registro de massa injetada para este evento/soluto."))

    do_fill <- isTRUE(input$up_fill_gap)
    obs_conc <- pmax(df$conc_ugL, 0)
    gf <- if (do_fill) {
      fill_pre_arrival_gap(df$time_since_release_s, obs_conc)
    } else {
      list(time = df$time_since_release_s, value = obs_conc,
           kind = rep("observed", nrow(df)), n_added = 0L, dt_used = NA_real_)
    }

    hydraulics <- c(D = hyd$D, alpha = hyd$alpha, As = hyd$As)
    fitU <- withProgress(message = "Ajustando uptake...", value = 0.3, {
      out <- tryCatch(
        fit_uptake(gf$time, gf$value, hyd$L, hyd$Q, hyd$A, hydraulics, mass_mg,
                   n_lhs = input$up_n_lhs, n_cells = input$up_n_cells),
        error = function(err) NULL
      )
      incProgress(0.6)
      out
    })
    validate(need(!is.null(fitU), "Falha no ajuste do uptake."))

    part <- tryCatch(
      partition_uptake(hyd$L, hyd$Q, hyd$A, hyd$D, hyd$alpha, hyd$As,
                        fitU$par[["lambda"]], fitU$par[["lambda_s"]], mass_mg, n_cells = input$up_n_cells),
      error = function(err) NULL
    )
    depth_main <- hyd$A / hyd$width
    depth_storage <- hyd$As / hyd$width
    Camb <- suppressWarnings(mean(df$background_ugL, na.rm = TRUE))
    met <- uptake_metrics(hyd$v, depth_main, depth_storage, fitU$par[["lambda"]],
                           fitU$par[["lambda_s"]], hyd$alpha, hyd$A, hyd$As, Camb)

    lambda_at_bound   <- near_bound(fitU$par[["lambda"]], c(-7, -1))
    lambda_s_at_bound <- near_bound(fitU$par[["lambda_s"]], c(-7, -1))
    uptake_significant <- is.finite(part$pct_total_uptake) && part$pct_total_uptake > 2

    rv$up_fit <- list(event_id = input$event_id, solute = input$solute, conc_col = input$conc_col,
                       lambda = fitU$par[["lambda"]], lambda_s = fitU$par[["lambda_s"]],
                       rmse = fitU$rmse, lhs = fitU$lhs, mass_mg = mass_mg,
                       excluded_times = excluded_times, obs = df, hyd = hyd,
                       plot_time = gf$time, plot_conc = gf$value, plot_kind = gf$kind,
                       n_gap_filled = gf$n_added, gap_fill_dt_s = gf$dt_used,
                       part = part, met = met, lambda_at_bound = lambda_at_bound,
                       lambda_s_at_bound = lambda_s_at_bound, uptake_significant = uptake_significant)
  })

  output$up_fit_plot <- renderPlot({
    req(rv$up_fit)
    u <- rv$up_fit
    plot_btc_fit(u$plot_time, u$plot_conc,
                 u$hyd$L, u$hyd$Q, u$hyd$A, u$hyd$D, u$hyd$alpha, u$hyd$As,
                 lambda = u$lambda, lambda_s = u$lambda_s, mass = u$mass_mg,
                 n_cells = input$up_n_cells, unit = "ug/L", obs_kind = u$plot_kind,
                 title = sprintf("%s -- %s (%s) RMSE=%.3f", u$event_id, u$solute, u$conc_col, u$rmse))
  })

  output$up_ident_plot <- renderPlot({
    req(rv$up_fit, rv$up_fit$lhs)
    plot_identifiability(rv$up_fit$lhs, input$up_ident_param, xlab = input$up_ident_param)
  })

  output$up_fit_summary <- renderPrint({
    req(rv$up_fit)
    u <- rv$up_fit
    cat(sprintf("lambda (canal) = %.4g 1/s %s\n", u$lambda,
                if (u$lambda_at_bound) "[NO LIMITE DA BUSCA -- checar identificabilidade]" else ""))
    cat(sprintf("lambda_s (zona de transição) = %.4g 1/s %s\n", u$lambda_s,
                if (u$lambda_s_at_bound) "[NO LIMITE DA BUSCA -- checar identificabilidade]" else ""))
    cat(sprintf("RMSE = %.4f | pontos excluídos = %d\n", u$rmse, length(u$excluded_times)))
    if (u$n_gap_filled > 0) {
      cat(sprintf("lacuna pré-chegada preenchida: %d pontos sintéticos (background=0) a cada %.0f s\n",
                  u$n_gap_filled, u$gap_fill_dt_s))
    }
    cat("\n")
    cat(sprintf("Uptake total do trecho: %.1f%%\n", u$part$pct_total_uptake))
    if (u$uptake_significant) {
      cat(sprintf("  -> canal principal: %.1f%% | zona de transição: %.1f%%\n",
                  u$part$pct_mainchannel, u$part$pct_storagezone))
    } else {
      cat("  -> uptake total muito baixo (<2%) -- partição não reportada (não significativa)\n")
    }
    if (!u$lambda_at_bound && u$uptake_significant) {
      cat(sprintf("\nSw (canal) = %.1f m | Sw (total, agregado) = %.1f m\n",
                  u$met$Sw_mainchannel_m, u$met$Sw_total_m))
    }
    cat(sprintf("U total = %.4g mg/m2/h\n", u$met$U_total * 3600))
    if (input$event_id %in% names(rv$accepted_hyd)) {
      hsrc <- rv$accepted_hyd[[input$event_id]]$conservative_source
      cat(sprintf("\n(hidráulica aceita usada: source=%s)\n", hsrc))
    }
  })

  observeEvent(input$up_accept, {
    req(rv$up_fit)
    u <- rv$up_fit
    key <- paste(u$event_id, u$solute, u$conc_col, sep = "|")
    rv$accepted_up[[key]] <- u
    append_override(list(
      timestamp = as.character(Sys.time()), event_id = u$event_id, stage = "uptake",
      solute = u$solute, conc_col = u$conc_col, conservative_source = u$hyd$conservative_source,
      excluded_times = paste(u$excluded_times, collapse = ";"),
      n_cells = input$up_n_cells, n_lhs = input$up_n_lhs,
      n_gap_filled = u$n_gap_filled, gap_fill_dt_s = u$gap_fill_dt_s,
      D_m2s = u$hyd$D, alpha_1s = u$hyd$alpha, As_m2 = u$hyd$As,
      lambda_1s = u$lambda, lambda_s_1s = u$lambda_s,
      rmse = u$rmse, notes = "accepted via app_tsm_review"
    ))
    showNotification(sprintf("Uptake de %s | %s | %s aceito e registrado.", u$event_id, u$solute, u$conc_col),
                      type = "message")
  })

  # ---- Overrides table / session status --------------------------------------
  output$overrides_table <- renderDT({
    input$hyd_accept; input$up_accept  # invalidate + redraw whenever something is accepted
    datatable(read_overrides(), rownames = FALSE, options = list(pageLength = 10))
  })

  output$session_status <- renderPrint({
    cat("Hidráulicas aceitas:", length(rv$accepted_hyd), "evento(s):\n")
    if (length(rv$accepted_hyd) > 0) cat(" ", paste(names(rv$accepted_hyd), collapse = "\n  "), "\n")
    cat("\nUptakes aceitos:", length(rv$accepted_up), "combinação(ões):\n")
    if (length(rv$accepted_up) > 0) cat(" ", paste(names(rv$accepted_up), collapse = "\n  "), "\n")
  })
}

shinyApp(ui, server)
