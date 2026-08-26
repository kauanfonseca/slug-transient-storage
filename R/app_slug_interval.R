library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(here)

hobo_with_meta <- readRDS(here("data", "hobo_with_meta.rds"))
arquivos_disponiveis <- sort(unique(hobo_with_meta$source_file))

pasta_saida <- here("data", "slug_trimmed")
if (!dir.exists(pasta_saida)) dir.create(pasta_saida, recursive = TRUE)

caminho_log <- here("data", "slug_trimmed_log.csv")

log_inicial <- if (file.exists(caminho_log)) {
  read_csv(caminho_log, show_col_types = FALSE)
} else {
  tibble::tibble(source_file = character(), stream = character(),
                 date = as.Date(character()), station = character(),
                 label = character(),
                 start_datetime = as.POSIXct(character()),
                 end_datetime = as.POSIXct(character()),
                 n_obs = integer(), output_file = character())
}

if (!"label" %in% names(log_inicial)) {
  log_inicial$label <- NA_character_
}

ui <- fluidPage(
  titlePanel("Seleção e exportação do intervalo do slug"),
  sidebarLayout(
    sidebarPanel(
      selectInput("arquivo", "Evento (source_file):", choices = arquivos_disponiveis),
      textInput("label", "Rótulo (opcional -- use para eventos com múltiplos slugs, ex: P, N):", value = ""),
      helpText("Arraste (brush) para selecionar o intervalo do slug.",
               "Dê duplo-clique numa área com brush para dar zoom;",
               "duplo-clique sem brush reseta o zoom."),
      actionButton("confirmar", "Confirmar e salvar", class = "btn-primary"),
      actionButton("reset_zoom", "Resetar zoom"),
      hr(),
      h4("Eventos já exportados:"),
      tableOutput("tabela_log")
    ),
    mainPanel(
      plotOutput("plot_serie",
                 brush = brushOpts(id = "brush", direction = "x", resetOnNew = FALSE),
                 dblclick = "plot_dblclick",
                 height = "500px"),
      verbatimTextOutput("brush_info")
    )
  )
)

server <- function(input, output, session) {
  
  log_exportacoes <- reactiveVal(log_inicial)
  zoom_range <- reactiveValues(x = NULL)
  
  dados_evento <- reactive({
    hobo_with_meta %>% filter(source_file == input$arquivo)
  })
  
  observeEvent(input$arquivo, {
    zoom_range$x <- NULL
  })
  
  observeEvent(input$plot_dblclick, {
    br <- input$brush
    if (!is.null(br)) {
      zoom_range$x <- c(br$xmin, br$xmax)
    } else {
      zoom_range$x <- NULL
    }
  })
  
  observeEvent(input$reset_zoom, {
    zoom_range$x <- NULL
  })
  
  output$plot_serie <- renderPlot({
    df <- dados_evento()
    
    p <- ggplot(df, aes(x = datetime, y = spc_uscm)) +
      geom_line() +
      geom_point(size = 0.6, alpha = 0.5) +
      geom_hline(aes(yintercept = background_spc), linetype = "dashed", color = "red") +
      labs(x = "Datetime", y = "SpC (µS/cm)", title = input$arquivo) +
      theme_minimal(base_size = 14)
    
    if (!is.null(zoom_range$x)) {
      x_lims <- as.POSIXct(zoom_range$x, origin = "1970-01-01", tz = "America/Sao_Paulo")
      p <- p + coord_cartesian(xlim = x_lims)
    }
    
    p
  })
  
  output$brush_info <- renderPrint({
    br <- input$brush
    if (is.null(br)) {
      cat("Nenhuma seleção feita ainda.")
    } else {
      cat("Início:", format(as.POSIXct(br$xmin, origin = "1970-01-01", tz = "America/Sao_Paulo")), "\n")
      cat("Fim:   ", format(as.POSIXct(br$xmax, origin = "1970-01-01", tz = "America/Sao_Paulo")), "\n")
    }
  })
  
  observeEvent(input$confirmar, {
    br <- input$brush
    req(br)
    
    start_dt <- as.POSIXct(br$xmin, origin = "1970-01-01", tz = "America/Sao_Paulo")
    end_dt   <- as.POSIXct(br$xmax, origin = "1970-01-01", tz = "America/Sao_Paulo")
    
    df_evento <- dados_evento()
    df_janela <- df_evento %>% filter(datetime >= start_dt, datetime <= end_dt)
    req(nrow(df_janela) > 0)
    
    stream_val  <- unique(df_evento$stream)
    date_val    <- unique(df_evento$date)
    station_val <- unique(df_evento$station)
    label_val   <- trimws(input$label)
    
    partes_nome <- c(
      stream_val,
      format(date_val, "%Y%m%d"),
      if (!is.na(station_val) && station_val != "") station_val else NULL,
      if (label_val != "") label_val else NULL
    )
    nome_base <- paste(partes_nome, collapse = "_")
    nome_arquivo <- paste0("slug_", nome_base, ".csv")
    caminho_completo <- file.path(pasta_saida, nome_arquivo)
    
    write_csv(df_janela, caminho_completo)
    
    nova_entrada <- tibble::tibble(
      source_file = input$arquivo, stream = stream_val, date = date_val,
      station = station_val, label = label_val,
      start_datetime = start_dt, end_datetime = end_dt,
      n_obs = nrow(df_janela), output_file = nome_arquivo
    )
    
    log_atual <- log_exportacoes() %>%
      filter(!(source_file == input$arquivo & label == label_val))
    novo_log <- bind_rows(log_atual, nova_entrada)
    log_exportacoes(novo_log)
    write_csv(novo_log, caminho_log)
    
    showNotification(paste("Salvo:", nome_arquivo, "-", nrow(df_janela), "observações"),
                     type = "message")
  })
  
  output$tabela_log <- renderTable({
    log_exportacoes() %>% arrange(stream, date, label)
  })
}

shinyApp(ui, server)