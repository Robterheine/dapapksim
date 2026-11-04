library(shiny)
library(ggplot2)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

datasets <- c(
  "NHANES (USA, algemene bevolking)" = "NHANES",
  "KNHANES (Zuid-Korea, algemene bevolking)" = "KNHANES",
  "CHNS (China, algemene bevolking)" = "CHNS",
  "HSE (Engeland, algemene bevolking)" = "HSE",
  "TANZANIA (DHS survey)" = "TANZANIA",
  "RUMC (Nederland, ziekenhuispopulatie)" = "RUMC",
  "RUMC_IC (Nederland, IC)" = "RUMC_IC"
)

ui <- fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
  titlePanel("Dapagliflozine farmacokinetiek in kinderen — Monte Carlo populatiesimulatie"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Simulatie-instellingen"),
      actionButton("run", "Simuleer", icon = icon("play"), class = "btn-primary btn-lg", width = "100%"),
      hr(),
      sliderInput("dose_mgkg", "Dosis (mg/kg)", min = 0.02, max = 0.5, value = 0.15, step = 0.01),
      sliderInput("age_range", "Leeftijdsrange (jaren)", min = 0, max = 18, value = c(10, 17), step = 0.5),
      selectInput("dataset", "Virtuele populatie (copulapop dataset)", choices = datasets, selected = "NHANES"),
      selectInput("sex", "Geslacht", choices = c("Beide" = "both", "Jongen" = "male", "Meisje" = "female")),
      sliderInput("n_subj", "Aantal kinderen (N)", min = 50, max = 2000, value = 500, step = 50),
      sliderInput("n_doses", "Aantal doses 1x/dag (tot steady state)", min = 3, max = 14, value = 7, step = 1),
      sliderInput("iov_cv", "Inter-occasion variabiliteit op CL (%CV per dosis)", min = 0, max = 100, value = 0, step = 5),
      checkboxInput("ka_bsv", "BSV op Ka meenemen (omega^2 = 8.52, vastgezet in publicatie)", value = TRUE),
      numericInput("seed", "Random seed", value = 42, min = 1),
      textInput("data_path", "copulapop data-pad (leeg = automatisch)", value = "",
                placeholder = "bv. /pad/naar/copulapop/data"),
      hr(),
      h4("Weergave"),
      checkboxInput("show_indiv", "Toon individuele curves (max. 30)", value = FALSE),
      checkboxInput("show_day1", "Toon mediaan dag 1 (eerste dosis)", value = FALSE),
      checkboxInput("log_y", "Logaritmische y-as", value = FALSE)
    ),
    mainPanel(
      width = 9,
      uiOutput("error_box"),
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "PK-curve (steady state)",
          plotOutput("pk_curve", height = "580px"),
          helpText("De band is een 90% predictie-interval over individuen (P5-P95). Alleen inter-individuele variabiliteit (IIV, plus optionele IOV per dosis) wordt gesimuleerd; residuele error ($SIGMA / EPS) wordt NIET gesimuleerd."),
          hr(),
          downloadButton("dl_curve", "Download curve-samenvatting (CSV)"),
          downloadButton("dl_metrics", "Download steady-state metrieken per kind (CSV)")
        ),
        tabPanel(
          "Steady-state metrieken (Cmax / Cmin / AUC)",
          h4("Populatiesamenvatting steady state: mediaan [90% predictie-interval P5–P95]"),
          tableOutput("ss_table"),
          plotOutput("metric_plots", height = "380px"),
          plotOutput("auc_age", height = "380px")
        ),
        tabPanel(
          "Virtuele populatie",
          verbatimTextOutput("pop_info"),
          plotOutput("pop_plots", height = "360px"),
          h4("Eerste 100 kinderen uit de virtuele populatie"),
          tableOutput("pop_table")
        ),
        tabPanel(
          "Model & assumpties",
          includeMarkdown("model_doc.md")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  sim_res <- reactiveVal(NULL)
  err_msg <- reactiveVal(NULL)

  output$error_box <- renderUI({
    if (is.null(err_msg())) return(NULL)
    div(class = "alert alert-danger", strong("Fout: "), err_msg())
  })

  # Pas de leeftijdsschuif aan op het bereik van de gekozen dataset
  observeEvent(input$dataset, {
    info <- copulapop_dataset_info(input$dataset)
    if (!is.null(info)) {
      lo <- max(0, floor(info$age_range[1]))
      hi <- min(18, ceiling(info$age_range[2]))
      if (lo < hi) {
        v_lo <- max(lo, min(10, hi))
        v_hi <- min(hi, max(17, lo))
        updateSliderInput(session, "age_range", min = lo, max = hi, value = c(v_lo, v_hi))
      }
    }
  })

  do_simulate <- function() {
    err_msg(NULL)
    vals <- isolate(reactiveValuesToList(input))
    tryCatch({
      if (nzchar(trimws(vals$data_path))) {
        copulapop::set_data_path(trimws(vals$data_path))
      }
      withProgress(message = "Virtuele populatie genereren (copulapop) ...", value = 0.2, {
        pop <- generate_virtual_population(vals$dataset, vals$n_subj,
                                           vals$age_range, vals$sex, vals$seed)
      })
      withProgress(message = "Monte Carlo PK-simulatie ...", value = 0.6, {
        res <- simulate_population_pk(pop, dose_mgkg = vals$dose_mgkg, tau = 24,
                                      n_doses = vals$n_doses, iov_cv = vals$iov_cv,
                                      ka_bsv = vals$ka_bsv,
                                      seed = vals$seed + 1000L)
      })
      sim_res(list(pop = pop, res = res))
      showNotification(
        sprintf("Simulatie voltooid: N = %d kinderen, dosis %.2f mg/kg, IOV %d%% CV",
                nrow(pop), vals$dose_mgkg, vals$iov_cv),
        type = "message", duration = 5
      )
    }, error = function(e) {
      err_msg(conditionMessage(e))
      showNotification(paste("Fout:", conditionMessage(e)), type = "error", duration = 10)
    })
  }

  observeEvent(input$run, do_simulate())
  do_simulate()  # eerste simulatie bij opstarten

  # ---------- PK-curve ----------
  output$pk_curve <- renderPlot({
    s <- sim_res(); req(s)
    st <- s$res$settings
    cs <- curve_summary(s$res$curve_ss, s$res$times)
    p <- ggplot(cs, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.25) +
      geom_line(aes(y = q95), color = "#2c7fb8", linetype = "dashed") +
      geom_line(aes(y = q05), color = "#2c7fb8", linetype = "dashed") +
      geom_line(color = "#08306b", linewidth = 1.2) +
      labs(
        x = "Tijd na dosis (uur)",
        y = "Dapagliflozine in plasma (ng/mL)",
        title = "Steady-state concentratie-tijdcurve (1x/dag)",
        subtitle = sprintf(
          "Mediaan + 90%% predictie-interval (P5-P95), N = %d, dosis %.2f mg/kg, leeftijd %.1f-%.1f jaar, IOV op CL %d%% CV, %d doses",
          st$n, st$dose_mgkg, min(s$pop$AGE), max(s$pop$AGE), st$iov_cv, st$n_doses
        )
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
    if (isTRUE(input$show_indiv)) {
      idx <- sample(seq_len(nrow(s$pop)), min(30, nrow(s$pop)))
      indiv <- data.frame(
        t = rep(s$res$times, each = length(idx)),
        C = as.vector(s$res$curve_ss[idx, , drop = FALSE]),
        ID = rep(idx, times = length(s$res$times))
      )
      p <- p + geom_line(data = indiv, aes(x = t, y = C, group = ID),
                         color = "grey55", alpha = 0.35, linewidth = 0.4)
    }
    if (isTRUE(input$show_day1)) {
      cs1 <- curve_summary(s$res$curve_day1, s$res$times)
      p <- p + geom_line(data = cs1, aes(x = t, y = median),
                         color = "#d94801", linetype = "dotdash", linewidth = 1) +
        annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
                 label = "oranje = mediaan dag 1", color = "#d94801", size = 3.5)
    }
    if (isTRUE(input$log_y)) p <- p + scale_y_log10()
    p
  })

  # ---------- Steady-state metrieken ----------
  output$ss_table <- renderTable({
    s <- sim_res(); req(s)
    summ <- summarize_metrics(s$res$metrics)
    summ[] <- lapply(summ, function(x) if (is.numeric(x)) round(x, 2) else x)
    summ
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$metric_plots <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    met_long <- data.frame(
      metriek = rep(c("Cmax,ss (ng/mL)", "Cmin,ss (ng/mL)", "AUCtau,ss (ng*h/mL)", "Cavg,ss (ng/mL)"),
                    each = nrow(m)),
      waarde = c(m$Cmax, m$Cmin, m$AUCtau, m$Cavg)
    )
    ggplot(met_long, aes(x = metriek, y = waarde)) +
      geom_boxplot(fill = "#9ecae1", outlier.alpha = 0.4) +
      facet_wrap(~metriek, scales = "free", ncol = 4) +
      labs(x = NULL, y = NULL,
           title = "Verdeling van steady-state expositiemetrieken over de virtuele populatie") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$auc_age <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    m$Geslacht <- factor(m$SEXISMALE, levels = c(0, 1), labels = c("meisje", "jongen"))
    ggplot(m, aes(x = AGE, y = AUCtau, color = Geslacht)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(x = "Leeftijd (jaren)", y = "AUCtau,ss (ng*h/mL)",
           title = "Steady-state AUC vs leeftijd") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  # ---------- Virtuele populatie ----------
  output$pop_info <- renderText({
    s <- sim_res(); req(s)
    pop <- s$pop
    info <- copulapop_dataset_info(input$dataset)
    lines <- c(
      sprintf("Dataset: %s", input$dataset),
      if (!is.null(info)) sprintf("Bron: %s - %s (N=%s, leeftijd %g-%g jaar)",
                                  info$country, info$population,
                                  format(info$n_total, big.mark = ","),
                                  info$age_range[1], info$age_range[2]) else "",
      sprintf("Gegenereerde virtuele populatie: N = %d", nrow(pop)),
      sprintf("Leeftijd: mediaan %.1f jaar (range %.1f-%.1f)",
              median(pop$AGE), min(pop$AGE), max(pop$AGE)),
      sprintf("Gewicht: mediaan %.1f kg (range %.1f-%.1f)",
              median(pop$WEIGHT), min(pop$WEIGHT), max(pop$WEIGHT)),
      sprintf("BMI: mediaan %.1f kg/m2", median(pop$BMI)),
      sprintf("eGFR: mediaan %.0f mL/min/1.73m2 (range %.0f-%.0f)",
              median(s$res$metrics$EGFR), min(s$res$metrics$EGFR), max(s$res$metrics$EGFR)),
      sprintf("eGFR-bron: %s", s$res$settings$egfr_source),
      if (!is.null(info)) sprintf("Dataset bevat creatinine: %s",
                                  if (isTRUE(info$has_creatinine)) "ja" else "nee") else "", 
      sprintf("Geslacht: %d%% jongen", round(100 * mean(pop$SEXISMALE)))
    )
    paste(lines[nzchar(lines)], collapse = "\n")
  })

  output$pop_plots <- renderPlot({
    s <- sim_res(); req(s)
    pop <- s$pop
    pl <- rbind(
      data.frame(variabele = "Leeftijd (jaren)", waarde = pop$AGE),
      data.frame(variabele = "Gewicht (kg)", waarde = pop$WEIGHT),
      data.frame(variabele = "BMI (kg/m2)", waarde = pop$BMI)
    )
    ggplot(pl, aes(x = waarde)) +
      geom_histogram(bins = 30, fill = "#2c7fb8", color = "white") +
      facet_wrap(~variabele, scales = "free") +
      labs(x = NULL, y = "Aantal", title = "Virtuele populatie gegenereerd met copulapop") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$pop_table <- renderTable({
    s <- sim_res(); req(s)
    pop <- s$pop
    out <- as.data.frame(lapply(pop, function(x) if (is.numeric(x)) round(x, 2) else x))
    head(out, 100)
  })

  # ---------- Downloads ----------
  output$dl_curve <- downloadHandler(
    filename = "dapagliflozine_ss_curve.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      write.csv(curve_summary(s$res$curve_ss, s$res$times), file, row.names = FALSE)
    }
  )

  output$dl_metrics <- downloadHandler(
    filename = "dapagliflozine_ss_metrics.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      write.csv(s$res$metrics, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
