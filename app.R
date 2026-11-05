required_pkgs <- c("shiny", "ggplot2", "markdown",
                   "VineCopula", "data.table", "jsonlite", "Matrix", "MASS",
                   "copulapop")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Ontbrekende packages: ", paste(missing_pkgs, collapse = ", "), "\n",
       "Installeer met:\n",
       "  install.packages(c(\"shiny\", \"ggplot2\", \"markdown\", \"VineCopula\", \"data.table\", \"jsonlite\", \"Matrix\", \"MASS\"))\n",
       "  devtools::install_github(\"Robterheine/copulapop\")")
}
library(shiny)
library(ggplot2)
library(markdown)
library(copulapop)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

mat_choices <- c(
  "Farhan et al. 2024 — doi:10.1002/jcph.2484 (in-vivo UGT1A9, Hill)" = "farhan2024",
  "Badee et al. 2019 via Jo et al. 2021 scenario 2 — doi:10.1002/psp4.12577 (sigmoidaal + exponentieel)" = "badee2019"
)

datasets <- c(
  "NHANES (USA, algemene bevolking)" = "NHANES",
  "HSE (Engeland, algemene bevolking)" = "HSE",
  "RUMC (Nederland, ziekenhuispopulatie)" = "RUMC"
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
      numericInput("dose_mgkg", "Dosis per gift (mg/kg) — keerdosis, niet de dagdosis", value = 0.15, min = 0.01, step = 0.01),
      sliderInput("tau", "Dosisinterval (uur) tussen twee doses", min = 0, max = 24, value = 24, step = 1),
      sliderInput("n_doses", "Aantal doses tot steady state", min = 3, max = 14, value = 7, step = 1),
      sliderInput("iov_cv", "Inter-occasion variabiliteit op F1 / toegediende dosis (%CV per dosis)", min = 0, max = 100, value = 0, step = 5),
      sliderInput("age_range", "Leeftijdsrange (jaren)", min = 0, max = 18, value = c(0, 5), step = 0.5),
      selectInput("dataset", "Virtuele populatie (copulapop dataset)", choices = datasets, selected = "RUMC"),
      selectInput("mat_option", "Maturatiefunctie op niet-renale CL (UGT1A9)", choices = mat_choices, selected = "farhan2024"),
      selectInput("sex", "Geslacht", choices = c("Beide" = "both", "Jongen" = "male", "Meisje" = "female")),
      sliderInput("n_subj", "Aantal kinderen (N)", min = 50, max = 2000, value = 1000, step = 50),
      numericInput("seed", "Random seed", value = 42, min = 1),
      hr(),
      h4("Weergave"),
      checkboxInput("show_indiv", "Toon individuele curves (max. 30)", value = FALSE),
      checkboxInput("show_day1", "Toon mediaan dag 1 (eerste dosis)", value = FALSE),
      checkboxInput("log_y", "Logaritmische y-as", value = FALSE),
      hr(),
      h4("Gemeten data (overlay)"),
      fileInput("obs_file", "Upload CSV/TSV (kolommen TAD of TIME + DV of CONC) — of plak hieronder", accept = c(".csv", ".tsv", ".txt")),
      textAreaInput("obs_text", "Of plak gemeten data hier (bv. TAD,DV per regel)", rows = 5,
                    placeholder = "TAD,DV\n0.5,12.3\n1,88.1\n2,101.5\n24,3.9"),
      checkboxInput("show_obs", "Toon gemeten data als overlay op de curve", value = TRUE)
    ),
    mainPanel(
      width = 9,
      uiOutput("error_box"),
      p("Deze app simuleert de farmacokinetiek van dapagliflozine in een virtuele kinderpopulatie (gegenereerd met copulapop) op basis van het popPK-model van Penland et al. 2023. Stel links dosis (mg/kg), leeftijdsrange en dataset in en klik op ", strong("Simuleer"), ". Zie het tabblad 'Model & assumpties' voor de onderbouwing en bronnen."),
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "PK-curve (steady state)",
          plotOutput("pk_curve", height = "580px"),
          helpText("De band is een 90% predictie-interval over individuen (P5-P95). Alleen inter-individuele variabiliteit (IIV, plus optionele IOV per dosis) wordt gesimuleerd; residuele error ($SIGMA / EPS) wordt NIET gesimuleerd."),
          hr(),
          wellPanel(
            h4("Publication-ready figuur downloaden"),
            selectInput("fig_format", "Formaat", choices = c("PNG (raster)" = "png", "PDF (vector)" = "pdf", "SVG (vector)" = "svg", "TIFF (raster)" = "tiff"), selected = "png"),
            numericInput("fig_width", "Breedte (inch)", value = 7, min = 3, max = 20, step = 0.5),
            numericInput("fig_height", "Hoogte (inch)", value = 5, min = 3, max = 20, step = 0.5),
            numericInput("fig_dpi", "Resolutie (dpi; rasterformaten)", value = 600, min = 72, max = 1200, step = 50),
            downloadButton("dl_fig", "Download figuur (incl. gemeten overlay indien aangevinkt)")
          ),
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
          "AUC & Cmax vs nierfunctie / leeftijd",
          plotOutput("exposure_covariates", height = "720px")
        ),
        tabPanel(
          "Virtuele populatie",
          verbatimTextOutput("pop_info"),
          plotOutput("pop_plots", height = "360px"),
          h4("Eerste 100 kinderen uit de virtuele populatie (incl. de in de simulatie gebruikte eGFR)"),
          tableOutput("pop_table")
        ),
        tabPanel(
          "Model & assumpties",
          includeMarkdown("model_doc.md")
        ),
        tabPanel(
          "Over deze app",
          includeMarkdown("about.md")
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
        v <- pmin(pmax(input$age_range, lo), hi)
        updateSliderInput(session, "age_range", min = lo, max = hi, value = v)
      }
    }
  })

  do_simulate <- function() {
    err_msg(NULL)
    vals <- isolate(reactiveValuesToList(input))
    seed <- suppressWarnings(as.integer(vals$seed))
    if (is.na(seed)) stop("Random seed moet een geheel getal zijn.")
    tryCatch({
      if (is.null(vals$dose_mgkg) || is.na(vals$dose_mgkg) || vals$dose_mgkg <= 0)
        stop("Dosis per gift (mg/kg) moet een positief getal zijn.")
      withProgress(message = "Virtuele populatie genereren (copulapop) ...", value = 0.2, {
        pop <- generate_virtual_population(vals$dataset, vals$n_subj,
                                           vals$age_range, vals$sex, seed)
      })
      withProgress(message = "Monte Carlo PK-simulatie ...", value = 0.6, {
        res <- simulate_population_pk(pop, dose_mgkg = vals$dose_mgkg, tau = vals$tau,
                                      n_doses = vals$n_doses, iov_f1_cv = vals$iov_cv,
                                      mat_option = vals$mat_option,
                                      seed = as.integer((as.numeric(seed) + 1000) %% 2147483647))
      })
      sim_res(list(pop = pop, res = res))
      showNotification(
        sprintf("Simulatie voltooid: N = %d kinderen, dosis %.2f mg/kg per gift, IOV %d%% CV, seed %d",
                nrow(pop), vals$dose_mgkg, vals$iov_cv, seed),
        type = "message", duration = 5
      )
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (!is.null(conditionCall(e)))
        msg <- paste0(msg, "  [call: ", paste(deparse(conditionCall(e)), collapse = " "), "]")
      err_msg(msg)
      showNotification(paste("Fout:", msg), type = "error", duration = 10)
    })
  }

  observeEvent(input$run, do_simulate())
  do_simulate()  # eerste simulatie bij opstarten

  # ---------- PK-curve ----------
  parse_obs_text <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
    first <- strsplit(txt, "\n", fixed = TRUE)[[1]][1]
    sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
    df <- utils::read.table(text = txt, header = TRUE, sep = sep, stringsAsFactors = FALSE,
                            na.strings = c("NA", "", "."), check.names = FALSE)
    time_col <- intersect(c("TAD", "TIME", "Time", "time", "tad"), names(df))[1]
    conc_col <- intersect(c("DV", "CONC", "Conc", "conc", "CPRED", "IPRED"), names(df))[1]
    if (is.na(time_col) || is.na(conc_col)) return(NULL)
    if ("MDV" %in% names(df)) df <- df[is.na(df$MDV) | df$MDV == 0, , drop = FALSE]
    out <- data.frame(time = as.numeric(df[[time_col]]), conc = as.numeric(df[[conc_col]]))
    out[!is.na(out$time) & !is.na(out$conc) & out$conc >= 0, , drop = FALSE]
  }

  obs_df <- reactive({
    if (!is.null(input$obs_text) && nzchar(trimws(input$obs_text))) {
      d <- parse_obs_text(input$obs_text)
      if (is.null(d)) showNotification("Handmatig geplakte data: geen kolommen TAD/TIME + DV/CONC gevonden.", type = "error", duration = 10)
      return(d)
    }
    f <- input$obs_file
    if (is.null(f)) return(NULL)
    txt <- paste(readLines(f$datapath, warn = FALSE), collapse = "\n")
    d <- parse_obs_text(txt)
    if (is.null(d)) showNotification("Geupload bestand: kolom niet gevonden. Vereist een tijd-kolom (TAD of TIME) en een concentratie-kolom (DV of CONC).", type = "error", duration = 10)
    d
  })

  build_curve_plot <- function(s, show_indiv, show_day1, log_y, obs = NULL, pub = FALSE) {
    st <- s$res$settings
    cs <- curve_summary(s$res$curve_ss, s$res$times)
    p <- ggplot(cs, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.25) +
      geom_line(aes(y = q95), color = "#2c7fb8", linetype = "dashed") +
      geom_line(aes(y = q05), color = "#2c7fb8", linetype = "dashed") +
      geom_line(color = "#08306b", linewidth = 1.2) +
      do.call(labs, if (pub) list(
        x = "Time after dose (h)",
        y = "Plasma dapagliflozine concentration (ng/mL)",
        title = sprintf("Steady-state concentration-time profile (dosing interval %g h)", st$tau),
        subtitle = sprintf(
          "Median + 90%% prediction interval (P5-P95), N = %d, dose %.2f mg/kg per dose, interval %g h, %d doses, IOV on F1 %d%% CV, age %.1f-%.1f y",
          st$n, st$dose_mgkg, st$tau, st$n_doses, st$iov_f1_cv, min(s$pop$AGE), max(s$pop$AGE)
        )
      ) else list(
        x = "Tijd na dosis (uur)",
        y = "Dapagliflozine in plasma (ng/mL)",
        title = sprintf("Steady-state concentratie-tijdcurve (dosisinterval %g h)", st$tau),
        subtitle = sprintf(
          "Mediaan + 90%% predictie-interval (P5-P95), N = %d, dosis %.2f mg/kg, interval %g h, %d doses, IOV op F1 %d%% CV, leeftijd %.1f-%.1f jaar, maturatie: %s",
          st$n, st$dose_mgkg, st$tau, st$n_doses, st$iov_f1_cv, min(s$pop$AGE), max(s$pop$AGE), st$mat_label
        )
      ))
    if (isTRUE(show_indiv)) {
      idx <- sample(seq_len(nrow(s$pop)), min(30, nrow(s$pop)))
      indiv <- data.frame(
        t = rep(s$res$times, each = length(idx)),
        C = as.vector(s$res$curve_ss[idx, , drop = FALSE]),
        ID = rep(idx, times = length(s$res$times))
      )
      p <- p + geom_line(data = indiv, aes(x = t, y = C, group = ID),
                         color = "grey55", alpha = 0.35, linewidth = 0.4)
    }
    if (isTRUE(show_day1)) {
      cs1 <- curve_summary(s$res$curve_day1, s$res$times)
      p <- p + geom_line(data = cs1, aes(x = t, y = median),
                         color = "#d94801", linetype = "dotdash", linewidth = 1) +
        annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
                 label = if (pub) "orange = median day 1" else "oranje = mediaan dag 1", color = "#d94801", size = 3.5)
    }
    if (!is.null(obs) && nrow(obs) > 0) {
      p <- p + geom_point(data = obs, aes(x = time, y = conc), inherit.aes = FALSE,
                          color = "#D55E00", size = 2.4, alpha = 0.9) +
        annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6,
                 label = sprintf(if (pub) "orange points = observed data (n = %d)" else "oranje punten = gemeten data (n = %d)", nrow(obs)),
                 color = "#D55E00", size = 3.5)
    }
    if (isTRUE(log_y)) p <- p + scale_y_log10()
    if (pub) {
      p + theme_bw(base_size = 14) +
        theme(plot.title = element_text(face = "bold"),
              panel.grid.minor = element_blank())
    } else {
      p + theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold"))
    }
  }

  output$pk_curve <- renderPlot({
    s <- sim_res(); req(s)
    obs <- if (isTRUE(input$show_obs)) obs_df() else NULL
    build_curve_plot(s, input$show_indiv, input$show_day1, input$log_y, obs = obs, pub = FALSE)
  })

  output$dl_fig <- downloadHandler(
    filename = function() paste0("dapagliflozine_pk_curve.", input$fig_format),
    content = function(file) {
      s <- sim_res(); req(s)
      obs <- if (isTRUE(input$show_obs)) obs_df() else NULL
      p <- build_curve_plot(s, input$show_indiv, input$show_day1, input$log_y, obs = obs, pub = TRUE)
      fmt <- input$fig_format
      w <- input$fig_width; h <- input$fig_height; dpi <- input$fig_dpi
      if (fmt == "png") {
        grDevices::png(file, width = w, height = h, units = "in", res = dpi)
      } else if (fmt == "tiff") {
        grDevices::tiff(file, width = w, height = h, units = "in", res = dpi)
      } else if (fmt == "pdf") {
        grDevices::pdf(file, width = w, height = h)
      } else {
        grDevices::svg(file, width = w, height = h)
      }
      print(p)
      grDevices::dev.off()
    }
  )

  # ---------- Steady-state metrieken ----------
  output$ss_table <- renderTable({
    s <- sim_res(); req(s)
    summ <- summarize_metrics(s$res$metrics, tau = s$res$settings$tau)
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

  output$exposure_covariates <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    m$Geslacht <- factor(m$SEXISMALE, levels = c(0, 1), labels = c("meisje", "jongen"))
    long <- rbind(
      data.frame(metriek = "AUCtau,ss (ng*h/mL)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$AUCtau, Geslacht = m$Geslacht),
      data.frame(metriek = "AUCtau,ss (ng*h/mL)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$AUCtau, Geslacht = m$Geslacht),
      data.frame(metriek = "Cmax,ss (ng/mL)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$Cmax, Geslacht = m$Geslacht),
      data.frame(metriek = "Cmax,ss (ng/mL)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$Cmax, Geslacht = m$Geslacht)
    )
    ggplot(long, aes(x = x, y = y, color = Geslacht)) +
      geom_point(alpha = 0.4, size = 1.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1) +
      facet_grid(metriek ~ xvariabele, scales = "free") +
      scale_color_manual(values = c("meisje" = "#E69F00", "jongen" = "#0072B2")) +
      labs(x = NULL, y = NULL,
           title = "Steady-state exposure vs nierfunctie (eGFR) en leeftijd",
           subtitle = "Punten = individuele kinderen; lijn + band = LOESS-trend met 95% CI per panel") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom")
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
      sprintf("Gegenereerde virtuele populatie: N = %d (gevraagd: %d)", nrow(pop), input$n_subj),
      sprintf("Leeftijd: mediaan %.1f jaar (range %.1f-%.1f)",
              median(pop$AGE), min(pop$AGE), max(pop$AGE)),
      sprintf("Gewicht: mediaan %.1f kg (range %.1f-%.1f)",
              median(pop$WEIGHT), min(pop$WEIGHT), max(pop$WEIGHT)),
      sprintf("BMI: mediaan %.1f kg/m2", median(pop$BMI)),
      sprintf("eGFR: mediaan %.0f mL/min/1.73m2 (range %.0f-%.0f)",
              median(s$res$metrics$EGFR), min(s$res$metrics$EGFR), max(s$res$metrics$EGFR)),
      sprintf("eGFR-bron: %s", s$res$settings$egfr_source),
      sprintf("eGFR-eenheid: mL/min/1.73m2 — %s", s$res$settings$egfr_unit_note),
      "Renale CL geschaald met BSA/2 (BSA_ref = 2 m2; absoluut GFR = eGFR x BSA/1.73 wordt NIET als covariaat gebruikt) — zie tabblad Model & assumpties",
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
      data.frame(variabele = "BMI (kg/m2)", waarde = pop$BMI),
      data.frame(variabele = "eGFR (mL/min/1.73m2)", waarde = s$res$metrics$EGFR)
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
    if ("EGFR" %in% names(out)) out$EGFR <- NULL
    egfr_used <- s$res$metrics$EGFR
    idx <- if ("ID" %in% names(out)) match(out$ID, s$res$metrics$ID) else seq_len(nrow(out))
    out[["EGFR (mL/min/1.73m2)"]] <- round(egfr_used[idx], 1)
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
