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
if (dir.exists("figures")) addResourcePath("figures", "figures")
if (dir.exists("docs")) addResourcePath("docs", "docs")

mat_choices <- c(
  "Farhan et al. 2024 — jcph.2484 (in-vivo UGT1A9, standaard)" = "farhan2024",
  "Badee et al. 2019 via Jo et al. 2021 — psp4.12577 (PBPK s2)" = "badee2019",
  "Simcyp standaard via Jo et al. 2021 — psp4.12577 (PBPK s1)" = "simcyp_s1"
)

ka_choices <- c(
  "Gekalibreerd klinisch (ω² = 0.36; standaard)" = "calibrated",
  "Penland 2023 bronmodel (ω² = 8.52 FIX)" = "penland_raw",
  "Geen Ka-variabiliteit (ω² = 0)" = "none"
)

datasets <- c(
  "RUMC (Nederland, ziekenhuispopulatie, 0–100 j, met CREAT)" = "RUMC",
  "NHANES (USA, algemene bevolking, 0–80 j)" = "NHANES",
  "HSE (Engeland, algemene bevolking, 2–90 j)" = "HSE"
)

ui <- fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
  titlePanel("Dapagliflozine farmacokinetiek & fysiologische extrapolatie — Monte Carlo populatiesimulatie"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Simulatie-instellingen"),
      actionButton("run", "Simuleer", icon = icon("play"), class = "btn-primary btn-lg", width = "100%"),
      hr(),
      radioButtons(
        "dose_mode", "Doseringswijze",
        choices = c("Op gewicht (mg/kg)" = "mgkg", "Vaste dosis (mg)" = "mg"),
        selected = "mgkg", inline = TRUE
      ),
      conditionalPanel(
        condition = "input.dose_mode == 'mgkg'",
        numericInput("dose_mgkg", "Dosis per gift (mg/kg) — keerdosis, zonder plafond",
                     value = 0.15, min = 0.001, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.dose_mode == 'mg'",
        numericInput("dose_mg", "Vaste dosis per gift (mg) — keerdosis, zonder plafond",
                     value = 10, min = 0.01, step = 0.5)
      ),
      sliderInput("tau", "Dosisinterval (uur) tussen twee doses", min = 1, max = 48, value = 24, step = 1),
      sliderInput("n_doses", "Aantal doses (1 = eenmalige gift)", min = 1, max = 28, value = 7, step = 1),
      checkboxInput("steady_state", "Directe analytische steady state (bij IOV = 0%)", value = FALSE),
      sliderInput("iov_cv", "Inter-occasion variabiliteit op F1 / toegediende dosis (%CV per dosis)",
                  min = 0, max = 100, value = 0, step = 5),
      tags$label("Leeftijdsrange (jaren; bv. 0.08 = 1 mnd)"),
      fluidRow(
        column(6, numericInput("age_min", "Minimum", value = 0, min = 0, max = 100, step = 0.1)),
        column(6, numericInput("age_max", "Maximum", value = 5, min = 0, max = 100, step = 0.5))
      ),
      uiOutput("dataset_age_hint"),
      selectInput("dataset", "Virtuele populatie (copulapop dataset)", choices = datasets, selected = "RUMC"),
      selectInput("mat_option", "Maturatiefunctie op hepatische + renale UGT1A9", choices = mat_choices, selected = "farhan2024"),
      selectInput("ka_mode", "Absorptie-variabiliteit (BSV op Ka)", choices = ka_choices, selected = "calibrated"),
      selectInput("sex", "Geslacht", choices = c("Beide" = "both", "Man / jongen" = "male", "Vrouw / meisje" = "female")),
      sliderInput("n_subj", "Aantal virtuele individuen (N)", min = 50, max = 2000, value = 1000, step = 50),
      numericInput("seed", "Random seed", value = 42, min = 1),
      hr(),
      h4("Weergave"),
      checkboxInput("show_indiv", "Toon individuele curves (max. 30)", value = FALSE),
      checkboxInput("show_day1", "Toon mediaan dag 1 (eerste dosis)", value = FALSE),
      checkboxInput("log_y", "Logaritmische y-as", value = FALSE),
      hr(),
      h4("Gemeten data (overlay)"),
      fileInput("obs_file", "Upload CSV/TSV (kolommen TAD of TIME + DV of CONC) — of plak hieronder",
                accept = c(".csv", ".tsv", ".txt")),
      textAreaInput("obs_text", "Of plak gemeten data hier (bv. TAD,DV per regel)", rows = 5,
                    placeholder = "TAD,DV\n0.5,12.3\n1,88.1\n2,101.5\n24,3.9"),
      checkboxInput("show_obs", "Toon gemeten data als overlay op de curve", value = TRUE)
    ),
    mainPanel(
      width = 9,
      uiOutput("error_box"),
      uiOutput("warn_box"),
      p("Deze app simuleert de farmacokinetiek van dapagliflozine in een virtuele populatie (gegenereerd met copulapop) op basis van het fysiologisch geëxtrapoleerde popPK-model van Penland et al. 2023 (met FFM/NFM-allometrie, albumine-ontogenie, Rhodin/CKiD-U25 nierfunctie en hepatische + renale UGT1A9-maturatie). Stel links dosis, leeftijdsrange (van pasgeborene tot oudere) en dataset in en klik op ", strong("Simuleer"), "."),
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "PK-curve",
          plotOutput("pk_curve", height = "580px"),
          helpText("De band is een 90% predictie-interval over individuen (P5-P95). Alleen inter-individuele variabiliteit (IIV, plus optionele IOV per dosis) wordt gesimuleerd; residuele error ($SIGMA / EPS) wordt NIET gesimuleerd."),
          hr(),
          wellPanel(
            h4("Publication-ready figuur downloaden"),
            selectInput("fig_format", "Formaat",
                        choices = c("PNG (raster)" = "png", "PDF (vector)" = "pdf", "SVG (vector)" = "svg", "TIFF (raster)" = "tiff"),
                        selected = "png"),
            numericInput("fig_width", "Breedte (inch)", value = 7, min = 2, max = 24, step = 0.5),
            numericInput("fig_height", "Hoogte (inch)", value = 5, min = 2, max = 24, step = 0.5),
            numericInput("fig_dpi", "Resolutie (dpi; rasterformaten)", value = 600, min = 72, max = 1200, step = 50),
            downloadButton("dl_fig", "Download figuur (incl. gemeten overlay indien aangevinkt)")
          ),
          hr(),
          downloadButton("dl_curve", "Download curve-samenvatting (CSV)"),
          downloadButton("dl_metrics", "Download individuele metrieken (CSV)")
        ),
        tabPanel(
          "Expositiemetrieken (Cmax / Cmin / AUC / AUCu)",
          h4("Populatiesamenvatting: mediaan, geometrisch gemiddelde en [90% predictie-interval P5–P95]"),
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
          h4("Eerste 100 individuen uit de virtuele populatie (incl. FFM en de in de simulatie gebruikte eGFR)"),
          tableOutput("pop_table")
        ),
        tabPanel(
          "Model & assumpties",
          withMathJax(div(class = "md-doc", includeMarkdown("model_doc.md")))
        ),
        tabPanel(
          "Validatie",
          div(
            class = "md-doc",
            wellPanel(
              strong("Volledig validatierapport (Engels, incl. alle tabellen en figuren): "),
              tags$a("Open PDF-rapport (VALIDATION_REPORT.pdf)", href = "docs/VALIDATION_REPORT.pdf",
                     target = "_blank", class = "btn btn-default btn-sm"),
              " ",
              tags$a("Open HTML-rapport (VALIDATION_REPORT.html)", href = "docs/VALIDATION_REPORT.html",
                     target = "_blank", class = "btn btn-default btn-sm")
            ),
            includeMarkdown("VALIDATION.md"),
            hr(),
            h3("Gedetailleerde validatiefiguren (uit het volledige rapport)"),
            tags$img(src = "figures/fig1_penland_egfr.svg", alt = "Figure 1: Penland 2023 CKD strata"),
            tags$img(src = "figures/fig2_tirucherai_dose.svg", alt = "Figure 2: Tirucherai 2016 paediatric T2DM"),
            tags$img(src = "figures/fig3_pbpk_age.svg", alt = "Figure 3: Jo 2021 paediatric PBPK"),
            tags$img(src = "figures/fig3_curves.svg", alt = "Figure 4: Paediatric PK profiles"),
            tags$img(src = "figures/fig5_maturation.svg", alt = "Figure 5: Developmental ontogeny")
          )
        ),
        tabPanel(
          "Over deze app",
          div(class = "md-doc", includeMarkdown("about.md"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  sim_res  <- reactiveVal(NULL)
  err_msg  <- reactiveVal(NULL)
  warn_msg <- reactiveVal(NULL)

  output$error_box <- renderUI({
    if (is.null(err_msg())) return(NULL)
    div(class = "alert alert-danger", strong("Fout: "), err_msg())
  })

  output$warn_box <- renderUI({
    w <- warn_msg()
    if (is.null(w) || length(w) == 0) return(NULL)
    div(class = "alert alert-warning",
        lapply(w, function(m) div(strong("Let op: "), m)))
  })

  output$dataset_age_hint <- renderUI({
    info <- copulapop_dataset_info(input$dataset)
    if (is.null(info)) return(NULL)
    helpText(sprintf("Beschikbaar in %s: %g tot %g jaar.",
                     input$dataset, info$age_range[1], info$age_range[2]))
  })

  # Pas leeftijdsgrenzen alleen aan als ze buiten het bereik van de nieuw gekozen dataset vallen
  observeEvent(input$dataset, ignoreInit = TRUE, {
    info <- copulapop_dataset_info(input$dataset)
    if (!is.null(info)) {
      ds_lo <- as.numeric(info$age_range[1])
      ds_hi <- as.numeric(info$age_range[2])
      cur_lo <- suppressWarnings(as.numeric(input$age_min))
      cur_hi <- suppressWarnings(as.numeric(input$age_max))
      if (is.na(cur_lo)) cur_lo <- ds_lo
      if (is.na(cur_hi)) cur_hi <- min(5, ds_hi)
      new_lo <- min(max(cur_lo, ds_lo), ds_hi)
      new_hi <- min(max(cur_hi, new_lo), ds_hi)
      if (new_lo == new_hi && new_hi < ds_hi) new_hi <- min(ds_hi, new_lo + 1)
      updateNumericInput(session, "age_min", min = ds_lo, max = ds_hi, value = new_lo)
      updateNumericInput(session, "age_max", min = ds_lo, max = ds_hi, value = new_hi)
    }
  })

  do_simulate <- function() {
    err_msg(NULL)
    warn_msg(NULL)
    vals <- isolate(reactiveValuesToList(input))
    seed <- suppressWarnings(as.integer(vals$seed))
    if (length(seed) != 1L || is.na(seed) || seed < 1L) {
      err_msg("Random seed moet een positief geheel getal zijn.")
      return(invisible(NULL))
    }

    warns <- character(0)
    tryCatch({
      d_mode <- if (length(vals$dose_mode) == 1L && nzchar(vals$dose_mode)) vals$dose_mode else "mgkg"
      d_mgkg <- NULL
      d_mg   <- NULL
      if (d_mode == "mg") {
        val_mg <- suppressWarnings(as.numeric(vals$dose_mg))
        if (length(val_mg) != 1L || is.na(val_mg) || val_mg <= 0) {
          stop("Vaste dosis per gift (mg) moet een positief getal zijn.")
        }
        d_mg <- val_mg
      } else {
        val_mgkg <- suppressWarnings(as.numeric(vals$dose_mgkg))
        if (length(val_mgkg) != 1L || is.na(val_mgkg) || val_mgkg <= 0) {
          stop("Dosis per gift (mg/kg) moet een positief getal zijn.")
        }
        d_mgkg <- val_mgkg
      }

      a_min <- suppressWarnings(as.numeric(vals$age_min))
      a_max <- suppressWarnings(as.numeric(vals$age_max))
      if (length(a_min) != 1L || length(a_max) != 1L || is.na(a_min) || is.na(a_max) || a_min < 0 || a_max < 0) {
        stop("Minimum- en maximumleeftijd moeten niet-negatieve getallen zijn.")
      }
      if (a_min > a_max) {
        stop("Minimumleeftijd mag niet groter zijn dan maximumleeftijd.")
      }

      ds_name <- if (length(vals$dataset) == 1L && nzchar(vals$dataset)) vals$dataset else "RUMC"
      n_subj  <- if (length(vals$n_subj) == 1L && !is.na(vals$n_subj)) as.integer(vals$n_subj) else 1000L
      sex_val <- if (length(vals$sex) == 1L && nzchar(vals$sex)) vals$sex else "both"
      tau_val <- if (length(vals$tau) == 1L && !is.na(vals$tau)) as.numeric(vals$tau) else 24
      nd_val  <- if (length(vals$n_doses) == 1L && !is.na(vals$n_doses)) as.integer(vals$n_doses) else 7L
      iov_val <- if (length(vals$iov_cv) == 1L && !is.na(vals$iov_cv)) as.numeric(vals$iov_cv) else 0
      mat_val <- if (length(vals$mat_option) == 1L && nzchar(vals$mat_option)) vals$mat_option else "farhan2024"
      ka_m    <- if (length(vals$ka_mode) == 1L && nzchar(vals$ka_mode)) vals$ka_mode else "calibrated"

      ds_info <- copulapop_dataset_info(ds_name)

      pop <- withCallingHandlers(
        withProgress(message = "Virtuele populatie genereren (copulapop) ...", value = 0.25, {
          generate_virtual_population(ds_name, n_subj, c(a_min, a_max), sex_val, seed)
        }),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      use_ss <- isTRUE(vals$steady_state) && iov_val == 0
      if (isTRUE(vals$steady_state) && iov_val > 0) {
        warns <- c(warns, sprintf(
          "Directe analytische steady state is alleen gedefinieerd bij IOV = 0%%; er is gesimuleerd met %d doses en %g%% IOV.",
          nd_val, iov_val
        ))
      }

      res <- withProgress(message = "Monte Carlo PK-simulatie ...", value = 0.70, {
        simulate_population_pk(
          pop,
          dose_mgkg    = d_mgkg,
          dose_mg      = d_mg,
          tau          = tau_val,
          n_doses      = nd_val,
          iov_f1_cv    = iov_val,
          ka_bsv       = (ka_m != "none"),
          ka_mode      = ka_m,
          mat_option   = mat_val,
          steady_state = use_ss,
          seed         = as.integer((as.numeric(seed) + 1000) %% 2147483647)
        )
      })

      if (!use_ss && nd_val > 1L && is.finite(res$settings$ss_attainment) && res$settings$ss_attainment < 0.90) {
        warns <- c(warns, sprintf(
          "Na %d doses (interval %g u) is mediaan %.0f%% van steady-state blootstelling bereikt (mediane t1/2,beta = %.1f u). Verhoog 'Aantal doses' of vink 'Directe analytische steady state' aan voor volledige steady state.",
          nd_val, tau_val, 100 * res$settings$ss_attainment, stats::median(res$metrics$t_half)
        ))
      }

      # Deterministische selectie van individuele curves zodat scherm en download identiek zijn
      set.seed(as.integer((as.numeric(seed) + 2026) %% 2147483647))
      indiv_idx <- sample(seq_len(nrow(pop)), min(30L, nrow(pop)))

      sim_res(list(
        pop          = pop,
        res          = res,
        dataset      = ds_name,
        dataset_info = ds_info,
        n_req        = n_subj,
        indiv_idx    = indiv_idx
      ))
      warn_msg(unique(warns))

      dose_txt <- if (d_mode == "mg") sprintf("%.3g mg", d_mg) else sprintf("%.3g mg/kg", d_mgkg)
      showNotification(
        sprintf("Simulatie voltooid: N = %d, dosis %s per gift, IOV %g%% CV, seed %d",
                nrow(pop), dose_txt, iov_val, seed),
        type = "message", duration = 5
      )
    }, error = function(e) {
      msg <- conditionMessage(e)
      err_msg(msg)
      showNotification(paste("Fout:", msg), type = "error", duration = 10)
    })
  }

  observeEvent(input$run, do_simulate())
  observeEvent(input$seed, {
    if (is.null(sim_res())) do_simulate()
  }, once = TRUE, ignoreInit = FALSE)

  # ---------- PK-curve & gemeten data ----------
  parse_obs_text <- function(txt) {
    if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
    tryCatch({
      lines <- strsplit(txt, "\r?\n")[[1]]
      lines <- lines[nzchar(trimws(lines))]
      if (length(lines) < 2) return(NULL)
      first <- lines[1]
      sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
      dec_txt <- if (sep %in% c(";", "\t")) gsub(",", ".", paste(lines, collapse = "\n"), fixed = TRUE) else paste(lines, collapse = "\n")
      df <- utils::read.table(
        text = dec_txt, header = TRUE, sep = sep, stringsAsFactors = FALSE,
        na.strings = c("NA", "", "."), check.names = FALSE, fill = TRUE, comment.char = "#"
      )
      time_col <- intersect(c("TAD", "TIME", "Time", "time", "tad"), names(df))[1]
      conc_col <- intersect(c("DV", "CONC", "Conc", "conc", "CPRED", "IPRED"), names(df))[1]
      if (is.na(time_col) || is.na(conc_col)) return(NULL)
      if ("MDV" %in% names(df)) df <- df[is.na(df$MDV) | df$MDV == 0, , drop = FALSE]
      out <- suppressWarnings(data.frame(
        time = as.numeric(df[[time_col]]),
        conc = as.numeric(df[[conc_col]])
      ))
      out <- out[!is.na(out$time) & !is.na(out$conc) & out$conc >= 0, , drop = FALSE]
      if (nrow(out) == 0) NULL else out
    }, error = function(e) NULL)
  }

  obs_df <- reactive({
    if (!is.null(input$obs_text) && nzchar(trimws(input$obs_text))) {
      d <- parse_obs_text(input$obs_text)
      if (is.null(d)) {
        showNotification("Handmatig geplakte data: geen geldige kolommen TAD/TIME + DV/CONC gevonden.",
                         type = "error", duration = 8)
      }
      return(d)
    }
    f <- input$obs_file
    if (is.null(f)) return(NULL)
    txt <- tryCatch(paste(readLines(f$datapath, warn = FALSE), collapse = "\n"), error = function(e) "")
    d <- parse_obs_text(txt)
    if (is.null(d)) {
      showNotification("Geupload bestand: vereist een tijd-kolom (TAD of TIME) en een concentratie-kolom (DV of CONC).",
                       type = "error", duration = 8)
    }
    d
  })

  build_curve_plot <- function(s, show_indiv, show_day1, log_y, obs = NULL, pub = FALSE) {
    st <- s$res$settings
    cs <- s$res$cs_ss
    if (isTRUE(log_y)) {
      cs <- cs[cs$median > 0 & cs$q95 > 0, , drop = FALSE]
      cs$q05 <- ifelse(cs$q05 > 0, cs$q05, pmin(cs$median * 0.1, 1e-3))
    }

    dose_lbl_nl <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg", st$dose_mgkg)
    dose_lbl_en <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg per dose", st$dose_mgkg)
    regimen_nl  <- if (isTRUE(st$steady_state)) "analytische steady state" else if (st$n_doses == 1) "eenmalige dosis" else sprintf("dosis %d", st$n_doses)
    regimen_en  <- if (isTRUE(st$steady_state)) "analytical steady state" else if (st$n_doses == 1) "single dose" else sprintf("dose %d", st$n_doses)

    p <- ggplot(cs, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.25) +
      geom_line(aes(y = q95), color = "#2c7fb8", linetype = "dashed") +
      geom_line(aes(y = q05), color = "#2c7fb8", linetype = "dashed") +
      geom_line(color = "#08306b", linewidth = 1.2) +
      do.call(labs, if (pub) list(
        x = "Time after dose (h)",
        y = "Plasma dapagliflozin concentration (ng/mL)",
        title = sprintf("Concentration-time profile (%s, interval %g h)", regimen_en, st$tau),
        subtitle = sprintf(
          "Median + 90%% PI (P5-P95), N = %d (%s), dose %s, IOV %d%% CV, age %.2f-%.2f y, UGT1A9: %s",
          st$n, s$dataset, dose_lbl_en, st$iov_f1_cv, min(s$pop$AGE), max(s$pop$AGE), st$mat_label
        )
      ) else list(
        x = "Tijd na dosis (uur)",
        y = "Dapagliflozine in plasma (ng/mL)",
        title = sprintf("Concentratie-tijdcurve (%s, interval %g u)", regimen_nl, st$tau),
        subtitle = sprintf(
          "Mediaan + 90%% PI (P5-P95), N = %d (%s), dosis %s, IOV %d%% CV, leeftijd %.2f-%.2f j, maturatie: %s",
          st$n, s$dataset, dose_lbl_nl, st$iov_f1_cv, min(s$pop$AGE), max(s$pop$AGE), st$mat_label
        )
      ))

    if (isTRUE(show_indiv)) {
      idx <- s$indiv_idx
      indiv <- data.frame(
        t  = rep(s$res$times, each = length(idx)),
        C  = as.vector(s$res$curve_ss[idx, , drop = FALSE]),
        ID = rep(idx, times = length(s$res$times))
      )
      if (isTRUE(log_y)) indiv <- indiv[indiv$C > 0, , drop = FALSE]
      if (nrow(indiv) > 0) {
        p <- p + geom_line(data = indiv, aes(x = t, y = C, group = ID),
                           color = "grey55", alpha = 0.35, linewidth = 0.4)
      }
    }
    if (isTRUE(show_day1)) {
      cs1 <- s$res$cs_day1
      if (isTRUE(log_y)) cs1 <- cs1[cs1$median > 0, , drop = FALSE]
      if (nrow(cs1) > 0) {
        p <- p + geom_line(data = cs1, aes(x = t, y = median),
                           color = "#d94801", linetype = "dotdash", linewidth = 1) +
          annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
                   label = if (pub) "orange = median day 1" else "oranje = mediaan dag 1",
                   color = "#d94801", size = 3.5)
      }
    }
    if (!is.null(obs) && nrow(obs) > 0) {
      obs_plot <- if (isTRUE(log_y)) obs[obs$conc > 0, , drop = FALSE] else obs
      if (nrow(obs_plot) > 0) {
        p <- p + geom_point(data = obs_plot, aes(x = time, y = conc), inherit.aes = FALSE,
                            color = "#D55E00", size = 2.4, alpha = 0.9) +
          annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6,
                   label = sprintf(if (pub) "orange points = observed data (n = %d)" else "oranje punten = gemeten data (n = %d)", nrow(obs_plot)),
                   color = "#D55E00", size = 3.5)
      }
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
      w   <- suppressWarnings(as.numeric(input$fig_width));  if (is.na(w)) w <- 7;  w <- min(max(w, 2), 24)
      h   <- suppressWarnings(as.numeric(input$fig_height)); if (is.na(h)) h <- 5;  h <- min(max(h, 2), 24)
      dpi <- suppressWarnings(as.numeric(input$fig_dpi));    if (is.na(dpi)) dpi <- 600; dpi <- min(max(dpi, 72), 1200)
      if (fmt == "png") {
        grDevices::png(file, width = w, height = h, units = "in", res = dpi)
      } else if (fmt == "tiff") {
        grDevices::tiff(file, width = w, height = h, units = "in", res = dpi)
      } else if (fmt == "pdf") {
        grDevices::pdf(file, width = w, height = h)
      } else {
        grDevices::svg(file, width = w, height = h)
      }
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      print(p)
    }
  )

  # ---------- Expositiemetrieken ----------
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
      metriek = factor(
        rep(c("Cmax (ng/mL)", "Cmin (ng/mL)", "AUCtau totaal (ng*h/mL)", "AUCu ongebonden (ng*h/mL)"),
            each = nrow(m)),
        levels = c("Cmax (ng/mL)", "Cmin (ng/mL)", "AUCtau totaal (ng*h/mL)", "AUCu ongebonden (ng*h/mL)")
      ),
      waarde = c(m$Cmax, m$Cmin, m$AUCtau, m$AUCu)
    )
    ggplot(met_long, aes(x = metriek, y = waarde)) +
      geom_boxplot(fill = "#9ecae1", outlier.alpha = 0.4) +
      facet_wrap(~metriek, scales = "free", ncol = 4) +
      labs(x = NULL, y = NULL,
           title = "Verdeling van expositiemetrieken (totaal en ongebonden) over de virtuele populatie") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  })

  sex_labels <- function(ages) {
    if (max(ages, na.rm = TRUE) >= 18) c("vrouw / meisje", "man / jongen") else c("meisje", "jongen")
  }

  output$auc_age <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    lbls <- sex_labels(m$AGE)
    m$Geslacht <- factor(m$SEXISMALE, levels = c(0, 1), labels = lbls)
    ggplot(m, aes(x = AGE, y = AUCtau, color = Geslacht)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(x = "Leeftijd (jaren)", y = "AUCtau (ng*h/mL)",
           title = "Blootstelling (AUCtau) vs leeftijd") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$exposure_covariates <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    lbls <- sex_labels(m$AGE)
    m$Geslacht <- factor(m$SEXISMALE, levels = c(0, 1), labels = lbls)
    long <- rbind(
      data.frame(metriek = "AUCtau (ng*h/mL)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$AUCtau, Geslacht = m$Geslacht),
      data.frame(metriek = "AUCtau (ng*h/mL)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$AUCtau, Geslacht = m$Geslacht),
      data.frame(metriek = "Cmax (ng/mL)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$Cmax, Geslacht = m$Geslacht),
      data.frame(metriek = "Cmax (ng/mL)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$Cmax, Geslacht = m$Geslacht)
    )
    ggplot(long, aes(x = x, y = y, color = Geslacht)) +
      geom_point(alpha = 0.4, size = 1.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1) +
      facet_grid(metriek ~ xvariabele, scales = "free") +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(x = NULL, y = NULL,
           title = "Blootstelling vs nierfunctie (eGFR) en leeftijd",
           subtitle = "Punten = virtuele individuen; lijn + band = LOESS-trend met 95% CI per panel") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom")
  })

  # ---------- Virtuele populatie ----------
  output$pop_info <- renderText({
    s <- sim_res(); req(s)
    pop  <- s$pop
    m    <- s$res$metrics
    info <- s$dataset_info
    lines <- c(
      sprintf("Dataset: %s", s$dataset),
      if (!is.null(info)) sprintf("Bron: %s - %s (N=%s, leeftijd %g-%g jaar)",
                                  info$country, info$population,
                                  format(info$n_total, big.mark = ","),
                                  info$age_range[1], info$age_range[2]) else "",
      sprintf("Gegenereerde virtuele populatie: N = %d (gevraagd: %d)", nrow(pop), s$n_req),
      sprintf("Leeftijd: mediaan %.2f jaar (range %.3f-%.2f)",
              stats::median(pop$AGE), min(pop$AGE), max(pop$AGE)),
      sprintf("Gewicht (BWT): mediaan %.1f kg (range %.1f-%.1f) | Vetvrije massa (FFM): mediaan %.1f kg",
              stats::median(pop$WEIGHT), min(pop$WEIGHT), max(pop$WEIGHT), stats::median(m$FFM)),
      sprintf("BMI: mediaan %.1f kg/m2 | Ongebonden fractie fu: mediaan %.3f (range %.3f-%.3f)",
              stats::median(pop$BMI), stats::median(m$FU), min(m$FU), max(m$FU)),
      sprintf("eGFR: mediaan %.0f mL/min/1.73m2 (range %.0f-%.0f)",
              stats::median(m$EGFR), min(m$EGFR), max(m$EGFR)),
      sprintf("eGFR-bron: %s", s$res$settings$egfr_source),
      sprintf("eGFR-eenheid: mL/min/1.73m2 — %s", s$res$settings$egfr_unit_note),
      if (!is.null(info)) sprintf("Dataset bevat serumcreatinine: %s",
                                  if (isTRUE(info$has_creatinine)) "ja" else "nee") else "",
      sprintf("Geslacht: %d%% mannelijk", round(100 * mean(pop$SEXISMALE)))
    )
    paste(lines[nzchar(lines)], collapse = "\n")
  })

  output$pop_plots <- renderPlot({
    s <- sim_res(); req(s)
    pop <- s$pop
    m   <- s$res$metrics
    pl <- rbind(
      data.frame(variabele = "Leeftijd (jaren)", waarde = pop$AGE),
      data.frame(variabele = "Gewicht (kg)", waarde = pop$WEIGHT),
      data.frame(variabele = "Vetvrije massa FFM (kg)", waarde = m$FFM),
      data.frame(variabele = "eGFR (mL/min/1.73m2)", waarde = m$EGFR)
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
    m   <- s$res$metrics
    out <- as.data.frame(lapply(pop, function(x) if (is.numeric(x)) round(x, 2) else x))
    if ("EGFR" %in% names(out)) out$EGFR <- NULL
    idx <- if ("ID" %in% names(out)) match(out$ID, m$ID) else seq_len(nrow(out))
    out[["FFM (kg)"]] <- round(m$FFM[idx], 2)
    out[["fu"]] <- round(m$FU[idx], 3)
    out[["EGFR (mL/min/1.73m2)"]] <- round(m$EGFR[idx], 1)
    head(out, 100)
  })

  # ---------- Downloads ----------
  output$dl_curve <- downloadHandler(
    filename = "dapagliflozine_pk_curve.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      write.csv(s$res$cs_ss, file, row.names = FALSE)
    }
  )

  output$dl_metrics <- downloadHandler(
    filename = "dapagliflozine_pk_metrics.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      write.csv(s$res$metrics, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
