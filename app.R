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
library(grid)

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)
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

presets_list <- hf_presets()
hf_preset_choices <- stats::setNames(names(presets_list), vapply(presets_list, `[[`, character(1), "label"))

side_by_side_plots <- function(p_left, p_right) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_left,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

ui <- fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
  titlePanel("Dapagliflozine PK/PD & fysiologische dosisextrapolatie — incl. pediatrisch hartfalen"),
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
      checkboxInput("steady_state", "Directe analytische steady state (bij IOV = 0%)", value = TRUE),
      sliderInput("iov_cv", "Inter-occasion variabiliteit op F1 / toegediende dosis (%CV per dosis)",
                  min = 0, max = 100, value = 0, step = 5),
      tags$label("Leeftijdsrange (jaren; bv. 0.08 = 1 mnd, 0.5 = 6 mnd)"),
      fluidRow(
        column(6, numericInput("age_min", "Minimum", value = 0.1, min = 0, max = 100, step = 0.05)),
        column(6, numericInput("age_max", "Maximum", value = 1.0, min = 0, max = 100, step = 0.25))
      ),
      uiOutput("dataset_age_hint"),
      selectInput("dataset", "Virtuele populatie (copulapop dataset)", choices = datasets, selected = "RUMC"),
      selectInput("mat_option", "Maturatiefunctie op hepatische + renale UGT1A9", choices = mat_choices, selected = "farhan2024"),
      selectInput("ka_mode", "Absorptie-variabiliteit (BSV op Ka)", choices = ka_choices, selected = "calibrated"),
      div(
        class = "hf-sidebar-box",
        h4("Hartfalen-fysiologie & sonde"),
        selectInput("hf_preset", "Klinisch hartfalenprofiel", choices = hf_preset_choices, selected = "standard"),
        uiOutput("hf_preset_note"),
        checkboxInput("show_hf_adv", "Toon fijnafstelling (cardiorenaal, stuwing, oedeem, albumine, F1)", value = FALSE),
        conditionalPanel(
          condition = "input.show_hf_adv || input.hf_preset != 'standard'",
          sliderInput("rf_pct", "Cardiorenaal GFR-behoud t.o.v. norm (%)", min = 25, max = 120, value = 100, step = 5),
          sliderInput("hep_pct", "Hepatische UGT1A9-functie bij stuwing (%)", min = 30, max = 100, value = 100, step = 5),
          sliderInput("edema_pct", "Overvulling / oedeem boven drooggewicht (%)", min = 0, max = 30, value = 0, step = 1),
          conditionalPanel(
            condition = "input.edema_pct > 0 && input.dose_mode == 'mgkg'",
            radioButtons(
              "dose_wt_basis", "Gewichtsbasis voor mg/kg-dosering",
              choices = c("Drooggewicht (aanbevolen)" = "dry", "Actueel oedemateus gewicht" = "actual"),
              selected = "dry"
            )
          ),
          radioButtons(
            "alb_mode", "Serumalbumine (eiwitbinding fu)",
            choices = c("Fysiologisch naar leeftijd (McNamara)" = "auto", "Handmatig (hypoalbuminemie / PLE)" = "manual"),
            selected = "auto"
          ),
          conditionalPanel(
            condition = "input.alb_mode == 'manual'",
            sliderInput("alb_val", "Serumalbumine (g/L; volw. ref = 44)", min = 15, max = 50, value = 30, step = 1)
          ),
          sliderInput("f1_pct", "Gemiddelde relatieve opname F1 (sonde/suspensie, %)", min = 30, max = 120, value = 100, step = 5)
        )
      ),
      selectInput("sex", "Geslacht", choices = c("Beide" = "both", "Man / jongen" = "male", "Vrouw / meisje" = "female")),
      sliderInput("n_subj", "Aantal virtuele individuen (N)", min = 50, max = 2000, value = 1000, step = 50),
      numericInput("seed", "Random seed", value = 42, min = 1),
      hr(),
      h4("Weergave"),
      radioButtons(
        "conc_scale", "Concentratie op PK-curve",
        choices = c("Totaal plasma C(t)" = "total", "Ongebonden vrij plasma Cu(t) + SGLT2-drempel" = "unbound"),
        selected = "total"
      ),
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
      p("Deze app simuleert de farmacokinetiek én semi-mechanistische farmacodynamiek van dapagliflozine in een virtuele populatie (", strong("copulapop"), ") van pasgeborene tot oudere, inclusief pathofysiologische aanpassingen en dosisvoorspelling voor ", strong("pediatrisch hartfalen"), ". Stel links dosis, leeftijdsvenster en eventueel een hartfalenprofiel in en klik op ", strong("Simuleer"), "."),
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
          downloadButton("dl_metrics", "Download individuele PK/PD-metrieken (CSV)")
        ),
        tabPanel(
          "Hartfalen: PK/PD & Dosisvoorspelling",
          h4("Dosisselectie & PK–PD analyse bij pediatrisch hartfalen"),
          uiOutput("hf_kpi_cards"),
          uiOutput("hf_alerts_ui"),
          div(
            class = "hf-subtabs",
            tabsetPanel(
              id = "hf_inner_tabs",
              type = "pills",
              tabPanel(
                "1. Therapeutisch venster & Referenties",
                div(
                  class = "hf-intro-box",
                  strong("Expositie-matching (ICH E11A) naar volwassen hartfalenstudies: "),
                  "In ", strong("DAPA-HF"), " en ", strong("DELIVER"), " is een vaste onderhoudsdosis van ",
                  strong("10 mg 1×/dag"), " effectief en veilig gebleken over een breed eGFR-spectrum. Bij gebrek aan pediatrische fase-III uitkomststudies richt dosisextrapolatie zich op ",
                  strong("blootstelling-matching"), " aan het ongebonden steady-state expositiebereik van deze volwassen referentiepopulaties.",
                  tags$br(), tags$br(),
                  strong("Waarom sturen op ongebonden blootstelling (AUCu,24h,ss) in plaats van totale plasmaconcentraties? "),
                  "Dapagliflozine is sterk gebonden aan serumalbumine (", code("91,4%"), "; ", code("fu,volw = 0,086"),
                  ") en gedraagt zich als een laag-extractiefarmacon (", code("ER < 0,05"), "). ",
                  "Bij fysiologische neonatale hypoalbuminemie, veneuze stuwing, cardiale cachexie of ", em("protein-losing enteropathy"), " (Fontan) stijgt ", code("fu"),
                  ": hierdoor neemt de totale klaring ", code("CL/F"), " proportioneel toe en daalt de totale plasmaconcentratie ", code("C(t)"),
                  ", terwijl de farmacologisch actieve ongebonden concentratie ", code("Cu(t)"), " (die glomerulair gefiltreerd wordt richting SGLT2) onveranderd blijft."
                ),
                plotOutput("hf_window_plot", height = "370px"),
                uiOutput("hf_window_cards_ui"),
                checkboxInput("show_custom_window", "Pas streefvenster en veiligheidsplafond handmatig aan (gevoeligheidsanalyse)", value = FALSE),
                conditionalPanel(
                  condition = "input.show_custom_window",
                  wellPanel(
                    fluidRow(
                      column(4, numericInput("hf_target_lo", "Ondergrens streefvenster (ng*h/mL; ref 20)", value = 20, min = 5, max = 60, step = 2)),
                      column(4, numericInput("hf_target_hi", "Bovengrens streefvenster (ng*h/mL; ref 75)", value = 75, min = 30, max = 150, step = 5)),
                      column(4, numericInput("hf_safety_max", "Veiligheidsplafond P95 (ng*h/mL; ref 90)", value = 90, min = 50, max = 200, step = 5))
                    ),
                    helpText("Aanpassingen werken alle PTA-curves, ontogenievergelijkingen en tabellen direct analytisch bij zonder herstart van de Monte Carlo-simulatie.")
                  )
                )
              ),
              tabPanel(
                "2. Dosisselectie (PTA) & Leeftijdstabel",
                div(
                  class = "hf-intro-box",
                  strong("Dosisoptimalisatie & gevoeligheid voor UGT1A9-ontogenie: "),
                  "Figuur A toont de ", strong("Probability of Target Attainment (PTA)"), " voor het gesimuleerde cohort: het percentage patiënten binnen het streefvenster (",
                  code("AUCu,24h,ss"), " tussen onder- en bovengrens) onder de randvoorwaarde dat ", code("P95"), " onder het veiligheidsplafond blijft. ",
                  "Figuur B kwantificeert de structurele modelonzekerheid in postnatale ", strong("UGT1A9-glucuronidering"), ": onder 6 maanden wijken snelle ", em("in-vivo"), " rijping (",
                  em("Farhan 2024"), ") en conservatieve PBPK-ontogenie (", em("Simcyp s1"), ") een factor 3 tot 4 uiteen, wat een voorzichtige startdosis met stapsgewijze optitratie onderbouwt."
                ),
                plotOutput("hf_pta_ontogeny_plot", height = "430px"),
                h5(strong("Expositie en dosisadvies per UGT1A9-ontogeniemodel voor het huidige cohort:")),
                tableOutput("hf_ontogeny_table"),
                hr(),
                h4("Leeftijdsgestratificeerd referentie-overzicht (0 mnd – 18 jaar)"),
                p("Onderstaande tabel evalueert standaardleeftijdsstrata onder het links ingestelde hartfalen- en formuleringsprofiel. ",
                  strong("Allometrische consequentie boven 1 jaar: "),
                  "doordat metabole UGT1A9-klaring schaalt met ", code("NFM^0,75"), ", is de gewichtsgenormaliseerde klaring bij peuters en kleuters hoger dan bij volwassenen en is ",
                  code("~0,18–0,20 mg/kg/dag"), " op euvolemisch drooggewicht nodig voor een volwassen ", code("10 mg"), "-equivalente blootstelling (", code("10 mg / 84 kg = 0,12 mg/kg"), ")."),
                tableOutput("hf_lifespan_table"),
                downloadButton("dl_hf_regimen", "Download leeftijdsgestratificeerde HF-referentietabel (CSV)")
              ),
              tabPanel(
                "3. PK–PD Relatie & Klinische Veiligheid",
                div(
                  class = "hf-explain-grid-3",
                  div(
                    class = "hf-explain-card",
                    div(class = "hf-explain-title", "1. Duale farmacodynamiek (luminaal vs. systemisch)"),
                    p("Dapagliflozine remt ", strong("luminaal SGLT2"), " in het S1/S2-segment na glomerulaire filtratie van ", code("Cu(t)"),
                      " (", code("IC50 = 1,1 nM"), "; tubulaire concentratiefactor ", code("2,25x"), " → ", code("Cu,50 = 0,20 ng/mL"),
                      "), naast ", strong("systemische / tubuloglomerulaire Emax-effecten"), " op intraglomerulaire druk, ketonmetabolisme en erytropoëse.")
                  ),
                  div(
                    class = "hf-explain-card hf-zone-high",
                    div(class = "hf-explain-title", "2. PK–PD dissociatie bij lage GFR (ceiling effect)"),
                    p("Bij onrijpe nefronen of cardiorenaal syndroom daalt de klaring (", code("AUCu,24h"), " stijgt), maar beperkt de lage ",
                      strong("gefiltreerde natrium- en glucosevracht"), " de osmotische natriurese terwijl SGLT2 al voor ", code("85–92%"), " verzadigd is. ",
                      strong("Opwaarts titreren om extra diurese af te dwingen is farmacologisch ineffectief."))
                  ),
                  div(
                    class = "hf-explain-card hf-zone-max",
                    div(class = "hf-explain-title", "3. Verschuiving therapeutische breedte (< 2 jaar)"),
                    p("Jonge kinderen hebben beperkte hepatische glycogeenreserves: normoglykemische glucosurie geeft bij vasten of intercurrente infectie risico op ",
                      strong("hypoglykemie en euglykemische ketoacidose"), " (", em("sick-day rules"), "). ",
                      "Bewaak daarnaast op prerenale AKI (naast lisdiuretica/RAS-remming) en renaal calorisch verlies bij ", em("failure to thrive"), ".")
                  )
                ),
                plotOutput("hf_pkpd_plot", height = "445px"),
                helpText("Toelichting bij figuur B: elk punt representeert een gesimuleerde patiënt, gekleurd naar eGFR. Patiënten met verlaagde eGFR (rood/oranje) vertonen hoge ongebonden blootstelling en maximale systemische Emax-respons (zwarte curves), maar een begrensde effectieve tubulaire natriuretische respons op de y-as door lage glomerulaire filtratie.")
              ),
              tabPanel(
                "4. Methodologische Onderbouwing & Bronnen",
                withMathJax(div(class = "md-doc", includeMarkdown("hf_pkpd_rationale.md")))
              )
            )
          )
        ),
        tabPanel(
          "Expositiemetrieken (PK & PD)",
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
          h4("Eerste 100 individuen uit de virtuele populatie (incl. drooggewicht, FFM, fu en effectieve eGFR)"),
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

  output$hf_preset_note <- renderUI({
    pr <- hf_presets()[[input$hf_preset]]
    if (is.null(pr) || is.null(pr$note)) return(NULL)
    helpText(pr$note)
  })

  # Synchroniseer hartfalen-fysiologie schuifjes wanneer een klinisch preset wordt gekozen
  observeEvent(input$hf_preset, ignoreInit = TRUE, {
    pr <- hf_presets()[[input$hf_preset]]
    if (!is.null(pr) && !is.null(pr$rf_factor)) {
      updateSliderInput(session, "rf_pct", value = round(100 * pr$rf_factor))
      updateSliderInput(session, "hep_pct", value = round(100 * pr$cl_hep_factor))
      updateSliderInput(session, "edema_pct", value = pr$edema_pct)
      updateRadioButtons(session, "alb_mode", selected = pr$alb_mode)
      updateSliderInput(session, "alb_val", value = pr$alb_val)
      updateSliderInput(session, "f1_pct", value = round(100 * pr$f1_rel))
    }
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

      rf_fac  <- if (length(vals$rf_pct) == 1L && !is.na(vals$rf_pct)) as.numeric(vals$rf_pct) / 100 else 1.0
      hep_fac <- if (length(vals$hep_pct) == 1L && !is.na(vals$hep_pct)) as.numeric(vals$hep_pct) / 100 else 1.0
      ed_pct  <- if (length(vals$edema_pct) == 1L && !is.na(vals$edema_pct)) as.numeric(vals$edema_pct) else 0.0
      wt_bas  <- if (ed_pct > 0 && length(vals$dose_wt_basis) == 1L && nzchar(vals$dose_wt_basis)) vals$dose_wt_basis else "dry"
      alb_ov  <- if (identical(vals$alb_mode, "manual") && length(vals$alb_val) == 1L && !is.na(vals$alb_val)) as.numeric(vals$alb_val) else NULL
      f1_rel  <- if (length(vals$f1_pct) == 1L && !is.na(vals$f1_pct)) as.numeric(vals$f1_pct) / 100 else 1.0

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

      res <- withProgress(message = "Monte Carlo PK/PD-simulatie ...", value = 0.70, {
        simulate_population_pk(
          pop,
          dose_mgkg         = d_mgkg,
          dose_mg           = d_mg,
          tau               = tau_val,
          n_doses           = nd_val,
          iov_f1_cv         = iov_val,
          ka_bsv            = (ka_m != "none"),
          ka_mode           = ka_m,
          mat_option        = mat_val,
          steady_state      = use_ss,
          seed              = as.integer((as.numeric(seed) + 1000) %% 2147483647),
          rf_factor         = rf_fac,
          cl_hep_factor     = hep_fac,
          edema_pct         = ed_pct,
          dose_weight_basis = wt_bas,
          alb_override      = alb_ov,
          f1_rel            = f1_rel
        )
      })

      if (!use_ss && nd_val > 1L && is.finite(res$settings$ss_attainment) && res$settings$ss_attainment < 0.90) {
        warns <- c(warns, sprintf(
          "Na %d doses (interval %g u) is mediaan %.0f%% van steady-state blootstelling bereikt (mediane t1/2,beta = %.1f u). Verhoog 'Aantal doses' of vink 'Directe analytische steady state' aan voor volledige steady state.",
          nd_val, tau_val, 100 * res$settings$ss_attainment, stats::median(res$metrics$t_half)
        ))
      }
      if (ed_pct > 0 && identical(wt_bas, "actual") && d_mode == "mgkg") {
        warns <- c(warns, sprintf(
          "Er wordt nu op actueel oedemateus gewicht (+%g%% overvulling) gedoseerd terwijl UGT1A9-klaring schaalt met het drooggewicht: dit verhoogt de blootstelling met +%g%%.",
          ed_pct, ed_pct
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
        sprintf("Simulatie voltooid: N = %d, dosis %s per gift, RF %d%%, seed %d",
                nrow(pop), dose_txt, round(100 * rf_fac), seed),
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

  build_curve_plot <- function(s, conc_scale, show_indiv, show_day1, log_y, obs = NULL, pub = FALSE) {
    st   <- s$res$settings
    is_u <- identical(conc_scale, "unbound")
    cs   <- if (is_u) s$res$cs_ss_u else s$res$cs_ss
    if (isTRUE(log_y)) {
      cs <- cs[cs$median > 0 & cs$q95 > 0, , drop = FALSE]
      cs$q05 <- ifelse(cs$q05 > 0, cs$q05, pmin(cs$median * 0.1, 1e-4))
    }

    dose_lbl_nl <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg", st$dose_mgkg)
    dose_lbl_en <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg per dose", st$dose_mgkg)
    regimen_nl  <- if (isTRUE(st$steady_state)) "analytische steady state" else if (st$n_doses == 1) "eenmalige dosis" else sprintf("dosis %d", st$n_doses)
    regimen_en  <- if (isTRUE(st$steady_state)) "analytical steady state" else if (st$n_doses == 1) "single dose" else sprintf("dose %d", st$n_doses)

    hf_tag_nl <- if (abs(st$rf_factor - 1) > 1e-3 || abs(st$cl_hep_factor - 1) > 1e-3 || !is.na(st$alb_override) || st$edema_pct > 0) {
      sprintf(" | HF: RF %.0f%%, Hep %.0f%%, fu med. %.3f", 100 * st$rf_factor, 100 * st$cl_hep_factor, stats::median(s$res$metrics$FU))
    } else ""
    hf_tag_en <- if (nzchar(hf_tag_nl)) {
      sprintf(" | HF: RF %.0f%%, Hep %.0f%%, med fu %.3f", 100 * st$rf_factor, 100 * st$cl_hep_factor, stats::median(s$res$metrics$FU))
    } else ""

    ylab_nl <- if (is_u) "Ongebonden dapagliflozine Cu(t) in plasma (ng/mL)" else "Totale dapagliflozine C(t) in plasma (ng/mL)"
    ylab_en <- if (is_u) "Unbound plasma dapagliflozin Cu(t) (ng/mL)" else "Total plasma dapagliflozin C(t) (ng/mL)"

    p <- ggplot(cs, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.25) +
      geom_line(aes(y = q95), color = "#2c7fb8", linetype = "dashed") +
      geom_line(aes(y = q05), color = "#2c7fb8", linetype = "dashed") +
      geom_line(color = "#08306b", linewidth = 1.2) +
      do.call(labs, if (pub) list(
        x = "Time after dose (h)",
        y = ylab_en,
        title = sprintf("%sconcentration-time profile (%s, interval %g h)",
                        if (is_u) "Unbound " else "", regimen_en, st$tau),
        subtitle = sprintf(
          "Median + 90%% PI (P5-P95), N = %d (%s), dose %s, age %.2f-%.2f y, UGT1A9: %s%s",
          st$n, s$dataset, dose_lbl_en, min(s$pop$AGE), max(s$pop$AGE), st$mat_label, hf_tag_en
        )
      ) else list(
        x = "Tijd na dosis (uur)",
        y = ylab_nl,
        title = sprintf("%sconcentratie-tijdcurve (%s, interval %g u)",
                        if (is_u) "Ongebonden " else "Totale ", regimen_nl, st$tau),
        subtitle = sprintf(
          "Mediaan + 90%% PI (P5-P95), N = %d (%s), dosis %s, leeftijd %.2f-%.2f j, maturatie: %s%s",
          st$n, s$dataset, dose_lbl_nl, min(s$pop$AGE), max(s$pop$AGE), st$mat_label, hf_tag_nl
        )
      ))

    if (is_u) {
      # Toon effectieve Cu,50 (0.20 ng/mL = 50% luminale SGLT2-remming) en Cu,90 (1.80 ng/mL = 90% remming)
      p <- p +
        geom_hline(yintercept = 0.20, color = "#1a7f37", linetype = "dotted", linewidth = 0.8) +
        geom_hline(yintercept = 1.80, color = "#1a7f37", linetype = "dotdash", linewidth = 0.7) +
        annotate("text", x = st$tau * 0.98, y = 0.22, hjust = 1, vjust = -0.3,
                 label = if (pub) "50% luminal SGLT2 inh. (Cu,50 = 0.20 ng/mL)" else "50% luminale SGLT2-remming (Cu,50 = 0.20 ng/mL)",
                 color = "#1a7f37", size = 3.2, fontface = "bold") +
        annotate("text", x = st$tau * 0.98, y = 1.85, hjust = 1, vjust = -0.3,
                 label = if (pub) "90% luminal SGLT2 inh. (Cu,90 = 1.80 ng/mL)" else "90% luminale SGLT2-remming (Cu,90 = 1.80 ng/mL)",
                 color = "#1a7f37", size = 3.2)
    }

    if (isTRUE(show_indiv)) {
      idx <- s$indiv_idx
      mat_c <- if (is_u) s$res$curve_ss_u else s$res$curve_ss
      indiv <- data.frame(
        t  = rep(s$res$times, each = length(idx)),
        C  = as.vector(mat_c[idx, , drop = FALSE]),
        ID = rep(idx, times = length(s$res$times))
      )
      if (isTRUE(log_y)) indiv <- indiv[indiv$C > 0, , drop = FALSE]
      if (nrow(indiv) > 0) {
        p <- p + geom_line(data = indiv, aes(x = t, y = C, group = ID),
                           color = "grey55", alpha = 0.35, linewidth = 0.4)
      }
    }
    if (isTRUE(show_day1)) {
      cs1 <- if (is_u) s$res$cs_day1_u else s$res$cs_day1
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
      med_fu <- stats::median(s$res$metrics$FU, na.rm = TRUE)
      obs_plot <- obs
      if (is_u) obs_plot$conc <- obs_plot$conc * med_fu
      if (isTRUE(log_y)) obs_plot <- obs_plot[obs_plot$conc > 0, , drop = FALSE]
      if (nrow(obs_plot) > 0) {
        obs_lbl <- if (is_u) {
          sprintf(if (pub) "orange = obs x median fu %.3f (n = %d)" else "oranje = gemeten data x mediane fu %.3f (n = %d)", med_fu, nrow(obs_plot))
        } else {
          sprintf(if (pub) "orange points = observed data (n = %d)" else "oranje punten = gemeten data (n = %d)", nrow(obs_plot))
        }
        p <- p + geom_point(data = obs_plot, aes(x = time, y = conc), inherit.aes = FALSE,
                            color = "#D55E00", size = 2.4, alpha = 0.9) +
          annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6,
                   label = obs_lbl, color = "#D55E00", size = 3.5)
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
    build_curve_plot(s, input$conc_scale, input$show_indiv, input$show_day1, input$log_y, obs = obs, pub = FALSE)
  })

  output$dl_fig <- downloadHandler(
    filename = function() paste0("dapagliflozine_pk_curve.", input$fig_format),
    content = function(file) {
      s <- sim_res(); req(s)
      obs <- if (isTRUE(input$show_obs)) obs_df() else NULL
      p <- build_curve_plot(s, input$conc_scale, input$show_indiv, input$show_day1, input$log_y, obs = obs, pub = TRUE)
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

  # ---------- Tab: Hartfalen PK/PD & Dosisvoorspelling ----------
  hf_targets_reactive <- reactive({
    lo <- suppressWarnings(as.numeric(input$hf_target_lo)); if (is.na(lo) || lo <= 0) lo <- 20
    hi <- suppressWarnings(as.numeric(input$hf_target_hi)); if (is.na(hi) || hi <= lo) hi <- max(lo + 10, 75)
    sf <- suppressWarnings(as.numeric(input$hf_safety_max)); if (is.na(sf) || sf < hi) sf <- max(hi + 10, 90)
    hf_reference_targets(target_lo = lo, target_hi = hi, safety_max = sf)
  })

  hf_pta_data <- reactive({
    s <- sim_res(); req(s)
    tg <- hf_targets_reactive()
    compute_pta_scan(s$res, target_lo = tg$target_lo, target_hi = tg$target_hi, safety_max = tg$safety_max)
  })

  hf_ontogeny_data <- reactive({
    s <- sim_res(); req(s)
    tg <- hf_targets_reactive()
    compare_ontogeny_dosing(s$pop, s$res$settings,
                            target_lo = tg$target_lo, target_hi = tg$target_hi, safety_max = tg$safety_max)
  })

  hf_lifespan_df <- reactive({
    s <- sim_res(); req(s)
    st <- s$res$settings
    tg <- hf_targets_reactive()
    generate_pediatric_hf_regimen_table(
      rf_factor     = st$rf_factor,
      cl_hep_factor = st$cl_hep_factor,
      edema_pct     = st$edema_pct,
      alb_override  = st$alb_override,
      f1_rel        = st$f1_rel,
      tau           = st$tau,
      safety_max    = tg$safety_max
    )
  })

  output$hf_kpi_cards <- renderUI({
    s   <- sim_res(); req(s)
    pta <- hf_pta_data(); req(pta)
    ont <- hf_ontogeny_data(); req(ont)
    tg  <- hf_targets_reactive()
    m   <- s$res$metrics
    st  <- s$res$settings

    med_wt     <- ont$med_wt_dry
    start_mg   <- min(5.0, ont$start_rec_mgkg * med_wt)
    targ_lo_mg <- min(10.0, ont$target_min_mgkg * med_wt)
    targ_hi_mg <- min(10.0, ont$target_max_mgkg * med_wt)

    targ_val_txt <- if (abs(ont$target_max_mgkg - ont$target_min_mgkg) / max(ont$target_min_mgkg, 1e-4) < 0.15) {
      sprintf("%.2f mg/kg", ont$target_max_mgkg)
    } else {
      sprintf("%.2f – %.2f mg/kg", ont$target_min_mgkg, ont$target_max_mgkg)
    }
    targ_sub_txt <- if (abs(targ_hi_mg - targ_lo_mg) < 0.08) {
      sprintf("~%.2f mg/gift bij mediaan drooggewicht (%.1f kg) — equivalent aan volwassen 10 mg DAPA-HF expositie", targ_hi_mg, med_wt)
    } else {
      sprintf("~%.2f – %.2f mg/gift (spreiding tussen Simcyp s1 en Farhan 2024 UGT1A9-ontogenie; plafond 10 mg)", targ_lo_mg, targ_hi_mg)
    }

    aucu_ss_med <- stats::median(m$FU * m$AUC0_inf * (24 / st$tau), na.rm = TRUE)
    cur_dose_lbl <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg", st$dose_mgkg)

    if (pta$cur_above_safety > 10) {
      status_cls <- "hf-kpi-value hf-status-danger"
      status_txt <- sprintf("%.0f%% > veiligheidsplafond", pta$cur_above_safety)
    } else if (pta$cur_below > 50) {
      status_cls <- "hf-kpi-value hf-status-warn"
      status_txt <- sprintf("Subtherapeutisch (%.0f%% < ondergrens)", pta$cur_below)
    } else {
      status_cls <- "hf-kpi-value hf-status-ok"
      status_txt <- sprintf("PTA %.0f%% in streefvenster", pta$cur_pta)
    }

    div(
      class = "hf-kpi-grid",
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label", "1. Conservatieve startdosis (~5 mg-eq)"),
        div(class = "hf-kpi-value", sprintf("%.2f mg/kg", ont$start_rec_mgkg)),
        div(class = "hf-kpi-sub",
            sprintf("~%.2f mg/gift — richt op AUCu,24h,ss ≈ 20 ng*h/mL (volwassen 5 mg-equivalent) onder het meest conservatieve UGT1A9-ontogeniemodel.", start_mg))
      ),
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label", "2. Onderhouds-streefdosis (~10 mg DAPA-HF-eq)"),
        div(class = "hf-kpi-value", targ_val_txt),
        div(class = "hf-kpi-sub", targ_sub_txt)
      ),
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label", sprintf("3. Expositiebeoordeling huidige dosis (%s)", cur_dose_lbl)),
        div(class = status_cls, status_txt),
        div(class = "hf-kpi-sub",
            sprintf("Mediane AUCu,24h,ss: %.1f ng*h/mL (streefvenster %g–%g | > %g ng*h/mL: %.0f%%).",
                    aucu_ss_med, tg$target_lo, tg$target_hi, tg$safety_max, pta$cur_above_safety))
      )
    )
  })

  output$hf_alerts_ui <- renderUI({
    ont <- hf_ontogeny_data(); req(ont)
    if (length(ont$alerts) == 0L) return(NULL)
    div(
      class = "hf-alert-panel",
      strong("Klinisch-farmacologische aandachtspunten voor het geselecteerde cohort:"),
      tags$ul(
        lapply(ont$alerts, function(a) {
          tags$li(strong(a$lead), " ", a$body)
        })
      )
    )
  })

  output$hf_window_plot <- renderPlot({
    s   <- sim_res(); req(s)
    ont <- hf_ontogeny_data(); req(ont)
    tg  <- hf_targets_reactive()
    build_window_explainer_plot(s, tg, ont)
  })

  output$hf_window_cards_ui <- renderUI({
    tg <- hf_targets_reactive()
    div(
      class = "hf-explain-grid-4",
      div(
        class = "hf-explain-card hf-zone-low",
        div(class = "hf-explain-title", sprintf("Ondergrens streefvenster: %g ng*h/mL", tg$target_lo)),
        p(strong("Farmacometrisch anker: "), "Mediane blootstelling bij ", strong("volwassen 5 mg 1×/dag"),
          " (", code("19,2 ng*h/mL"), " bij eGFR 90) en nabij het ", code("P5"), " van volwassenen op 10 mg (", code("23,7 ng*h/mL"), "; Penland 2023 Tabel S4)."),
        p(strong("PD-correlatie: "), "Geeft reeds ~79% gemiddelde 24-uurs SGLT2-bezetting en ~80% van het maximale TGF-effect; rationeel startanker bij initiatie.")
      ),
      div(
        class = "hf-explain-card hf-zone-ok",
        div(class = "hf-explain-title", "Streefmidden onderhoud: 38 – 48 ng*h/mL"),
        p(strong("Farmacometrisch anker: "), "Mediane blootstelling bij ", strong("volwassen 10 mg 1×/dag"),
          ": ", code("38,5 ng*h/mL"), " bij normale nierfunctie (eGFR 90) tot ", code("48,2 ng*h/mL"), " voor het gemiddelde ", strong("DAPA-HF"), " profiel (66 j, eGFR 66)."),
        p(strong("Klinische interpretatie: "), "Primair expositiedoel bij optitratie naar de pediatrische onderhoudsdosering (optimalisatie-anker 42 ng*h/mL).")
      ),
      div(
        class = "hf-explain-card hf-zone-ok",
        div(class = "hf-explain-title", sprintf("Bovengrens streefvenster: %g ng*h/mL", tg$target_hi)),
        p(strong("Farmacometrisch anker: "), "Komt overeen met het ", strong("95e percentiel (74,5 ng*h/mL)"),
          " bij volwassenen op 10 mg met milde tot matige nierfunctiestoornis (eGFR 60–89; ", code("866 × 0,086"), ", Penland 2023 Tabel S4)."),
        p(strong("Klinische interpretatie: "), sprintf("Exposities van %g–%g ng*h/mL vallen binnen de reguliere inter-individuele spreiding van de geregistreerde 10 mg HFrEF-dosering.", tg$target_lo, tg$target_hi))
      ),
      div(
        class = "hf-explain-card hf-zone-max",
        div(class = "hf-explain-title", sprintf("Veiligheidsplafond (P95): %g ng*h/mL", tg$safety_max)),
        p(strong("Farmacometrisch anker: "), "Gelijk aan de ", strong("mediane vrije blootstelling (87,7 ng*h/mL)"),
          " in ernstige chronische nierschade (eGFR 15–29) op 10 mg/dag (", code("1020 × 0,086"), "; Penland 2023 / DAPA-CKD)."),
        p(strong("Klinische interpretatie: "), "Bovengrens voor P95 bij dosisselectie; boven dit niveau ontbreken langetermijn-veiligheidsdata zonder dosisreductie.")
      )
    )
  })

  output$hf_pta_ontogeny_plot <- renderPlot({
    s   <- sim_res(); req(s)
    pta <- hf_pta_data(); req(pta)
    ont <- hf_ontogeny_data(); req(ont)
    tg  <- hf_targets_reactive()
    st  <- s$res$settings

    gdf <- pta$grid_df
    dose_unit_lbl <- if (st$dose_mode == "mg") "Vaste keerdosis (mg)" else "Gewichtsgedoseerde keerdosis (mg/kg)"
    mat_short_map <- c(
      farhan2024 = "Farhan 2024, in-vivo UGT1A9",
      badee2019  = "Badee 2019, PBPK s2",
      simcyp_s1  = "Simcyp s1, conservatief"
    )
    mat_short_lbl <- mat_short_map[[st$mat_option]]

    pta_long <- rbind(
      data.frame(dose = gdf$dose, pct = gdf$pct_in_target,
                 cat = sprintf("Binnen streefvenster (%g–%g ng*h/mL)", tg$target_lo, tg$target_hi)),
      data.frame(dose = gdf$dose, pct = gdf$pct_below,
                 cat = sprintf("Subtherapeutisch (< %g ng*h/mL)", tg$target_lo)),
      data.frame(dose = gdf$dose, pct = gdf$pct_above_safety,
                 cat = sprintf("> Veiligheidsplafond (> %g ng*h/mL)", tg$safety_max))
    )
    lvls <- unique(pta_long$cat)
    pta_long$cat <- factor(pta_long$cat, levels = lvls)
    cols_pta <- stats::setNames(c("#1a7f37", "#2c7fb8", "#b2182b"), lvls)
    ltys_pta <- stats::setNames(c("solid", "dotted", "dashed"), lvls)

    p_left <- ggplot(pta_long, aes(x = dose, y = pct, color = cat, linetype = cat)) +
      geom_line(linewidth = 1.15) +
      geom_vline(xintercept = pta$cur_dose, color = "#08306b", linewidth = 0.8) +
      geom_vline(xintercept = pta$dose_opt_pta, color = "#1a7f37", linetype = "dotdash", linewidth = 0.8) +
      annotate("text", x = pta$cur_dose, y = 97, hjust = -0.08, size = 3.2, fontface = "bold", color = "#08306b",
               label = sprintf("Huidig: %.3g", pta$cur_dose)) +
      scale_color_manual(values = cols_pta) +
      scale_linetype_manual(values = ltys_pta) +
      coord_cartesian(ylim = c(0, 103)) +
      labs(
        title = sprintf("A. Dosis–PTA curve (%s)", mat_short_lbl),
        subtitle = sprintf("PTA-optimum = %.3g %s (PTA %.0f%% bij P95 <= %g; groene streeplijn)",
                           pta$dose_opt_pta, if (st$dose_mode == "mg") "mg" else "mg/kg", pta$max_pta_pct, tg$safety_max),
        x = dose_unit_lbl,
        y = "Proportie virtueel cohort (%)",
        color = NULL, linetype = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title.position = "plot",
            plot.title = element_text(face = "bold", color = "#08306b", size = 12.2),
            plot.subtitle = element_text(size = 10.2, color = "#3d4f63"),
            plot.margin = margin(6, 12, 6, 6),
            legend.position = "bottom",
            legend.direction = "vertical")

    ddf <- ont$dist_df
    y_cap <- max(tg$safety_max * 1.25, stats::quantile(ddf$AUCu_24h, 0.97, na.rm = TRUE))
    p_right <- ggplot(ddf, aes(x = ontogenie, y = pmin(AUCu_24h, y_cap), fill = ontogenie)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = tg$target_lo, ymax = tg$target_hi,
               fill = "#e6f4ea", alpha = 0.75) +
      geom_hline(yintercept = tg$adult_10mg_hf, color = "#1a7f37", linetype = "dotted", linewidth = 0.8) +
      geom_hline(yintercept = tg$safety_max, color = "#b2182b", linetype = "dashed", linewidth = 0.85) +
      geom_boxplot(width = 0.46, outlier.alpha = 0.20, alpha = 0.85, color = "#1f2d3d") +
      scale_fill_manual(values = c("#9ecae1", "#fdd0a2", "#c7e9c0")) +
      coord_cartesian(ylim = c(0, y_cap)) +
      annotate("text", x = 3.45, y = tg$safety_max, hjust = 1, vjust = -0.4, size = 3.1,
               color = "#b2182b", fontface = "bold", label = sprintf("Veiligheidsplafond (%g)", tg$safety_max)) +
      annotate("text", x = 3.45, y = tg$adult_10mg_hf, hjust = 1, vjust = -0.4, size = 3.0,
               color = "#1a7f37", label = "DAPA-HF 10 mg mediaan (48,2)") +
      labs(
        title = "B. Gevoeligheid voor UGT1A9-ontogenie",
        subtitle = "Steady-state AUCu,24h bij huidige dosis (groene zone = streefvenster)",
        x = NULL,
        y = "Ongebonden blootstelling AUCu,24h,ss (ng*h/mL)"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title.position = "plot",
            plot.title = element_text(face = "bold", color = "#08306b", size = 12.2),
            plot.subtitle = element_text(size = 10.2, color = "#3d4f63"),
            plot.margin = margin(6, 12, 6, 8),
            legend.position = "none")

    side_by_side_plots(p_left, p_right)
  })

  output$hf_ontogeny_table <- renderTable({
    ont <- hf_ontogeny_data(); req(ont)
    ont$summary_df
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$hf_pkpd_plot <- renderPlot({
    s  <- sim_res(); req(s)
    st <- s$res$settings
    m  <- s$res$metrics
    tg <- hf_targets_reactive()

    cs_occ  <- s$res$cs_sglt2
    ref_occ <- adult_hf_reference_curves(tau = st$tau, out_by = 0.1)

    p_occ <- ggplot(cs_occ, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.22) +
      geom_line(color = "#08306b", linewidth = 1.2) +
      geom_line(data = ref_occ, aes(x = t, y = SGLT2_inh, linetype = ref, color = ref),
                inherit.aes = FALSE, linewidth = 0.95) +
      scale_color_manual(values = c(
        "Volwassen DAPA-HF ref: 10 mg 1x/dag" = "#b2182b",
        "Volwassen startdosis ref: 5 mg 1x/dag" = "#d95f02"
      )) +
      scale_linetype_manual(values = c(
        "Volwassen DAPA-HF ref: 10 mg 1x/dag" = "dashed",
        "Volwassen startdosis ref: 5 mg 1x/dag" = "dotted"
      )) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = "A. Proximaal tubulaire SGLT2-remming over het interval",
        subtitle = "Blauw = cohort (mediaan + 90% PI); referenties = volwassen 10 & 5 mg QD",
        x = "Tijd na toediening (uur)",
        y = "Luminale SGLT2-receptorbezetting (%)",
        color = NULL, linetype = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title.position = "plot",
            plot.title = element_text(face = "bold", color = "#08306b", size = 12.2),
            plot.subtitle = element_text(size = 10.2, color = "#3d4f63"),
            plot.margin = margin(6, 12, 6, 6),
            legend.position = "bottom",
            legend.direction = "vertical")

    # Rechts: blootstelling-respons curves vs GFR-gelimiteerde renale PD per individu
    x_max <- max(110, min(220, stats::quantile(m$AUCu_24h, 0.96, na.rm = TRUE) * 1.15))
    x_seq <- seq(0.5, x_max, length.out = 200)
    th    <- nonmem_theta()
    er_curves <- rbind(
      data.frame(aucu = x_seq, resp = 100 * x_seq / (th$eauc50_u_tgf + x_seq),
                 as_type = "Glomerulaire hemodynamiek / TGF Emax (EAUCu,50 = 5,2 ng*h/mL)"),
      data.frame(aucu = x_seq, resp = 100 * x_seq / (th$eauc50_u_hct + x_seq),
                 as_type = "Systemisch / hematocriet Emax (EAUCu,50 = 20,7 ng*h/mL)")
    )

    pts <- m[m$AUCu_24h <= x_max * 1.05, , drop = FALSE]
    if (nrow(pts) > 600L) {
      set.seed(42)
      pts <- pts[sample(nrow(pts), 600L), , drop = FALSE]
    }

    p_er <- ggplot() +
      annotate("rect", xmin = tg$target_lo, xmax = tg$target_hi, ymin = -Inf, ymax = Inf,
               fill = "#e6f4ea", alpha = 0.65) +
      geom_point(data = pts, aes(x = AUCu_24h, y = Renal_PD_eff, color = EGFR),
                 alpha = 0.65, size = 1.8) +
      geom_line(data = er_curves, aes(x = aucu, y = resp, linetype = as_type),
                color = "#111111", linewidth = 1.0) +
      scale_color_gradientn(
        colors = c("#b2182b", "#ef8a62", "#67a9cf", "#08306b"),
        limits = c(15, 130),
        oob = scales::squish,
        name = "eGFR\n(mL/min/1.73m2)"
      ) +
      coord_cartesian(xlim = c(0, x_max), ylim = c(0, 115)) +
      labs(
        title = "B. Systemische Emax vs. GFR-gelimiteerde tubulaire PD",
        subtitle = "Punten = tubulaire solute-respons (kleur = eGFR); lijnen = systemische Emax",
        x = "Ongebonden blootstelling AUCu,24h (ng*h/mL)",
        y = "Relatieve PD-respons (%)",
        linetype = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title.position = "plot",
            plot.title = element_text(face = "bold", color = "#08306b", size = 12.2),
            plot.subtitle = element_text(size = 10.2, color = "#3d4f63"),
            plot.margin = margin(6, 12, 6, 8),
            legend.position = "bottom",
            legend.box = "vertical")

    side_by_side_plots(p_occ, p_er)
  })

  output$hf_lifespan_table <- renderTable({
    hf_lifespan_df()
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$dl_hf_regimen <- downloadHandler(
    filename = "dapagliflozine_pediatrisch_hartfalen_doseringstabel.csv",
    content = function(file) {
      df <- hf_lifespan_df(); req(df)
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  # ---------- Expositiemetrieken ----------
  output$ss_table <- renderTable({
    s <- sim_res(); req(s)
    summ <- summarize_metrics(s$res$metrics, tau = s$res$settings$tau, include_pd = TRUE)
    summ[] <- lapply(summ, function(x) if (is.numeric(x)) round(x, 2) else x)
    summ
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$metric_plots <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    met_long <- data.frame(
      metriek = factor(
        rep(c("Cmax totaal (ng/mL)", "AUCtau totaal (ng*h/mL)", "AUCu ongebonden (ng*h/mL)", "SGLT2-remming gem. (%)"),
            each = nrow(m)),
        levels = c("Cmax totaal (ng/mL)", "AUCtau totaal (ng*h/mL)", "AUCu ongebonden (ng*h/mL)", "SGLT2-remming gem. (%)")
      ),
      waarde = c(m$Cmax, m$AUCtau, m$AUCu, m$SGLT2_inh_avg)
    )
    ggplot(met_long, aes(x = metriek, y = waarde)) +
      geom_boxplot(fill = "#9ecae1", outlier.alpha = 0.4) +
      facet_wrap(~metriek, scales = "free", ncol = 4) +
      labs(x = NULL, y = NULL,
           title = "Verdeling van PK- en PD-metrieken (totaal, ongebonden & luminale SGLT2-remming)") +
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
    ggplot(m, aes(x = AGE, y = AUCu_24h, color = Geslacht)) +
      geom_hline(yintercept = c(20, 75), linetype = "dashed", color = "#1a7f37", linewidth = 0.6) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(x = "Leeftijd (jaren)", y = "Ongebonden AUCu,24h (ng*h/mL)",
           title = "Ongebonden blootstelling (AUCu,24h) vs leeftijd (groene lijnen = 20-75 ng*h/mL referenties)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$exposure_covariates <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    lbls <- sex_labels(m$AGE)
    m$Geslacht <- factor(m$SEXISMALE, levels = c(0, 1), labels = lbls)
    long <- rbind(
      data.frame(metriek = "AUCu ongebonden (ng*h/mL)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$AUCu, Geslacht = m$Geslacht),
      data.frame(metriek = "AUCu ongebonden (ng*h/mL)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$AUCu, Geslacht = m$Geslacht),
      data.frame(metriek = "GFR-gekoppelde renale PD (%)", xvariabele = "eGFR (mL/min/1.73m2)",
                 x = m$EGFR, y = m$Renal_PD_eff, Geslacht = m$Geslacht),
      data.frame(metriek = "GFR-gekoppelde renale PD (%)", xvariabele = "Leeftijd (jaren)",
                 x = m$AGE, y = m$Renal_PD_eff, Geslacht = m$Geslacht)
    )
    ggplot(long, aes(x = x, y = y, color = Geslacht)) +
      geom_point(alpha = 0.4, size = 1.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1) +
      facet_grid(metriek ~ xvariabele, scales = "free") +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(x = NULL, y = NULL,
           title = "Ongebonden PK (AUCu) én renale PD vs nierfunctie (eGFR) en leeftijd",
           subtitle = "Toont direct de tegengestelde richting bij lage GFR: AUCu stijgt, terwijl tubulaire osmotische PD daalt") +
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
    st   <- s$res$settings
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
      sprintf("Actueel gewicht (BWT): mediaan %.1f kg | Drooggewicht: mediaan %.1f kg | Vetvrije massa (FFM): mediaan %.1f kg",
              stats::median(pop$WEIGHT), stats::median(m$BWT_DRY), stats::median(m$FFM)),
      sprintf("Serumalbumine: mediaan %.1f g/L | Ongebonden fractie fu: mediaan %.3f (range %.3f-%.3f)",
              stats::median(m$ALB), stats::median(m$FU), min(m$FU), max(m$FU)),
      sprintf("Effectieve eGFR: mediaan %.0f mL/min/1.73m2 (range %.0f-%.0f) | Relatieve nierfunctie RF: mediaan %.2f",
              stats::median(m$EGFR), min(m$EGFR), max(m$EGFR), stats::median(m$RF)),
      sprintf("Hartfalen- & sondeparameters: cardiorenaal RF x%.2f | hepatische UGT1A9 x%.2f | oedeem +%g%% | F1 %.0f%%",
              st$rf_factor, st$cl_hep_factor, st$edema_pct, 100 * st$f1_rel),
      sprintf("eGFR-bron: %s", st$egfr_source),
      sprintf("eGFR-eenheid: mL/min/1.73m2 — %s", st$egfr_unit_note),
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
      data.frame(variabele = "Drooggewicht (kg)", waarde = m$BWT_DRY),
      data.frame(variabele = "Vetvrije massa FFM (kg)", waarde = m$FFM),
      data.frame(variabele = "Effectieve eGFR (mL/min/1.73m2)", waarde = m$EGFR)
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
    out[["Drooggewicht (kg)"]] <- round(m$BWT_DRY[idx], 2)
    out[["FFM (kg)"]] <- round(m$FFM[idx], 2)
    out[["Alb (g/L)"]] <- round(m$ALB[idx], 1)
    out[["fu"]] <- round(m$FU[idx], 3)
    out[["EGFR (mL/min/1.73m2)"]] <- round(m$EGFR[idx], 1)
    head(out, 100)
  })

  # ---------- Downloads ----------
  output$dl_curve <- downloadHandler(
    filename = "dapagliflozine_pk_curve.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      df_out <- cbind(
        s$res$cs_ss,
        Cu_q05    = s$res$cs_ss_u$q05,
        Cu_median = s$res$cs_ss_u$median,
        Cu_q95    = s$res$cs_ss_u$q95,
        SGLT2_inh_median = s$res$cs_sglt2$median
      )
      utils::write.csv(df_out, file, row.names = FALSE)
    }
  )

  output$dl_metrics <- downloadHandler(
    filename = "dapagliflozine_pkpd_metrics.csv",
    content = function(file) {
      s <- sim_res(); req(s)
      utils::write.csv(s$res$metrics, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
