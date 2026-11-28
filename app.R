required_pkgs <- c("shiny", "ggplot2", "markdown",
                   "VineCopula", "data.table", "jsonlite", "Matrix", "MASS",
                   "copulapop")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages / Ontbrekende packages: ", paste(missing_pkgs, collapse = ", "), "\n",
       "Install with / Installeer met:\n",
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

get_mat_choices <- function(lang = "en") {
  if (identical(lang, "nl")) {
    c(
      "Farhan et al. 2024 — jcph.2484 (in-vivo UGT1A9, standaard)" = "farhan2024",
      "Badee et al. 2019 via Jo et al. 2021 — psp4.12577 (PBPK s2)" = "badee2019",
      "Simcyp standaard via Jo et al. 2021 — psp4.12577 (PBPK s1)" = "simcyp_s1"
    )
  } else {
    c(
      "Farhan et al. 2024 — jcph.2484 (in-vivo UGT1A9, default)" = "farhan2024",
      "Badee et al. 2019 via Jo et al. 2021 — psp4.12577 (PBPK s2)" = "badee2019",
      "Simcyp default via Jo et al. 2021 — psp4.12577 (PBPK s1)" = "simcyp_s1"
    )
  }
}

get_ka_choices <- function(lang = "en") {
  if (identical(lang, "nl")) {
    c(
      "Gekalibreerd klinisch (ω² = 0.36; standaard)" = "calibrated",
      "Penland 2023 bronmodel (ω² = 8.52 FIX)" = "penland_raw",
      "Geen Ka-variabiliteit (ω² = 0)" = "none"
    )
  } else {
    c(
      "Clinically calibrated (ω² = 0.36; default)" = "calibrated",
      "Penland 2023 source model (ω² = 8.52 FIX)" = "penland_raw",
      "No Ka variability (ω² = 0)" = "none"
    )
  }
}

get_datasets <- function(lang = "en") {
  if (identical(lang, "nl")) {
    c(
      "RUMC (Nederland, ziekenhuispopulatie, 0–100 j, met CREAT)" = "RUMC",
      "NHANES (USA, algemene bevolking, 0–80 j)" = "NHANES",
      "HSE (Engeland, algemene bevolking, 2–90 j)" = "HSE"
    )
  } else {
    c(
      "RUMC (Netherlands, hospital population, 0–100 y, with CREAT)" = "RUMC",
      "NHANES (USA, general population, 0–80 y)" = "NHANES",
      "HSE (England, general population, 2–90 y)" = "HSE"
    )
  }
}

get_hf_preset_choices <- function(lang = "en") {
  pr <- hf_presets(lang)
  stats::setNames(names(pr), vapply(pr, `[[`, character(1), "label"))
}

get_sex_choices <- function(lang = "en") {
  if (identical(lang, "nl")) {
    c("Beide" = "both", "Man / jongen" = "male", "Vrouw / meisje" = "female")
  } else {
    c("Both" = "both", "Male / boy" = "male", "Female / girl" = "female")
  }
}

get_fig_format_choices <- function(lang = "en") {
  c("PNG (raster)" = "png", "PDF (vector)" = "pdf", "SVG (vector)" = "svg", "TIFF (raster)" = "tiff")
}

i18n_span <- function(en, nl) {
  tags$span(class = "i18n-txt", `data-en` = en, `data-nl` = nl, en)
}

i18n_div <- function(en_ui, nl_ui) {
  tagList(
    div(class = "i18n-en-block", en_ui),
    div(class = "i18n-nl-block", style = "display: none;", nl_ui)
  )
}

side_by_side_plots <- function(p_left, p_right) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_left,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

# Gebruik MathJax CommonHTML i.p.v. Safari's gebrekkige NativeMML-fallback (voorkomt ontbrekende symbolen)
options(shiny.mathjax.config = "config=TeX-AMS_CHTML")

# Render Markdown waarbij zowel display-formules ($$...$$, \[...\]) als inline-formules (\(...\))
# eerst tegen CommonMark-escaping worden beschermd zodat \(, \), _ en \\ intact aan MathJax worden doorgegeven.
include_math_md <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  mask_pattern <- function(s, pat, store_env, prefix) {
    m <- gregexpr(pat, s, perl = TRUE)
    hits <- regmatches(s, m)[[1]]
    if (length(hits) == 0L) return(s)
    start_idx <- length(store_env$items)
    tokens <- sprintf("%s%04dEND", prefix, start_idx + seq_along(hits))
    store_env$items <- c(store_env$items, hits)
    regmatches(s, m) <- list(tokens)
    s
  }

  c_env <- new.env(parent = emptyenv()); c_env$items <- character(0)
  txt <- mask_pattern(txt, "(?s)```.*?```", c_env, "MDCODETOKEN")
  txt <- mask_pattern(txt, "`[^`\\n]+`",   c_env, "MDCODETOKEN")

  m_env <- new.env(parent = emptyenv()); m_env$items <- character(0)
  txt <- mask_pattern(txt, "(?s)\\$\\$.+?\\$\\$",     m_env, "MDMATHTOKEN")
  txt <- mask_pattern(txt, "(?s)\\\\\\[.+?\\\\\\]",   m_env, "MDMATHTOKEN")
  txt <- mask_pattern(txt, "(?s)\\\\\\(.+?\\\\\\)",   m_env, "MDMATHTOKEN")

  for (i in seq_along(c_env$items)) {
    txt <- gsub(sprintf("MDCODETOKEN%04dEND", i), c_env$items[[i]], txt, fixed = TRUE)
  }

  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(enc2utf8(txt), tmp, useBytes = TRUE)
  html <- as.character(shiny::includeMarkdown(tmp))

  for (i in seq_along(m_env$items)) {
    esc_math <- m_env$items[[i]]
    esc_math <- gsub("&", "&amp;", esc_math, fixed = TRUE)
    esc_math <- gsub("<", "&lt;",  esc_math, fixed = TRUE)
    esc_math <- gsub(">", "&gt;",  esc_math, fixed = TRUE)
    html <- gsub(sprintf("MDMATHTOKEN%04dEND", i), esc_math, html, fixed = TRUE)
  }

  htmltools::HTML(html)
}

# Pre-render alle Engelstalige en Nederlandstalige Markdown-documenten eenmalig bij opstarten
md_cache <- list(
  hf_rationale = list(
    en = include_math_md("hf_pkpd_rationale_en.md"),
    nl = include_math_md("hf_pkpd_rationale.md")
  ),
  model_rationale = list(
    en = include_math_md("model_rationale_en.md"),
    nl = include_math_md("model_rationale.md")
  ),
  model_doc = list(
    en = include_math_md("model_doc_en.md"),
    nl = include_math_md("model_doc.md")
  ),
  validation = list(
    en = includeMarkdown("VALIDATION_en.md"),
    nl = includeMarkdown("VALIDATION.md")
  ),
  about = list(
    en = includeMarkdown("about_en.md"),
    nl = includeMarkdown("about.md")
  )
)

retypeset_script <- function() {
  tags$script(HTML("
    setTimeout(function() {
      if (typeof initCopyButtons === 'function') initCopyButtons();
      if (window.MathJax && MathJax.Hub) {
        MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
      }
    }, 0);
  "))
}

ui <- fluidPage(
  title = "Dapagliflozin PK/PD & physiological dose extrapolation — incl. paediatric heart failure",
  withMathJax(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
    tags$script(HTML("
      window.currentAppLang = 'en';
      var copyBtnI18n = {
        en: { copy: 'Copy code', copied: 'Copied \\u2713' },
        nl: { copy: 'Kopieer code', copied: 'Gekopieerd \\u2713' }
      };

      function getCopyLabels() {
        return copyBtnI18n[window.currentAppLang] || copyBtnI18n.en;
      }

      function initCopyButtons() {
        var lbls = getCopyLabels();
        $('.md-doc pre').each(function() {
          var $pre = $(this);
          if ($pre.parent().hasClass('pre-copy-wrap')) return;
          $pre.wrap('<div class=\"pre-copy-wrap\"></div>');
          var $btn = $('<button type=\"button\" class=\"btn btn-default btn-xs copy-code-btn\"></button>');
          $btn.text(lbls.copy);
          $btn.on('click', function() {
            var codeEl = $pre.find('code')[0] || $pre[0];
            var txt = codeEl.innerText || codeEl.textContent || '';
            var done = function() {
              var cur = getCopyLabels();
              $btn.text(cur.copied).addClass('copied');
              setTimeout(function() {
                $btn.text(getCopyLabels().copy).removeClass('copied');
              }, 1800);
            };
            if (navigator.clipboard && navigator.clipboard.writeText) {
              navigator.clipboard.writeText(txt).then(done).catch(function() {
                var ta = document.createElement('textarea');
                ta.value = txt; ta.style.position = 'fixed'; ta.style.opacity = '0';
                document.body.appendChild(ta); ta.select();
                try { document.execCommand('copy'); done(); } catch(e) {}
                document.body.removeChild(ta);
              });
            } else {
              var ta = document.createElement('textarea');
              ta.value = txt; ta.style.position = 'fixed'; ta.style.opacity = '0';
              document.body.appendChild(ta); ta.select();
              try { document.execCommand('copy'); done(); } catch(e) {}
              document.body.removeChild(ta);
            }
          });
          $pre.parent().append($btn);
        });
      }

      function applyClientLang(lang) {
        lang = (lang === 'nl') ? 'nl' : 'en';
        window.currentAppLang = lang;
        document.documentElement.lang = lang;
        document.title = (lang === 'nl')
          ? 'Dapagliflozine PK/PD & fysiologische dosisextrapolatie \\u2014 incl. pediatrisch hartfalen'
          : 'Dapagliflozin PK/PD & physiological dose extrapolation \\u2014 incl. paediatric heart failure';

        $('.i18n-txt').each(function() {
          var val = $(this).attr('data-' + lang);
          if (typeof val !== 'undefined') $(this).text(val);
        });
        $('.i18n-en-block').toggle(lang === 'en');
        $('.i18n-nl-block').toggle(lang === 'nl');

        var lbls = getCopyLabels();
        $('.copy-code-btn').not('.copied').text(lbls.copy);

        var $fileWrap = $('#obs_file').closest('.shiny-input-container');
        if ($fileWrap.length) {
          var $btnSpan = $fileWrap.find('.btn-file');
          if ($btnSpan.length) {
            $btnSpan.contents().filter(function() { return this.nodeType === 3; }).first()
              .replaceWith(lang === 'nl' ? 'Bladeren...' : 'Browse...');
          }
          var $roInput = $fileWrap.find('input.form-control[readonly]');
          if ($roInput.length) {
            $roInput.attr('placeholder', lang === 'nl' ? 'Geen bestand gekozen' : 'No file selected');
          }
        }
      }

      $(document).ready(function() {
        applyClientLang('en');
        initCopyButtons();
      });
      $(document).on('change', '#lang input[name=\"lang\"]', function() {
        applyClientLang($(this).val());
      });
      $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"], a[data-toggle=\"pill\"]', function() {
        initCopyButtons();
        if (window.MathJax && MathJax.Hub) {
          MathJax.Hub.Queue(['Typeset', MathJax.Hub]);
        }
      });
      function scrollToCodeSection() {
        var $target = $('#model_doc_container h2').filter(function() {
          return $(this).text().indexOf('6.') === 0;
        }).first();
        if ($target.length) {
          $('html, body').animate({ scrollTop: $target.offset().top - 20 }, 250);
        }
      }
    "))
  ),
  div(
    class = "app-header-bar",
    h2(
      class = "app-main-title",
      i18n_span(
        "Dapagliflozin PK/PD & physiological dose extrapolation — incl. paediatric heart failure",
        "Dapagliflozine PK/PD & fysiologische dosisextrapolatie — incl. pediatrisch hartfalen"
      )
    ),
    div(
      class = "lang-switch-box",
      title = "Switch language / Wissel taal (English / Nederlands)",
      span(class = "lang-switch-icon", icon("globe")),
      radioButtons(
        "lang",
        label = NULL,
        choices = c("English" = "en", "Nederlands" = "nl"),
        selected = "en",
        inline = TRUE
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4(i18n_span("Simulation settings", "Simulatie-instellingen")),
      actionButton("run", i18n_span("Simulate", "Simuleer"),
                   icon = icon("play"), class = "btn-primary btn-lg", width = "100%"),
      hr(),
      radioButtons(
        "dose_mode",
        label = i18n_span("Dosing mode", "Doseringswijze"),
        choiceNames = list(
          i18n_span("Weight-based (mg/kg)", "Op gewicht (mg/kg)"),
          i18n_span("Fixed dose (mg)", "Vaste dosis (mg)")
        ),
        choiceValues = list("mgkg", "mg"),
        selected = "mgkg",
        inline = TRUE
      ),
      conditionalPanel(
        condition = "input.dose_mode == 'mgkg'",
        numericInput(
          "dose_mgkg",
          i18n_span("Dose per administration (mg/kg) — unit dose, uncapped",
                    "Dosis per gift (mg/kg) — keerdosis, zonder plafond"),
          value = 0.15, min = 0.001, step = 0.01
        )
      ),
      conditionalPanel(
        condition = "input.dose_mode == 'mg'",
        numericInput(
          "dose_mg",
          i18n_span("Fixed dose per administration (mg) — unit dose, uncapped",
                    "Vaste dosis per gift (mg) — keerdosis, zonder plafond"),
          value = 10, min = 0.01, step = 0.5
        )
      ),
      sliderInput(
        "tau",
        i18n_span("Dosing interval (h) between doses", "Dosisinterval (uur) tussen twee doses"),
        min = 1, max = 48, value = 24, step = 1
      ),
      sliderInput(
        "n_doses",
        i18n_span("Number of doses (1 = single dose)", "Aantal doses (1 = eenmalige gift)"),
        min = 1, max = 28, value = 7, step = 1
      ),
      checkboxInput(
        "steady_state",
        i18n_span("Direct analytical steady state (when IOV = 0%)",
                  "Directe analytische steady state (bij IOV = 0%)"),
        value = TRUE
      ),
      sliderInput(
        "iov_cv",
        i18n_span("Inter-occasion variability on F1 / administered dose (%CV per dose)",
                  "Inter-occasion variabiliteit op F1 / toegediende dosis (%CV per dosis)"),
        min = 0, max = 100, value = 0, step = 5
      ),
      tags$label(i18n_span("Age range (years; e.g. 0.08 = 1 mo, 0.5 = 6 mo)",
                           "Leeftijdsrange (jaren; bv. 0.08 = 1 mnd, 0.5 = 6 mnd)")),
      fluidRow(
        column(6, numericInput("age_min", i18n_span("Minimum", "Minimum"), value = 0.1, min = 0, max = 100, step = 0.05)),
        column(6, numericInput("age_max", i18n_span("Maximum", "Maximum"), value = 1.0, min = 0, max = 100, step = 0.25))
      ),
      uiOutput("dataset_age_hint"),
      selectInput(
        "dataset",
        i18n_span("Virtual population (copulapop dataset)", "Virtuele populatie (copulapop dataset)"),
        choices = get_datasets("en"), selected = "RUMC"
      ),
      selectInput(
        "mat_option",
        i18n_span("Maturation function on hepatic + renal UGT1A9",
                  "Maturatiefunctie op hepatische + renale UGT1A9"),
        choices = get_mat_choices("en"), selected = "farhan2024"
      ),
      selectInput(
        "ka_mode",
        i18n_span("Absorption variability (BSV on Ka)", "Absorptie-variabiliteit (BSV op Ka)"),
        choices = get_ka_choices("en"), selected = "calibrated"
      ),
      div(
        class = "hf-sidebar-box",
        h4(i18n_span("Heart failure physiology & enteral tube", "Hartfalen-fysiologie & sonde")),
        selectInput(
          "hf_preset",
          i18n_span("Clinical heart failure profile", "Klinisch hartfalenprofiel"),
          choices = get_hf_preset_choices("en"), selected = "standard"
        ),
        uiOutput("hf_preset_note"),
        checkboxInput(
          "show_hf_adv",
          i18n_span("Show fine-tuning (cardiorenal, congestion, oedema, albumin, F1)",
                    "Toon fijnafstelling (cardiorenaal, stuwing, oedeem, albumine, F1)"),
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.show_hf_adv || input.hf_preset != 'standard'",
          sliderInput(
            "rf_pct",
            i18n_span("Cardiorenal GFR preservation vs norm (%)", "Cardiorenaal GFR-behoud t.o.v. norm (%)"),
            min = 25, max = 120, value = 100, step = 5
          ),
          sliderInput(
            "hep_pct",
            i18n_span("Hepatic UGT1A9 function in congestion (%)", "Hepatische UGT1A9-functie bij stuwing (%)"),
            min = 30, max = 100, value = 100, step = 5
          ),
          sliderInput(
            "edema_pct",
            i18n_span("Fluid overload / oedema above dry weight (%)", "Overvulling / oedeem boven drooggewicht (%)"),
            min = 0, max = 30, value = 0, step = 1
          ),
          conditionalPanel(
            condition = "input.edema_pct > 0 && input.dose_mode == 'mgkg'",
            radioButtons(
              "dose_wt_basis",
              label = i18n_span("Weight basis for mg/kg dosing", "Gewichtsbasis voor mg/kg-dosering"),
              choiceNames = list(
                i18n_span("Dry weight (recommended)", "Drooggewicht (aanbevolen)"),
                i18n_span("Actual oedematous weight", "Actueel oedemateus gewicht")
              ),
              choiceValues = list("dry", "actual"),
              selected = "dry"
            )
          ),
          radioButtons(
            "alb_mode",
            label = i18n_span("Serum albumin (protein binding fu)", "Serumalbumine (eiwitbinding fu)"),
            choiceNames = list(
              i18n_span("Physiological by age (McNamara)", "Fysiologisch naar leeftijd (McNamara)"),
              i18n_span("Manual (hypoalbuminaemia / PLE)", "Handmatig (hypoalbuminemie / PLE)")
            ),
            choiceValues = list("auto", "manual"),
            selected = "auto"
          ),
          conditionalPanel(
            condition = "input.alb_mode == 'manual'",
            sliderInput(
              "alb_val",
              i18n_span("Serum albumin (g/L; adult ref = 44)", "Serumalbumine (g/L; volw. ref = 44)"),
              min = 15, max = 50, value = 30, step = 1
            )
          ),
          sliderInput(
            "f1_pct",
            i18n_span("Mean relative bioavailability F1 (tube/suspension, %)",
                      "Gemiddelde relatieve opname F1 (sonde/suspensie, %)"),
            min = 30, max = 120, value = 100, step = 5
          )
        )
      ),
      selectInput("sex", i18n_span("Sex", "Geslacht"), choices = get_sex_choices("en")),
      sliderInput(
        "n_subj",
        i18n_span("Number of virtual individuals (N)", "Aantal virtuele individuen (N)"),
        min = 50, max = 2000, value = 1000, step = 50
      ),
      numericInput("seed", i18n_span("Random seed", "Random seed"), value = 42, min = 1),
      hr(),
      h4(i18n_span("Display", "Weergave")),
      radioButtons(
        "conc_scale",
        label = i18n_span("Concentration on PK curve", "Concentratie op PK-curve"),
        choiceNames = list(
          i18n_span("Total plasma C(t)", "Totaal plasma C(t)"),
          i18n_span("Unbound plasma Cu(t) + SGLT2 threshold", "Ongebonden vrij plasma Cu(t) + SGLT2-drempel")
        ),
        choiceValues = list("total", "unbound"),
        selected = "total"
      ),
      checkboxInput("show_indiv", i18n_span("Show individual curves (max. 30)", "Toon individuele curves (max. 30)"), value = FALSE),
      checkboxInput("show_day1", i18n_span("Show median day 1 (first dose)", "Toon mediaan dag 1 (eerste dosis)"), value = FALSE),
      checkboxInput("log_y", i18n_span("Logarithmic y-axis", "Logaritmische y-as"), value = FALSE),
      hr(),
      h4(i18n_span("Observed data (overlay)", "Gemeten data (overlay)")),
      fileInput(
        "obs_file",
        i18n_span("Upload CSV/TSV (columns TAD or TIME + DV or CONC) — or paste below",
                  "Upload CSV/TSV (kolommen TAD of TIME + DV of CONC) — of plak hieronder"),
        accept = c(".csv", ".tsv", ".txt")
      ),
      textAreaInput(
        "obs_text",
        i18n_span("Or paste observed data here (e.g. TAD,DV per line)",
                  "Of plak gemeten data hier (bv. TAD,DV per regel)"),
        rows = 5,
        placeholder = "TAD,DV\n0.5,12.3\n1,88.1\n2,101.5\n24,3.9"
      ),
      checkboxInput(
        "show_obs",
        i18n_span("Show observed data as overlay on curve", "Toon gemeten data als overlay op de curve"),
        value = TRUE
      )
    ),
    mainPanel(
      width = 9,
      uiOutput("error_box"),
      uiOutput("warn_box"),
      i18n_div(
        p("This app simulates the pharmacokinetics and semi-mechanistic pharmacodynamics of dapagliflozin in a virtual population (",
          strong("copulapop"), ") from neonates to older adults, including pathophysiological adjustments and dose prediction for ",
          strong("paediatric heart failure"), ". Configure dose, age window, and an optional heart failure profile on the left and click ",
          strong("Simulate"), "."),
        p("Deze app simuleert de farmacokinetiek én semi-mechanistische farmacodynamiek van dapagliflozine in een virtuele populatie (",
          strong("copulapop"), ") van pasgeborene tot oudere, inclusief pathofysiologische aanpassingen en dosisvoorspelling voor ",
          strong("pediatrisch hartfalen"), ". Stel links dosis, leeftijdsvenster en eventueel een hartfalenprofiel in en klik op ",
          strong("Simuleer"), ".")
      ),
      tabsetPanel(
        id = "tabs",
        tabPanel(
          title = i18n_span("PK curve", "PK-curve"),
          value = "tab_pk",
          plotOutput("pk_curve", height = "580px"),
          helpText(i18n_span(
            "The shaded band shows the 90% prediction interval across individuals (P5-P95). Only inter-individual variability (IIV, plus optional per-dose IOV) is simulated; residual error ($SIGMA / EPS) is NOT simulated.",
            "De band is een 90% predictie-interval over individuen (P5-P95). Alleen inter-individuele variabiliteit (IIV, plus optionele IOV per dosis) wordt gesimuleerd; residuele error ($SIGMA / EPS) wordt NIET gesimuleerd."
          )),
          hr(),
          wellPanel(
            h4(i18n_span("Download publication-ready figure", "Publication-ready figuur downloaden")),
            selectInput("fig_format", i18n_span("Format", "Formaat"),
                        choices = get_fig_format_choices("en"), selected = "png"),
            numericInput("fig_width", i18n_span("Width (inches)", "Breedte (inch)"), value = 7, min = 2, max = 24, step = 0.5),
            numericInput("fig_height", i18n_span("Height (inches)", "Hoogte (inch)"), value = 5, min = 2, max = 24, step = 0.5),
            numericInput("fig_dpi", i18n_span("Resolution (dpi; raster formats)", "Resolutie (dpi; rasterformaten)"), value = 600, min = 72, max = 1200, step = 50),
            downloadButton("dl_fig", i18n_span("Download figure (incl. observed overlay if checked)",
                                               "Download figuur (incl. gemeten overlay indien aangevinkt)"))
          ),
          hr(),
          downloadButton("dl_curve", i18n_span("Download curve summary (CSV)", "Download curve-samenvatting (CSV)")),
          downloadButton("dl_metrics", i18n_span("Download individual PK/PD metrics (CSV)", "Download individuele PK/PD-metrieken (CSV)"))
        ),
        tabPanel(
          title = i18n_span("Heart failure: PK/PD & Dose prediction", "Hartfalen: PK/PD & Dosisvoorspelling"),
          value = "tab_hf",
          h4(i18n_span("Dose selection & PK–PD analysis in paediatric heart failure",
                       "Dosisselectie & PK–PD analyse bij pediatrisch hartfalen")),
          uiOutput("hf_kpi_cards"),
          uiOutput("hf_alerts_ui"),
          div(
            class = "hf-subtabs",
            tabsetPanel(
              id = "hf_inner_tabs",
              type = "pills",
              tabPanel(
                title = i18n_span("1. Therapeutic window & References", "1. Therapeutisch venster & Referenties"),
                value = "hf_sub_1",
                i18n_div(
                  div(
                    class = "hf-intro-box",
                    strong("Exposure matching (ICH E11A) to adult heart failure trials: "),
                    "In ", strong("DAPA-HF"), " and ", strong("DELIVER"), ", a fixed maintenance dose of ",
                    strong("10 mg once daily"), " proved efficacious and safe across a broad eGFR spectrum. In the absence of paediatric phase-III outcome trials, dose extrapolation targets ",
                    strong("exposure matching"), " to the unbound steady-state exposure range of these adult reference populations.",
                    tags$br(), tags$br(),
                    strong("Why guide dosing by unbound exposure (AUCu,24h,ss) rather than total plasma concentrations? "),
                    "Dapagliflozin is extensively bound to serum albumin (", code("91.4%"), "; ", code("fu,adult = 0.086"),
                    ") and behaves as a low-extraction drug (", code("ER < 0.05"), "). ",
                    "In physiological neonatal hypoalbuminaemia, venous congestion, cardiac cachexia, or ", em("protein-losing enteropathy"), " (Fontan), ", code("fu"),
                    " increases: apparent total clearance ", code("CL/F"), " rises proportionally and total plasma concentration ", code("C(t)"),
                    " falls, whereas pharmacologically active unbound concentration ", code("Cu(t)"), " (filtered at the glomerulus toward luminal SGLT2) remains unchanged."
                  ),
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
                  )
                ),
                plotOutput("hf_window_plot", height = "370px"),
                uiOutput("hf_window_cards_ui"),
                checkboxInput(
                  "show_custom_window",
                  i18n_span("Manually adjust target window and safety ceiling (sensitivity analysis)",
                            "Pas streefvenster en veiligheidsplafond handmatig aan (gevoeligheidsanalyse)"),
                  value = FALSE
                ),
                conditionalPanel(
                  condition = "input.show_custom_window",
                  wellPanel(
                    fluidRow(
                      column(4, numericInput("hf_target_lo",
                                             i18n_span("Target lower bound (ng*h/mL; ref 20)", "Ondergrens streefvenster (ng*h/mL; ref 20)"),
                                             value = 20, min = 5, max = 60, step = 2)),
                      column(4, numericInput("hf_target_hi",
                                             i18n_span("Target upper bound (ng*h/mL; ref 75)", "Bovengrens streefvenster (ng*h/mL; ref 75)"),
                                             value = 75, min = 30, max = 150, step = 5)),
                      column(4, numericInput("hf_safety_max",
                                             i18n_span("Safety ceiling P95 (ng*h/mL; ref 90)", "Veiligheidsplafond P95 (ng*h/mL; ref 90)"),
                                             value = 90, min = 50, max = 200, step = 5))
                    ),
                    helpText(i18n_span(
                      "Adjustments analytically update all PTA curves, ontogeny comparisons, and tables instantly without restarting the Monte Carlo simulation.",
                      "Aanpassingen werken alle PTA-curves, ontogenievergelijkingen en tabellen direct analytisch bij zonder herstart van de Monte Carlo-simulatie."
                    ))
                  )
                )
              ),
              tabPanel(
                title = i18n_span("2. Dose selection (PTA) & Age table", "2. Dosisselectie (PTA) & Leeftijdstabel"),
                value = "hf_sub_2",
                i18n_div(
                  div(
                    class = "hf-intro-box",
                    strong("Dose optimisation & UGT1A9 ontogeny sensitivity: "),
                    "Panel A shows the ", strong("Probability of Target Attainment (PTA)"), " for the simulated cohort: the percentage of patients inside the target window (",
                    code("AUCu,24h,ss"), " between lower and upper bounds) subject to keeping ", code("P95"), " below the safety ceiling. ",
                    "Panel B quantifies structural model uncertainty in postnatal ", strong("UGT1A9 glucuronidation"), ": below 6 months, rapid ", em("in-vivo"), " maturation (",
                    em("Farhan 2024"), ") and conservative PBPK ontogeny (", em("Simcyp s1"), ") differ by a factor of 3 to 4, supporting a cautious starting dose with stepwise up-titration."
                  ),
                  div(
                    class = "hf-intro-box",
                    strong("Dosisoptimalisatie & gevoeligheid voor UGT1A9-ontogenie: "),
                    "Figuur A toont de ", strong("Probability of Target Attainment (PTA)"), " voor het gesimuleerde cohort: het percentage patiënten binnen het streefvenster (",
                    code("AUCu,24h,ss"), " tussen onder- en bovengrens) onder de randvoorwaarde dat ", code("P95"), " onder het veiligheidsplafond blijft. ",
                    "Figuur B kwantificeert de structurele modelonzekerheid in postnatale ", strong("UGT1A9-glucuronidering"), ": onder 6 maanden wijken snelle ", em("in-vivo"), " rijping (",
                    em("Farhan 2024"), ") en conservatieve PBPK-ontogenie (", em("Simcyp s1"), ") een factor 3 tot 4 uiteen, wat een voorzichtige startdosis met stapsgewijze optitratie onderbouwt."
                  )
                ),
                plotOutput("hf_pta_ontogeny_plot", height = "430px"),
                h5(strong(i18n_span(
                  "Exposure and dosing guidance per UGT1A9 ontogeny model for the current cohort:",
                  "Expositie en dosisadvies per UGT1A9-ontogeniemodel voor het huidige cohort:"
                ))),
                tableOutput("hf_ontogeny_table"),
                hr(),
                h4(i18n_span("Age-stratified reference overview (0 mo – 18 years)",
                             "Leeftijdsgestratificeerd referentie-overzicht (0 mnd – 18 jaar)")),
                i18n_div(
                  p("The table below evaluates standard paediatric age strata under the heart failure and formulation profile selected on the left. ",
                    strong("Allometric implication above 1 year: "),
                    "because metabolic UGT1A9 clearance scales with ", code("NFM^0.75"), ", weight-normalised clearance in toddlers and young children exceeds adult clearance, requiring ",
                    code("~0.18–0.20 mg/kg/day"), " on euvolaemic dry weight to match adult ", code("10 mg"), " exposure (", code("10 mg / 84 kg = 0.12 mg/kg"), ")."),
                  p("Onderstaande tabel evalueert standaardleeftijdsstrata onder het links ingestelde hartfalen- en formuleringsprofiel. ",
                    strong("Allometrische consequentie boven 1 jaar: "),
                    "doordat metabole UGT1A9-klaring schaalt met ", code("NFM^0,75"), ", is de gewichtsgenormaliseerde klaring bij peuters en kleuters hoger dan bij volwassenen en is ",
                    code("~0,18–0,20 mg/kg/dag"), " op euvolemisch drooggewicht nodig voor een volwassen ", code("10 mg"), "-equivalente blootstelling (", code("10 mg / 84 kg = 0,12 mg/kg"), ").")
                ),
                tableOutput("hf_lifespan_table"),
                downloadButton("dl_hf_regimen",
                               i18n_span("Download age-stratified HF reference table (CSV)",
                                         "Download leeftijdsgestratificeerde HF-referentietabel (CSV)"))
              ),
              tabPanel(
                title = i18n_span("3. PK–PD Relationship & Clinical Safety", "3. PK–PD Relatie & Klinische Veiligheid"),
                value = "hf_sub_3",
                i18n_div(
                  div(
                    class = "hf-explain-grid-3",
                    div(
                      class = "hf-explain-card",
                      div(class = "hf-explain-title", "1. Dual pharmacodynamics (luminal vs. systemic)"),
                      p("Dapagliflozin inhibits ", strong("luminal SGLT2"), " in the S1/S2 proximal tubule following glomerular filtration of ", code("Cu(t)"),
                        " (", code("IC50 = 1.1 nM"), "; tubular concentration factor ", code("2.25x"), " → ", code("Cu,50 = 0.20 ng/mL"),
                        "), alongside ", strong("systemic / tubuloglomerular Emax effects"), " on intraglomerular pressure, ketone metabolism, and erythropoiesis.")
                    ),
                    div(
                      class = "hf-explain-card hf-zone-high",
                      div(class = "hf-explain-title", "2. PK–PD dissociation at low GFR (ceiling effect)"),
                      p("In immature nephrons or cardiorenal syndrome, clearance declines (", code("AUCu,24h"), " rises), yet low ",
                        strong("filtered sodium and glucose load"), " limits osmotic natriuresis while luminal SGLT2 is already ", code("85–92%"), " saturated. ",
                        strong("Up-titrating to force additional diuresis is pharmacologically ineffective."))
                    ),
                    div(
                      class = "hf-explain-card hf-zone-max",
                      div(class = "hf-explain-title", "3. Shifted therapeutic index (< 2 years)"),
                      p("Young children have limited hepatic glycogen stores: normoglycaemic glucosuria during fasting or intercurrent infection carries a risk of ",
                        strong("hypoglycaemia and euglycaemic ketoacidosis"), " (", em("sick-day rules"), "). ",
                        "Also monitor for prerenal AKI (alongside loop diuretics/RAS blockade) and renal caloric loss in ", em("failure to thrive"), ".")
                    )
                  ),
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
                  )
                ),
                plotOutput("hf_pkpd_plot", height = "445px"),
                helpText(i18n_span(
                  "Note on Panel B: each point represents a simulated patient coloured by eGFR. Patients with reduced eGFR (red/orange) exhibit high unbound exposure and maximal systemic Emax response (black curves), yet a blunted effective tubular natriuretic response on the y-axis due to low glomerular filtration.",
                  "Toelichting bij figuur B: elk punt representeert een gesimuleerde patiënt, gekleurd naar eGFR. Patiënten met verlaagde eGFR (rood/oranje) vertonen hoge ongebonden blootstelling en maximale systemische Emax-respons (zwarte curves), maar een begrensde effectieve tubulaire natriuretische respons op de y-as door lage glomerulaire filtratie."
                ))
              ),
              tabPanel(
                title = i18n_span("4. Methodological Rationale & References", "4. Methodologische Onderbouwing & Bronnen"),
                value = "hf_sub_4",
                uiOutput("hf_rationale_md_ui")
              )
            )
          )
        ),
        tabPanel(
          title = i18n_span("Exposure metrics (PK & PD)", "Expositiemetrieken (PK & PD)"),
          value = "tab_metrics",
          h4(i18n_span("Population summary: median, geometric mean, and [90% prediction interval P5–P95]",
                       "Populatiesamenvatting: mediaan, geometrisch gemiddelde en [90% predictie-interval P5–P95]")),
          tableOutput("ss_table"),
          plotOutput("metric_plots", height = "380px"),
          plotOutput("auc_age", height = "380px")
        ),
        tabPanel(
          title = i18n_span("AUC & Cmax vs renal function / age", "AUC & Cmax vs nierfunctie / leeftijd"),
          value = "tab_cov",
          plotOutput("exposure_covariates", height = "720px")
        ),
        tabPanel(
          title = i18n_span("Virtual population", "Virtuele populatie"),
          value = "tab_pop",
          verbatimTextOutput("pop_info"),
          plotOutput("pop_plots", height = "360px"),
          h4(i18n_span("First 100 individuals from the virtual population (incl. dry weight, FFM, fu, and effective eGFR)",
                       "Eerste 100 individuen uit de virtuele populatie (incl. drooggewicht, FFM, fu en effectieve eGFR)")),
          tableOutput("pop_table")
        ),
        tabPanel(
          title = i18n_span("Model & assumptions", "Model & assumpties"),
          value = "tab_model",
          div(
            class = "hf-subtabs",
            tabsetPanel(
              id = "model_inner_tabs",
              type = "pills",
              tabPanel(
                title = i18n_span("1. Step-by-Step Physiological Derivation & Calibration",
                                  "1. Stapsgewijze Fysiologische Afleiding & Kalibratie"),
                value = "model_sub_1",
                uiOutput("model_rationale_md_ui")
              ),
              tabPanel(
                title = i18n_span("2. Complete Mathematical Equations & Reproduction (R / NONMEM)",
                                  "2. Volledige Wiskundige Formules & Reproductie (R / NONMEM)"),
                value = "model_sub_2",
                wellPanel(
                  strong(i18n_span("Jump directly to R & NONMEM reproduction code (§6): ",
                                   "Direct naar R- & NONMEM-reproductiecode (§6): ")),
                  tags$button(
                    type = "button",
                    class = "btn btn-primary btn-sm",
                    onclick = "scrollToCodeSection()",
                    i18n_span("Jump to §6: R & NONMEM code (copy-paste)",
                              "Spring naar §6: R- & NONMEM-code (copy-paste)")
                  ),
                  " ",
                  downloadButton("dl_nm_ctl",
                                 i18n_span("Download NONMEM control stream (.ctl)",
                                           "Download NONMEM control stream (.ctl)"),
                                 class = "btn-sm"),
                  " ",
                  downloadButton("dl_nm_csv",
                                 i18n_span("Download 1-child example dataset (.csv)",
                                           "Download 1-kind voorbeelddataset (.csv)"),
                                 class = "btn-sm")
                ),
                uiOutput("model_doc_md_ui")
              )
            )
          )
        ),
        tabPanel(
          title = i18n_span("Validation", "Validatie"),
          value = "tab_val",
          div(
            class = "md-doc",
            wellPanel(
              strong(i18n_span("Full validation report (English, incl. all tables and figures): ",
                               "Volledig validatierapport (Engels, incl. alle tabellen en figuren): ")),
              tags$a(i18n_span("Open PDF report (VALIDATION_REPORT.pdf)",
                               "Open PDF-rapport (VALIDATION_REPORT.pdf)"),
                     href = "docs/VALIDATION_REPORT.pdf",
                     target = "_blank", class = "btn btn-default btn-sm"),
              " ",
              tags$a(i18n_span("Open HTML report (VALIDATION_REPORT.html)",
                               "Open HTML-rapport (VALIDATION_REPORT.html)"),
                     href = "docs/VALIDATION_REPORT.html",
                     target = "_blank", class = "btn btn-default btn-sm")
            ),
            uiOutput("validation_md_ui"),
            hr(),
            h3(i18n_span("Detailed validation figures (from the full report)",
                         "Gedetailleerde validatiefiguren (uit het volledige rapport)")),
            tags$img(src = "figures/fig1_penland_egfr.svg", alt = "Figure 1: Penland 2023 CKD strata"),
            tags$img(src = "figures/fig2_tirucherai_dose.svg", alt = "Figure 2: Tirucherai 2016 paediatric T2DM"),
            tags$img(src = "figures/fig3_pbpk_age.svg", alt = "Figure 3: Jo 2021 paediatric PBPK"),
            tags$img(src = "figures/fig3_curves.svg", alt = "Figure 4: Paediatric PK profiles"),
            tags$img(src = "figures/fig5_maturation.svg", alt = "Figure 5: Developmental ontogeny")
          )
        ),
        tabPanel(
          title = i18n_span("About", "Over deze app"),
          value = "tab_about",
          uiOutput("about_md_ui")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  sim_res  <- reactiveVal(NULL)
  err_msg  <- reactiveVal(NULL)
  warn_msg <- reactiveVal(NULL)

  cur_lang <- reactive({
    if (identical(input$lang, "nl")) "nl" else "en"
  })

  prev_hf_preset <- reactiveVal("standard")
  prev_dataset   <- reactiveVal("RUMC")

  # Werk selectInput-keuzelijsten bij wanneer de taal wisselt, met behoud van de huidige selectie
  observeEvent(input$lang, ignoreInit = TRUE, {
    lg <- cur_lang()
    updateSelectInput(session, "dataset", choices = get_datasets(lg), selected = isolate(input$dataset))
    updateSelectInput(session, "mat_option", choices = get_mat_choices(lg), selected = isolate(input$mat_option))
    updateSelectInput(session, "ka_mode", choices = get_ka_choices(lg), selected = isolate(input$ka_mode))
    updateSelectInput(session, "hf_preset", choices = get_hf_preset_choices(lg), selected = isolate(input$hf_preset))
    updateSelectInput(session, "sex", choices = get_sex_choices(lg), selected = isolate(input$sex))
    updateSelectInput(session, "fig_format", choices = get_fig_format_choices(lg), selected = isolate(input$fig_format))
  })

  # Markdown-documenten op basis van geselecteerde taal
  output$hf_rationale_md_ui <- renderUI({
    div(class = "md-doc", md_cache$hf_rationale[[cur_lang()]], retypeset_script())
  })
  output$model_rationale_md_ui <- renderUI({
    div(class = "md-doc", md_cache$model_rationale[[cur_lang()]], retypeset_script())
  })
  output$model_doc_md_ui <- renderUI({
    div(id = "model_doc_container", class = "md-doc", md_cache$model_doc[[cur_lang()]], retypeset_script())
  })
  output$validation_md_ui <- renderUI({
    md_cache$validation[[cur_lang()]]
  })
  output$about_md_ui <- renderUI({
    div(class = "md-doc", md_cache$about[[cur_lang()]])
  })

  pick_bilingual <- function(item, lg) {
    if (is.list(item) && !is.null(item[[lg]])) item[[lg]]
    else if (is.character(item)) item
    else ""
  }

  output$error_box <- renderUI({
    e <- err_msg()
    if (is.null(e)) return(NULL)
    lg <- cur_lang()
    div(class = "alert alert-danger",
        strong(if (lg == "nl") "Fout: " else "Error: "),
        pick_bilingual(e, lg))
  })

  output$warn_box <- renderUI({
    w <- warn_msg()
    if (is.null(w) || length(w) == 0) return(NULL)
    lg <- cur_lang()
    lead_lbl <- if (lg == "nl") "Let op: " else "Note: "
    msgs <- unique(vapply(w, pick_bilingual, character(1), lg = lg))
    msgs <- msgs[nzchar(msgs)]
    if (length(msgs) == 0) return(NULL)
    div(class = "alert alert-warning",
        lapply(msgs, function(m) div(strong(lead_lbl), m)))
  })

  output$dataset_age_hint <- renderUI({
    info <- copulapop_dataset_info(input$dataset)
    if (is.null(info)) return(NULL)
    if (cur_lang() == "nl") {
      helpText(sprintf("Beschikbaar in %s: %g tot %g jaar.",
                       input$dataset, info$age_range[1], info$age_range[2]))
    } else {
      helpText(sprintf("Available in %s: %g to %g years.",
                       input$dataset, info$age_range[1], info$age_range[2]))
    }
  })

  output$hf_preset_note <- renderUI({
    pr <- hf_presets(cur_lang())[[input$hf_preset]]
    if (is.null(pr) || is.null(pr$note)) return(NULL)
    helpText(pr$note)
  })

  # Synchroniseer hartfalen-fysiologie schuifjes alleen wanneer daadwerkelijk een ander klinisch preset wordt gekozen
  observeEvent(input$hf_preset, ignoreInit = TRUE, {
    if (identical(input$hf_preset, prev_hf_preset())) return()
    prev_hf_preset(input$hf_preset)
    pr <- hf_presets("en")[[input$hf_preset]]
    if (!is.null(pr) && !is.null(pr$rf_factor)) {
      updateSliderInput(session, "rf_pct", value = round(100 * pr$rf_factor))
      updateSliderInput(session, "hep_pct", value = round(100 * pr$cl_hep_factor))
      updateSliderInput(session, "edema_pct", value = pr$edema_pct)
      updateRadioButtons(session, "alb_mode", selected = pr$alb_mode)
      updateSliderInput(session, "alb_val", value = pr$alb_val)
      updateSliderInput(session, "f1_pct", value = round(100 * pr$f1_rel))
    }
  })

  # Pas leeftijdsgrenzen alleen aan als een andere dataset wordt gekozen en de grenzen buiten bereik vallen
  observeEvent(input$dataset, ignoreInit = TRUE, {
    if (identical(input$dataset, prev_dataset())) return()
    prev_dataset(input$dataset)
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

  raise_bilingual_stop <- function(en, nl) {
    cond <- structure(
      list(message = en, en = en, nl = nl, call = NULL),
      class = c("bilingual_error", "error", "condition")
    )
    stop(cond)
  }

  do_simulate <- function() {
    err_msg(NULL)
    warn_msg(NULL)
    vals <- isolate(reactiveValuesToList(input))
    lg   <- if (identical(vals$lang, "nl")) "nl" else "en"
    seed <- suppressWarnings(as.integer(vals$seed))
    if (length(seed) != 1L || is.na(seed) || seed < 1L) {
      err_msg(list(
        en = "Random seed must be a positive integer.",
        nl = "Random seed moet een positief geheel getal zijn."
      ))
      return(invisible(NULL))
    }

    warns <- list()
    tryCatch({
      d_mode <- if (length(vals$dose_mode) == 1L && nzchar(vals$dose_mode)) vals$dose_mode else "mgkg"
      d_mgkg <- NULL
      d_mg   <- NULL
      if (d_mode == "mg") {
        val_mg <- suppressWarnings(as.numeric(vals$dose_mg))
        if (length(val_mg) != 1L || is.na(val_mg) || val_mg <= 0) {
          raise_bilingual_stop(
            "Fixed dose per administration (mg) must be a positive number.",
            "Vaste dosis per gift (mg) moet een positief getal zijn."
          )
        }
        d_mg <- val_mg
      } else {
        val_mgkg <- suppressWarnings(as.numeric(vals$dose_mgkg))
        if (length(val_mgkg) != 1L || is.na(val_mgkg) || val_mgkg <= 0) {
          raise_bilingual_stop(
            "Dose per administration (mg/kg) must be a positive number.",
            "Dosis per gift (mg/kg) moet een positief getal zijn."
          )
        }
        d_mgkg <- val_mgkg
      }

      a_min <- suppressWarnings(as.numeric(vals$age_min))
      a_max <- suppressWarnings(as.numeric(vals$age_max))
      if (length(a_min) != 1L || length(a_max) != 1L || is.na(a_min) || is.na(a_max) || a_min < 0 || a_max < 0) {
        raise_bilingual_stop(
          "Minimum and maximum age must be non-negative numbers.",
          "Minimum- en maximumleeftijd moeten niet-negatieve getallen zijn."
        )
      }
      if (a_min > a_max) {
        raise_bilingual_stop(
          "Minimum age cannot be greater than maximum age.",
          "Minimumleeftijd mag niet groter zijn dan maximumleeftijd."
        )
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

      prog_pop <- if (lg == "nl") "Virtuele populatie genereren (copulapop) ..." else "Generating virtual population (copulapop) ..."
      pop <- withCallingHandlers(
        withProgress(message = prog_pop, value = 0.25, {
          generate_virtual_population(ds_name, n_subj, c(a_min, a_max), sex_val, seed, lang = lg)
        }),
        warning = function(w) {
          m_txt <- conditionMessage(w)
          warns[[length(warns) + 1L]] <<- list(en = m_txt, nl = m_txt)
          invokeRestart("muffleWarning")
        }
      )

      use_ss <- isTRUE(vals$steady_state) && iov_val == 0
      if (isTRUE(vals$steady_state) && iov_val > 0) {
        warns[[length(warns) + 1L]] <- list(
          en = sprintf(
            "Direct analytical steady state is only defined at IOV = 0%%; simulated %d doses with %g%% IOV instead.",
            nd_val, iov_val
          ),
          nl = sprintf(
            "Directe analytische steady state is alleen gedefinieerd bij IOV = 0%%; er is gesimuleerd met %d doses en %g%% IOV.",
            nd_val, iov_val
          )
        )
      }

      prog_sim <- if (lg == "nl") "Monte Carlo PK/PD-simulatie ..." else "Monte Carlo PK/PD simulation ..."
      res <- withProgress(message = prog_sim, value = 0.70, {
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
        med_th <- stats::median(res$metrics$t_half)
        warns[[length(warns) + 1L]] <- list(
          en = sprintf(
            "After %d doses (interval %g h), median %.0f%% of steady-state exposure has been reached (median t1/2,beta = %.1f h). Increase 'Number of doses' or check 'Direct analytical steady state' for full steady state.",
            nd_val, tau_val, 100 * res$settings$ss_attainment, med_th
          ),
          nl = sprintf(
            "Na %d doses (interval %g u) is mediaan %.0f%% van steady-state blootstelling bereikt (mediane t1/2,beta = %.1f u). Verhoog 'Aantal doses' of vink 'Directe analytische steady state' aan voor volledige steady state.",
            nd_val, tau_val, 100 * res$settings$ss_attainment, med_th
          )
        )
      }
      if (ed_pct > 0 && identical(wt_bas, "actual") && d_mode == "mgkg") {
        warns[[length(warns) + 1L]] <- list(
          en = sprintf(
            "Dosing is currently based on actual oedematous weight (+%g%% fluid overload) while UGT1A9 clearance scales with dry weight: this increases exposure by +%g%%.",
            ed_pct, ed_pct
          ),
          nl = sprintf(
            "Er wordt nu op actueel oedemateus gewicht (+%g%% overvulling) gedoseerd terwijl UGT1A9-klaring schaalt met het drooggewicht: dit verhoogt de blootstelling met +%g%%.",
            ed_pct, ed_pct
          )
        )
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
      warn_msg(warns)

      dose_txt <- if (d_mode == "mg") sprintf("%.3g mg", d_mg) else sprintf("%.3g mg/kg", d_mgkg)
      notif_txt <- if (lg == "nl") {
        sprintf("Simulatie voltooid: N = %d, dosis %s per gift, RF %d%%, seed %d",
                nrow(pop), dose_txt, round(100 * rf_fac), seed)
      } else {
        sprintf("Simulation completed: N = %d, dose %s per administration, RF %d%%, seed %d",
                nrow(pop), dose_txt, round(100 * rf_fac), seed)
      }
      showNotification(notif_txt, type = "message", duration = 5)
    }, error = function(e) {
      if (inherits(e, "bilingual_error")) {
        err_msg(list(en = e$en, nl = e$nl))
        msg_show <- if (lg == "nl") e$nl else e$en
      } else {
        msg_show <- conditionMessage(e)
        err_msg(list(en = msg_show, nl = msg_show))
      }
      showNotification(paste(if (lg == "nl") "Fout:" else "Error:", msg_show),
                       type = "error", duration = 10)
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
    lg <- isolate(cur_lang())
    if (!is.null(input$obs_text) && nzchar(trimws(input$obs_text))) {
      d <- parse_obs_text(input$obs_text)
      if (is.null(d)) {
        showNotification(
          if (lg == "nl") {
            "Handmatig geplakte data: geen geldige kolommen TAD/TIME + DV/CONC gevonden."
          } else {
            "Manually pasted data: no valid TAD/TIME + DV/CONC columns found."
          },
          type = "error", duration = 8
        )
      }
      return(d)
    }
    f <- input$obs_file
    if (is.null(f)) return(NULL)
    txt <- tryCatch(paste(readLines(f$datapath, warn = FALSE), collapse = "\n"), error = function(e) "")
    d <- parse_obs_text(txt)
    if (is.null(d)) {
      showNotification(
        if (lg == "nl") {
          "Geupload bestand: vereist een tijd-kolom (TAD of TIME) en een concentratie-kolom (DV of CONC)."
        } else {
          "Uploaded file: requires a time column (TAD or TIME) and a concentration column (DV or CONC)."
        },
        type = "error", duration = 8
      )
    }
    d
  })

  build_curve_plot <- function(s, conc_scale, show_indiv, show_day1, log_y, obs = NULL, pub = FALSE, lang = "en") {
    use_en <- isTRUE(pub) || !identical(lang, "nl")
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
      do.call(labs, if (use_en) list(
        x = "Time after dose (h)",
        y = ylab_en,
        title = sprintf("%sconcentration-time profile (%s, interval %g h)",
                        if (is_u) "Unbound " else "Total ", regimen_en, st$tau),
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
      p <- p +
        geom_hline(yintercept = 0.20, color = "#1a7f37", linetype = "dotted", linewidth = 0.8) +
        geom_hline(yintercept = 1.80, color = "#1a7f37", linetype = "dotdash", linewidth = 0.7) +
        annotate("text", x = st$tau * 0.98, y = 0.22, hjust = 1, vjust = -0.3,
                 label = if (use_en) "50% luminal SGLT2 inh. (Cu,50 = 0.20 ng/mL)" else "50% luminale SGLT2-remming (Cu,50 = 0.20 ng/mL)",
                 color = "#1a7f37", size = 3.2, fontface = "bold") +
        annotate("text", x = st$tau * 0.98, y = 1.85, hjust = 1, vjust = -0.3,
                 label = if (use_en) "90% luminal SGLT2 inh. (Cu,90 = 1.80 ng/mL)" else "90% luminale SGLT2-remming (Cu,90 = 1.80 ng/mL)",
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
                   label = if (use_en) "orange = median day 1" else "oranje = mediaan dag 1",
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
          sprintf(if (use_en) "orange = obs x median fu %.3f (n = %d)" else "oranje = gemeten data x mediane fu %.3f (n = %d)", med_fu, nrow(obs_plot))
        } else {
          sprintf(if (use_en) "orange points = observed data (n = %d)" else "oranje punten = gemeten data (n = %d)", nrow(obs_plot))
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
    build_curve_plot(s, input$conc_scale, input$show_indiv, input$show_day1, input$log_y,
                     obs = obs, pub = FALSE, lang = cur_lang())
  })

  output$dl_fig <- downloadHandler(
    filename = function() {
      paste0(if (cur_lang() == "nl") "dapagliflozine_pk_curve." else "dapagliflozin_pk_curve.", input$fig_format)
    },
    content = function(file) {
      s <- sim_res(); req(s)
      obs <- if (isTRUE(input$show_obs)) obs_df() else NULL
      p <- build_curve_plot(s, input$conc_scale, input$show_indiv, input$show_day1, input$log_y,
                            obs = obs, pub = TRUE, lang = cur_lang())
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
                            target_lo = tg$target_lo, target_hi = tg$target_hi, safety_max = tg$safety_max,
                            lang = cur_lang())
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
      safety_max    = tg$safety_max,
      lang          = cur_lang()
    )
  })

  output$hf_kpi_cards <- renderUI({
    s   <- sim_res(); req(s)
    pta <- hf_pta_data(); req(pta)
    ont <- hf_ontogeny_data(); req(ont)
    tg  <- hf_targets_reactive()
    m   <- s$res$metrics
    st  <- s$res$settings
    is_nl <- (cur_lang() == "nl")

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
      if (is_nl) {
        sprintf("~%.2f mg/gift bij mediaan drooggewicht (%.1f kg) — equivalent aan volwassen 10 mg DAPA-HF expositie", targ_hi_mg, med_wt)
      } else {
        sprintf("~%.2f mg/dose at median dry weight (%.1f kg) — equivalent to adult 10 mg DAPA-HF exposure", targ_hi_mg, med_wt)
      }
    } else {
      if (is_nl) {
        sprintf("~%.2f – %.2f mg/gift (spreiding tussen Simcyp s1 en Farhan 2024 UGT1A9-ontogenie; plafond 10 mg)", targ_lo_mg, targ_hi_mg)
      } else {
        sprintf("~%.2f – %.2f mg/dose (spread between Simcyp s1 and Farhan 2024 UGT1A9 ontogeny; 10 mg cap)", targ_lo_mg, targ_hi_mg)
      }
    }

    aucu_ss_med <- stats::median(m$FU * m$AUC0_inf * (24 / st$tau), na.rm = TRUE)
    cur_dose_lbl <- if (st$dose_mode == "mg") sprintf("%.3g mg", st$dose_mg) else sprintf("%.3g mg/kg", st$dose_mgkg)

    if (pta$cur_above_safety > 10) {
      status_cls <- "hf-kpi-value hf-status-danger"
      status_txt <- sprintf(if (is_nl) "%.0f%% > veiligheidsplafond" else "%.0f%% > safety ceiling", pta$cur_above_safety)
    } else if (pta$cur_below > 50) {
      status_cls <- "hf-kpi-value hf-status-warn"
      status_txt <- sprintf(if (is_nl) "Subtherapeutisch (%.0f%% < ondergrens)" else "Subtherapeutic (%.0f%% < lower bound)", pta$cur_below)
    } else {
      status_cls <- "hf-kpi-value hf-status-ok"
      status_txt <- sprintf(if (is_nl) "PTA %.0f%% in streefvenster" else "PTA %.0f%% in target window", pta$cur_pta)
    }

    div(
      class = "hf-kpi-grid",
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label",
            if (is_nl) "1. Conservatieve startdosis (~5 mg-eq)" else "1. Conservative starting dose (~5 mg-eq)"),
        div(class = "hf-kpi-value", sprintf("%.2f mg/kg", ont$start_rec_mgkg)),
        div(class = "hf-kpi-sub",
            if (is_nl) {
              sprintf("~%.2f mg/gift — richt op AUCu,24h,ss ≈ 20 ng*h/mL (volwassen 5 mg-equivalent) onder het meest conservatieve UGT1A9-ontogeniemodel.", start_mg)
            } else {
              sprintf("~%.2f mg/dose — targets AUCu,24h,ss ≈ 20 ng*h/mL (adult 5 mg-equivalent) under the most conservative UGT1A9 ontogeny model.", start_mg)
            })
      ),
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label",
            if (is_nl) "2. Onderhouds-streefdosis (~10 mg DAPA-HF-eq)" else "2. Maintenance target dose (~10 mg DAPA-HF-eq)"),
        div(class = "hf-kpi-value", targ_val_txt),
        div(class = "hf-kpi-sub", targ_sub_txt)
      ),
      div(
        class = "hf-kpi-card",
        div(class = "hf-kpi-label",
            sprintf(if (is_nl) "3. Expositiebeoordeling huidige dosis (%s)" else "3. Exposure assessment of current dose (%s)", cur_dose_lbl)),
        div(class = status_cls, status_txt),
        div(class = "hf-kpi-sub",
            sprintf(if (is_nl) {
              "Mediane AUCu,24h,ss: %.1f ng*h/mL (streefvenster %g–%g | > %g ng*h/mL: %.0f%%)."
            } else {
              "Median AUCu,24h,ss: %.1f ng*h/mL (target window %g–%g | > %g ng*h/mL: %.0f%%)."
            }, aucu_ss_med, tg$target_lo, tg$target_hi, tg$safety_max, pta$cur_above_safety))
      )
    )
  })

  output$hf_alerts_ui <- renderUI({
    ont <- hf_ontogeny_data(); req(ont)
    if (length(ont$alerts) == 0L) return(NULL)
    is_nl <- (cur_lang() == "nl")
    div(
      class = "hf-alert-panel",
      strong(if (is_nl) {
        "Klinisch-farmacologische aandachtspunten voor het geselecteerde cohort:"
      } else {
        "Clinical-pharmacological considerations for the selected cohort:"
      }),
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
    build_window_explainer_plot(s, tg, ont, lang = cur_lang())
  })

  output$hf_window_cards_ui <- renderUI({
    tg <- hf_targets_reactive()
    if (cur_lang() == "nl") {
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
    } else {
      div(
        class = "hf-explain-grid-4",
        div(
          class = "hf-explain-card hf-zone-low",
          div(class = "hf-explain-title", sprintf("Target lower bound: %g ng*h/mL", tg$target_lo)),
          p(strong("Pharmacometric anchor: "), "Median exposure at ", strong("adult 5 mg once daily"),
            " (", code("19.2 ng*h/mL"), " at eGFR 90) and close to ", code("P5"), " of adults on 10 mg (", code("23.7 ng*h/mL"), "; Penland 2023 Table S4)."),
          p(strong("PD correlation: "), "Already yields ~79% mean 24-hour SGLT2 occupancy and ~80% of maximal TGF response; rational starting anchor for initiation.")
        ),
        div(
          class = "hf-explain-card hf-zone-ok",
          div(class = "hf-explain-title", "Maintenance midpoint: 38 – 48 ng*h/mL"),
          p(strong("Pharmacometric anchor: "), "Median exposure at ", strong("adult 10 mg once daily"),
            ": ", code("38.5 ng*h/mL"), " in normal renal function (eGFR 90) to ", code("48.2 ng*h/mL"), " for the typical ", strong("DAPA-HF"), " profile (66 y, eGFR 66)."),
          p(strong("Clinical interpretation: "), "Primary exposure target when up-titrating to paediatric maintenance dosing (optimisation anchor 42 ng*h/mL).")
        ),
        div(
          class = "hf-explain-card hf-zone-ok",
          div(class = "hf-explain-title", sprintf("Target upper bound: %g ng*h/mL", tg$target_hi)),
          p(strong("Pharmacometric anchor: "), "Corresponds to the ", strong("95th percentile (74.5 ng*h/mL)"),
            " in adults on 10 mg with mild-to-moderate renal impairment (eGFR 60–89; ", code("866 × 0.086"), ", Penland 2023 Table S4)."),
          p(strong("Clinical interpretation: "), sprintf("Exposures of %g–%g ng*h/mL lie within the expected inter-individual spread of the approved 10 mg HFrEF regimen.", tg$target_lo, tg$target_hi))
        ),
        div(
          class = "hf-explain-card hf-zone-max",
          div(class = "hf-explain-title", sprintf("Safety ceiling (P95): %g ng*h/mL", tg$safety_max)),
          p(strong("Pharmacometric anchor: "), "Matches the ", strong("median unbound exposure (87.7 ng*h/mL)"),
            " in severe chronic kidney disease (eGFR 15–29) on 10 mg/day (", code("1020 × 0.086"), "; Penland 2023 / DAPA-CKD)."),
          p(strong("Clinical interpretation: "), "Upper P95 cap for dose selection; above this threshold, long-term safety data without dose reduction are lacking.")
        )
      )
    }
  })

  output$hf_pta_ontogeny_plot <- renderPlot({
    s   <- sim_res(); req(s)
    pta <- hf_pta_data(); req(pta)
    ont <- hf_ontogeny_data(); req(ont)
    tg  <- hf_targets_reactive()
    st  <- s$res$settings
    is_nl <- (cur_lang() == "nl")

    gdf <- pta$grid_df
    dose_unit_lbl <- if (is_nl) {
      if (st$dose_mode == "mg") "Vaste keerdosis (mg)" else "Gewichtsgedoseerde keerdosis (mg/kg)"
    } else {
      if (st$dose_mode == "mg") "Fixed unit dose (mg)" else "Weight-based unit dose (mg/kg)"
    }
    mat_short_map <- if (is_nl) {
      c(
        farhan2024 = "Farhan 2024, in-vivo UGT1A9",
        badee2019  = "Badee 2019, PBPK s2",
        simcyp_s1  = "Simcyp s1, conservatief"
      )
    } else {
      c(
        farhan2024 = "Farhan 2024, in-vivo UGT1A9",
        badee2019  = "Badee 2019, PBPK s2",
        simcyp_s1  = "Simcyp s1, conservative"
      )
    }
    mat_short_lbl <- mat_short_map[[st$mat_option]]

    cat_in  <- sprintf(if (is_nl) "Binnen streefvenster (%g–%g ng*h/mL)" else "Within target window (%g–%g ng*h/mL)", tg$target_lo, tg$target_hi)
    cat_lo  <- sprintf(if (is_nl) "Subtherapeutisch (< %g ng*h/mL)" else "Subtherapeutic (< %g ng*h/mL)", tg$target_lo)
    cat_max <- sprintf(if (is_nl) "> Veiligheidsplafond (> %g ng*h/mL)" else "> Safety ceiling (> %g ng*h/mL)", tg$safety_max)

    pta_long <- rbind(
      data.frame(dose = gdf$dose, pct = gdf$pct_in_target,    cat = cat_in),
      data.frame(dose = gdf$dose, pct = gdf$pct_below,        cat = cat_lo),
      data.frame(dose = gdf$dose, pct = gdf$pct_above_safety, cat = cat_max)
    )
    lvls <- c(cat_in, cat_lo, cat_max)
    pta_long$cat <- factor(pta_long$cat, levels = lvls)
    cols_pta <- stats::setNames(c("#1a7f37", "#2c7fb8", "#b2182b"), lvls)
    ltys_pta <- stats::setNames(c("solid", "dotted", "dashed"), lvls)

    p_left <- ggplot(pta_long, aes(x = dose, y = pct, color = cat, linetype = cat)) +
      geom_line(linewidth = 1.15) +
      geom_vline(xintercept = pta$cur_dose, color = "#08306b", linewidth = 0.8) +
      geom_vline(xintercept = pta$dose_opt_pta, color = "#1a7f37", linetype = "dotdash", linewidth = 0.8) +
      annotate("text", x = pta$cur_dose, y = 97, hjust = -0.08, size = 3.2, fontface = "bold", color = "#08306b",
               label = sprintf(if (is_nl) "Huidig: %.3g" else "Current: %.3g", pta$cur_dose)) +
      scale_color_manual(values = cols_pta) +
      scale_linetype_manual(values = ltys_pta) +
      coord_cartesian(ylim = c(0, 103)) +
      labs(
        title = sprintf(if (is_nl) "A. Dosis–PTA curve (%s)" else "A. Dose–PTA curve (%s)", mat_short_lbl),
        subtitle = sprintf(
          if (is_nl) {
            "PTA-optimum = %.3g %s (PTA %.0f%% bij P95 <= %g; groene streeplijn)"
          } else {
            "PTA optimum = %.3g %s (PTA %.0f%% at P95 <= %g; green dash-dot line)"
          },
          pta$dose_opt_pta, if (st$dose_mode == "mg") "mg" else "mg/kg", pta$max_pta_pct, tg$safety_max
        ),
        x = dose_unit_lbl,
        y = if (is_nl) "Proportie virtueel cohort (%)" else "Proportion of virtual cohort (%)",
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
               color = "#b2182b", fontface = "bold",
               label = sprintf(if (is_nl) "Veiligheidsplafond (%g)" else "Safety ceiling (%g)", tg$safety_max)) +
      annotate("text", x = 3.45, y = tg$adult_10mg_hf, hjust = 1, vjust = -0.4, size = 3.0,
               color = "#1a7f37",
               label = if (is_nl) "DAPA-HF 10 mg mediaan (48,2)" else "DAPA-HF 10 mg median (48.2)") +
      labs(
        title = if (is_nl) "B. Gevoeligheid voor UGT1A9-ontogenie" else "B. Sensitivity to UGT1A9 ontogeny",
        subtitle = if (is_nl) {
          "Steady-state AUCu,24h bij huidige dosis (groene zone = streefvenster)"
        } else {
          "Steady-state AUCu,24h at current dose (green band = target window)"
        },
        x = NULL,
        y = if (is_nl) "Ongebonden blootstelling AUCu,24h,ss (ng*h/mL)" else "Unbound exposure AUCu,24h,ss (ng*h/mL)"
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
    lg <- cur_lang()
    is_nl <- (lg == "nl")

    cs_occ  <- s$res$cs_sglt2
    ref_occ <- adult_hf_reference_curves(tau = st$tau, out_by = 0.1, lang = lg)
    ref_lvls <- unique(ref_occ$ref)
    ref_occ$ref <- factor(ref_occ$ref, levels = ref_lvls)

    p_occ <- ggplot(cs_occ, aes(x = t, y = median)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2c7fb8", alpha = 0.22) +
      geom_line(color = "#08306b", linewidth = 1.2) +
      geom_line(data = ref_occ, aes(x = t, y = SGLT2_inh, linetype = ref, color = ref),
                inherit.aes = FALSE, linewidth = 0.95) +
      scale_color_manual(values = stats::setNames(c("#b2182b", "#d95f02"), ref_lvls)) +
      scale_linetype_manual(values = stats::setNames(c("dashed", "dotted"), ref_lvls)) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = if (is_nl) {
          "A. Proximaal tubulaire SGLT2-remming over het interval"
        } else {
          "A. Proximal tubular SGLT2 inhibition across the interval"
        },
        subtitle = if (is_nl) {
          "Blauw = cohort (mediaan + 90% PI); referenties = volwassen 10 & 5 mg QD"
        } else {
          "Blue = cohort (median + 90% PI); references = adult 10 & 5 mg QD"
        },
        x = if (is_nl) "Tijd na toediening (uur)" else "Time after dose (h)",
        y = if (is_nl) "Luminale SGLT2-receptorbezetting (%)" else "Luminal SGLT2 receptor occupancy (%)",
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
    lbl_tgf <- if (is_nl) {
      "Glomerulaire hemodynamiek / TGF Emax (EAUCu,50 = 5,2 ng*h/mL)"
    } else {
      "Glomerular haemodynamics / TGF Emax (EAUCu,50 = 5.2 ng*h/mL)"
    }
    lbl_hct <- if (is_nl) {
      "Systemisch / hematocriet Emax (EAUCu,50 = 20,7 ng*h/mL)"
    } else {
      "Systemic / haematocrit Emax (EAUCu,50 = 20.7 ng*h/mL)"
    }
    er_curves <- rbind(
      data.frame(aucu = x_seq, resp = 100 * x_seq / (th$eauc50_u_tgf + x_seq), as_type = lbl_tgf),
      data.frame(aucu = x_seq, resp = 100 * x_seq / (th$eauc50_u_hct + x_seq), as_type = lbl_hct)
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
        title = if (is_nl) {
          "B. Systemische Emax vs. GFR-gelimiteerde tubulaire PD"
        } else {
          "B. Systemic Emax vs. GFR-limited tubular PD"
        },
        subtitle = if (is_nl) {
          "Punten = tubulaire solute-respons (kleur = eGFR); lijnen = systemische Emax"
        } else {
          "Points = tubular solute response (colour = eGFR); curves = systemic Emax"
        },
        x = if (is_nl) "Ongebonden blootstelling AUCu,24h (ng*h/mL)" else "Unbound exposure AUCu,24h (ng*h/mL)",
        y = if (is_nl) "Relatieve PD-respons (%)" else "Relative PD response (%)",
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
    filename = function() {
      if (cur_lang() == "nl") {
        "dapagliflozine_pediatrisch_hartfalen_doseringstabel.csv"
      } else {
        "dapagliflozin_paediatric_hf_dosing_table.csv"
      }
    },
    content = function(file) {
      df <- hf_lifespan_df(); req(df)
      utils::write.csv(df, file, row.names = FALSE)
    }
  )

  # ---------- Expositiemetrieken ----------
  output$ss_table <- renderTable({
    s <- sim_res(); req(s)
    summ <- summarize_metrics(s$res$metrics, tau = s$res$settings$tau, include_pd = TRUE, lang = cur_lang())
    summ[] <- lapply(summ, function(x) if (is.numeric(x)) round(x, 2) else x)
    summ
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$metric_plots <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    is_nl <- (cur_lang() == "nl")
    facet_lvls <- if (is_nl) {
      c("Cmax totaal (ng/mL)", "AUCtau totaal (ng*h/mL)", "AUCu ongebonden (ng*h/mL)", "SGLT2-remming gem. (%)")
    } else {
      c("Cmax total (ng/mL)", "AUCtau total (ng*h/mL)", "AUCu unbound (ng*h/mL)", "Mean SGLT2 inhibition (%)")
    }
    met_long <- data.frame(
      metriek = factor(rep(facet_lvls, each = nrow(m)), levels = facet_lvls),
      waarde  = c(m$Cmax, m$AUCtau, m$AUCu, m$SGLT2_inh_avg)
    )
    ggplot(met_long, aes(x = metriek, y = waarde)) +
      geom_boxplot(fill = "#9ecae1", outlier.alpha = 0.4) +
      facet_wrap(~metriek, scales = "free", ncol = 4) +
      labs(
        x = NULL, y = NULL,
        title = if (is_nl) {
          "Verdeling van PK- en PD-metrieken (totaal, ongebonden & luminale SGLT2-remming)"
        } else {
          "Distribution of PK and PD metrics (total, unbound & luminal SGLT2 inhibition)"
        }
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank(), plot.title = element_text(face = "bold"))
  })

  sex_labels <- function(ages, lang = "en") {
    if (identical(lang, "nl")) {
      if (max(ages, na.rm = TRUE) >= 18) c("vrouw / meisje", "man / jongen") else c("meisje", "jongen")
    } else {
      if (max(ages, na.rm = TRUE) >= 18) c("female / girl", "male / boy") else c("girl", "boy")
    }
  }

  output$auc_age <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    lg <- cur_lang()
    is_nl <- (lg == "nl")
    lbls <- sex_labels(m$AGE, lg)
    m$SexGrp <- factor(m$SEXISMALE, levels = c(0, 1), labels = lbls)
    ggplot(m, aes(x = AGE, y = AUCu_24h, color = SexGrp)) +
      geom_hline(yintercept = c(20, 75), linetype = "dashed", color = "#1a7f37", linewidth = 0.6) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(
        x = if (is_nl) "Leeftijd (jaren)" else "Age (years)",
        y = if (is_nl) "Ongebonden AUCu,24h (ng*h/mL)" else "Unbound AUCu,24h (ng*h/mL)",
        color = if (is_nl) "Geslacht" else "Sex",
        title = if (is_nl) {
          "Ongebonden blootstelling (AUCu,24h) vs leeftijd (groene lijnen = 20-75 ng*h/mL referenties)"
        } else {
          "Unbound exposure (AUCu,24h) vs age (green dashed lines = 20-75 ng*h/mL references)"
        }
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$exposure_covariates <- renderPlot({
    s <- sim_res(); req(s)
    m <- s$res$metrics
    lg <- cur_lang()
    is_nl <- (lg == "nl")
    lbls <- sex_labels(m$AGE, lg)
    m$SexGrp <- factor(m$SEXISMALE, levels = c(0, 1), labels = lbls)

    m_aucu <- if (is_nl) "AUCu ongebonden (ng*h/mL)" else "Unbound AUCu (ng*h/mL)"
    m_rpd  <- if (is_nl) "GFR-gekoppelde renale PD (%)" else "GFR-coupled renal PD (%)"
    x_egfr <- "eGFR (mL/min/1.73m2)"
    x_age  <- if (is_nl) "Leeftijd (jaren)" else "Age (years)"

    long <- rbind(
      data.frame(metriek = m_aucu, xvariabele = x_egfr, x = m$EGFR, y = m$AUCu,         SexGrp = m$SexGrp),
      data.frame(metriek = m_aucu, xvariabele = x_age,  x = m$AGE,  y = m$AUCu,         SexGrp = m$SexGrp),
      data.frame(metriek = m_rpd,  xvariabele = x_egfr, x = m$EGFR, y = m$Renal_PD_eff, SexGrp = m$SexGrp),
      data.frame(metriek = m_rpd,  xvariabele = x_age,  x = m$AGE,  y = m$Renal_PD_eff, SexGrp = m$SexGrp)
    )
    long$metriek    <- factor(long$metriek, levels = c(m_aucu, m_rpd))
    long$xvariabele <- factor(long$xvariabele, levels = c(x_egfr, x_age))

    ggplot(long, aes(x = x, y = y, color = SexGrp)) +
      geom_point(alpha = 0.4, size = 1.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1) +
      facet_grid(metriek ~ xvariabele, scales = "free") +
      scale_color_manual(values = stats::setNames(c("#E69F00", "#0072B2"), lbls)) +
      labs(
        x = NULL, y = NULL,
        color = if (is_nl) "Geslacht" else "Sex",
        title = if (is_nl) {
          "Ongebonden PK (AUCu) én renale PD vs nierfunctie (eGFR) en leeftijd"
        } else {
          "Unbound PK (AUCu) and renal PD vs renal function (eGFR) and age"
        },
        subtitle = if (is_nl) {
          "Toont direct de tegengestelde richting bij lage GFR: AUCu stijgt, terwijl tubulaire osmotische PD daalt"
        } else {
          "Illustrates the opposing trajectories at low GFR: AUCu rises, whereas tubular osmotic PD declines"
        }
      ) +
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
    is_nl <- (cur_lang() == "nl")

    egfr_src <- if (is_nl && !is.null(st$egfr_source_nl)) st$egfr_source_nl
                else if (!is.null(st$egfr_source_en)) st$egfr_source_en
                else st$egfr_source
    egfr_unt <- if (is_nl && !is.null(st$egfr_unit_note_nl)) st$egfr_unit_note_nl
                else if (!is.null(st$egfr_unit_note_en)) st$egfr_unit_note_en
                else st$egfr_unit_note

    lines <- if (is_nl) {
      c(
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
        sprintf("eGFR-bron: %s", egfr_src),
        sprintf("eGFR-eenheid: mL/min/1.73m2 — %s", egfr_unt),
        if (!is.null(info)) sprintf("Dataset bevat serumcreatinine: %s",
                                    if (isTRUE(info$has_creatinine)) "ja" else "nee") else "",
        sprintf("Geslacht: %d%% mannelijk", round(100 * mean(pop$SEXISMALE)))
      )
    } else {
      c(
        sprintf("Dataset: %s", s$dataset),
        if (!is.null(info)) sprintf("Source: %s - %s (N=%s, age %g-%g years)",
                                    info$country, info$population,
                                    format(info$n_total, big.mark = ","),
                                    info$age_range[1], info$age_range[2]) else "",
        sprintf("Generated virtual population: N = %d (requested: %d)", nrow(pop), s$n_req),
        sprintf("Age: median %.2f years (range %.3f-%.2f)",
                stats::median(pop$AGE), min(pop$AGE), max(pop$AGE)),
        sprintf("Actual weight (BWT): median %.1f kg | Dry weight: median %.1f kg | Fat-free mass (FFM): median %.1f kg",
                stats::median(pop$WEIGHT), stats::median(m$BWT_DRY), stats::median(m$FFM)),
        sprintf("Serum albumin: median %.1f g/L | Unbound fraction fu: median %.3f (range %.3f-%.3f)",
                stats::median(m$ALB), stats::median(m$FU), min(m$FU), max(m$FU)),
        sprintf("Effective eGFR: median %.0f mL/min/1.73m2 (range %.0f-%.0f) | Relative renal function RF: median %.2f",
                stats::median(m$EGFR), min(m$EGFR), max(m$EGFR), stats::median(m$RF)),
        sprintf("Heart failure & enteral parameters: cardiorenal RF x%.2f | hepatic UGT1A9 x%.2f | oedema +%g%% | F1 %.0f%%",
                st$rf_factor, st$cl_hep_factor, st$edema_pct, 100 * st$f1_rel),
        sprintf("eGFR source: %s", egfr_src),
        sprintf("eGFR unit: mL/min/1.73m2 — %s", egfr_unt),
        if (!is.null(info)) sprintf("Dataset contains serum creatinine: %s",
                                    if (isTRUE(info$has_creatinine)) "yes" else "no") else "",
        sprintf("Sex: %d%% male", round(100 * mean(pop$SEXISMALE)))
      )
    }
    paste(lines[nzchar(lines)], collapse = "\n")
  })

  output$pop_plots <- renderPlot({
    s <- sim_res(); req(s)
    pop <- s$pop
    m   <- s$res$metrics
    is_nl <- (cur_lang() == "nl")
    v_lvls <- if (is_nl) {
      c("Leeftijd (jaren)", "Drooggewicht (kg)", "Vetvrije massa FFM (kg)", "Effectieve eGFR (mL/min/1.73m2)")
    } else {
      c("Age (years)", "Dry weight (kg)", "Fat-free mass FFM (kg)", "Effective eGFR (mL/min/1.73m2)")
    }
    pl <- rbind(
      data.frame(variabele = v_lvls[1], waarde = pop$AGE),
      data.frame(variabele = v_lvls[2], waarde = m$BWT_DRY),
      data.frame(variabele = v_lvls[3], waarde = m$FFM),
      data.frame(variabele = v_lvls[4], waarde = m$EGFR)
    )
    pl$variabele <- factor(pl$variabele, levels = v_lvls)
    ggplot(pl, aes(x = waarde)) +
      geom_histogram(bins = 30, fill = "#2c7fb8", color = "white") +
      facet_wrap(~variabele, scales = "free") +
      labs(
        x = NULL,
        y = if (is_nl) "Aantal" else "Count",
        title = if (is_nl) "Virtuele populatie gegenereerd met copulapop" else "Virtual population generated with copulapop"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$pop_table <- renderTable({
    s <- sim_res(); req(s)
    pop <- s$pop
    m   <- s$res$metrics
    is_nl <- (cur_lang() == "nl")
    out <- as.data.frame(lapply(pop, function(x) if (is.numeric(x)) round(x, 2) else x))
    if ("EGFR" %in% names(out)) out$EGFR <- NULL
    idx <- if ("ID" %in% names(out)) match(out$ID, m$ID) else seq_len(nrow(out))
    dry_col <- if (is_nl) "Drooggewicht (kg)" else "Dry weight (kg)"
    out[[dry_col]]              <- round(m$BWT_DRY[idx], 2)
    out[["FFM (kg)"]]           <- round(m$FFM[idx], 2)
    out[["Alb (g/L)"]]          <- round(m$ALB[idx], 1)
    out[["fu"]]                 <- round(m$FU[idx], 3)
    out[["EGFR (mL/min/1.73m2)"]] <- round(m$EGFR[idx], 1)
    head(out, 100)
  })

  # ---------- Downloads ----------
  output$dl_curve <- downloadHandler(
    filename = function() {
      if (cur_lang() == "nl") "dapagliflozine_pk_curve.csv" else "dapagliflozin_pk_curve.csv"
    },
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
    filename = function() {
      if (cur_lang() == "nl") "dapagliflozine_pkpd_metrics.csv" else "dapagliflozin_pkpd_metrics.csv"
    },
    content = function(file) {
      s <- sim_res(); req(s)
      utils::write.csv(s$res$metrics, file, row.names = FALSE)
    }
  )

  output$dl_nm_ctl <- downloadHandler(
    filename = "dapa_ped_hf_sim.ctl",
    content = function(file) {
      file.copy("nonmem/dapa_ped_hf_sim.ctl", file, overwrite = TRUE)
    }
  )

  output$dl_nm_csv <- downloadHandler(
    filename = "example_child_hf.csv",
    content = function(file) {
      file.copy("nonmem/example_child_hf.csv", file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
