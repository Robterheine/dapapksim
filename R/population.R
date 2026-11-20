# =====================================================================
# Wrapper rond het copulapop package (Robterheine/copulapop)
# Genereert een virtuele patiëntenpopulatie met vine-copulas over elke
# gewenste leeftijdsrange (neonaten t/m volwassenen).
# =====================================================================

ensure_copulapop <- function() {
  if (!requireNamespace("copulapop", quietly = TRUE)) {
    stop("Package 'copulapop' is niet geinstalleerd. Installeer met:\n",
         "  install.packages(c(\"VineCopula\", \"data.table\", \"jsonlite\", \"Matrix\", \"MASS\"))\n",
         "  devtools::install_github(\"Robterheine/copulapop\")")
  }
  if (is.null(copulapop::get_data_path())) {
    candidates <- c(
      system.file("extdata", package = "copulapop"),
      system.file("data", package = "copulapop")
    )
    candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
    has_json <- sapply(candidates, function(p) length(list.files(p, pattern = "\\.json$")) > 0)
    candidates <- candidates[has_json]
    if (length(candidates) > 0) {
      suppressMessages(copulapop::set_data_path(candidates[1]))
    } else {
      stop("Geen copulapop JSON-databestanden gevonden in de package-installatie. ",
           "Gebruik copulapop::set_data_path() om het pad naar de JSON-bestanden in te stellen.")
    }
  }
  invisible(TRUE)
}

copulapop_dataset_info <- function(dataset) {
  tryCatch({
    ensure_copulapop()
    copulapop::get_dataset_info(dataset)
  }, error = function(e) NULL)
}

# Genereer n virtuele individuen in elke gewenste leeftijdsrange (van pasgeborene tot oudere).
# Bevat twee belangrijke optimalisaties rond copulapop::generate_population():
#   1. Vermijd de has_adult-drempel op exact 18.0 jaar (waarbij copulapop anders 93.5% van de
#      trekkingen in het volwassen stratum 18-100 j plaatst en daarna wegfiltert).
#   2. Adaptieve twee-staps oversampling met accumulatie over onafhankelijke seeds, zodat ook
#      zeer smalle neonatale/infantiele leeftijdsvensters (bv. 0-0.1 jaar) snel N halen.
generate_virtual_population <- function(dataset, n, age_range, sex = "both", seed = 1) {
  ensure_copulapop()
  n <- max(1L, as.integer(round(as.numeric(n))))
  if (length(age_range) != 2 || any(is.na(age_range))) {
    stop("Ongeldige leeftijdsrange: geef een minimum- en maximumleeftijd op.")
  }
  lo <- min(as.numeric(age_range))
  hi <- max(as.numeric(age_range))

  info <- copulapop_dataset_info(dataset)
  if (!is.null(info) && !is.null(info$age_range)) {
    ds_lo <- as.numeric(info$age_range[1])
    ds_hi <- as.numeric(info$age_range[2])
    if (hi < ds_lo || lo > ds_hi) {
      stop(sprintf(
        "Dataset %s bevat alleen leeftijden van %g tot %g jaar (gevraagd: %g-%g jaar). Kies een passende dataset of leeftijdsrange.",
        dataset, ds_lo, ds_hi, lo, hi
      ))
    }
    lo <- max(lo, ds_lo)
    hi <- min(hi, ds_hi)
    # Als onder- en bovengrens samenvallen (bv. exact 0.5-0.5 jaar), verbreed minimaal (+/-0.05 j)
    if (hi - lo < 1e-4) {
      lo <- max(ds_lo, lo - 0.05)
      hi <- min(ds_hi, hi + 0.05)
    }
  } else if (hi - lo < 1e-4) {
    lo <- max(0, lo - 0.05)
    hi <- hi + 0.05
  }

  # Vermijd copulapop's has_adult = (age_range[2] >= 18) sprong bij puur pediatrische selectie t/m 18 j
  query_hi <- if (lo < 18 && abs(hi - 18) < 1e-6) 17.9999 else hi
  query_range <- c(lo, query_hi)

  # Schat benodigde eerste steekproefgrootte op basis van relatieve breedte van het leeftijdsvenster
  span <- max(query_hi - lo, 0.02)
  ds_max <- if (!is.null(info)) as.numeric(info$age_range[2]) else 100
  stratum_span <- if (query_hi < 18) {
    18 - (if (!is.null(info)) as.numeric(info$age_range[1]) else 0)
  } else if (lo >= 18) {
    max(1, ds_max - 18)
  } else {
    # Gemengde pediatrisch + volwassen trekking: houd rekening met laag pediatrisch aandeel in RUMC
    ped_frac <- if (!is.null(info) && !is.null(info$n_pediatric) && !is.null(info$n_total) && info$n_total > 0) {
      max(0.05, info$n_pediatric / info$n_total)
    } else 0.25
    ds_max / min(1, ped_frac * 2.5)
  }

  est_mult <- max(1.0, (stratum_span / span) * 1.15)
  n_first  <- as.integer(min(200000, max(n, ceiling(n * est_mult))))

  base_seed <- as.integer(abs(as.numeric(seed)) %% 2147400000) + 1L
  acc <- NULL
  total_requested <- 0
  first_err <- NULL

  for (batch in 1:4) {
    if (batch == 1L) {
      n_req <- n_first
    } else {
      got <- if (is.null(acc)) 0L else nrow(acc)
      rem <- n - got
      if (rem <= 0L) break
      emp_rate <- max(got / max(total_requested, 1), 1e-4)
      n_req <- as.integer(min(200000, max(ceiling(rem / emp_rate * 1.35), n * 4L)))
    }
    batch_seed <- as.integer((base_seed + (batch - 1L) * 104729L) %% 2147483647L)
    chunk <- tryCatch(
      suppressWarnings(copulapop::generate_population(
        dataset, n = n_req, age_range = query_range, sex = sex, seed = batch_seed
      )),
      error = function(e) {
        if (is.null(first_err)) first_err <<- e
        NULL
      }
    )
    total_requested <- total_requested + n_req
    if (!is.null(chunk) && nrow(chunk) > 0) {
      acc <- if (is.null(acc)) chunk else rbind(acc, chunk)
      if (nrow(acc) >= n) break
    }
  }

  if (is.null(acc) || nrow(acc) == 0) {
    if (!is.null(first_err)) {
      stop("copulapop::generate_population faalde: ", conditionMessage(first_err))
    }
    stop("copulapop kon geen populatie genereren voor deze instellingen.")
  }

  if (nrow(acc) >= n) {
    set.seed(as.integer((base_seed + 1L) %% 2147483647L))
    out <- acc[sample(nrow(acc), n), , drop = FALSE]
    out$ID <- seq_len(nrow(out))
    rownames(out) <- NULL
    return(out)
  }

  warning(sprintf(
    "Slechts %d van %d gevraagde individuen gegenereerd in de leeftijdsrange %g-%g jaar.",
    nrow(acc), n, lo, hi
  ))
  acc$ID <- seq_len(nrow(acc))
  rownames(acc) <- NULL
  acc
}
