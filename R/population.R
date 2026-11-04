# =====================================================================
# Wrapper rond het copulapop package (Robterheine/copulapop)
# Genereert een virtuele patientenpopulatie met vine-copulas.
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
      copulapop::set_data_path(candidates[1])
    } else {
      stop("Geen copulapop JSON-databestanden gevonden in de package-installatie. ",
           "Geef het pad op in het veld 'copulapop data-pad' of gebruik copulapop::set_data_path().")
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

# Genereer n kinderen in de gevraagde leeftijdsrange.
# copulapop::generate_population() kan minder dan n leveren bij een smalle
# leeftijdsrange (interne oversampling is beperkt), daarom wordt stapsgewijs
# met een grotere request-n gewerkt todat n kinderen beschikbaar zijn.
generate_virtual_population <- function(dataset, n, age_range, sex = "both", seed = 1) {
  ensure_copulapop()
  best <- NULL
  for (k in c(1, 2, 4, 8, 16, 32, 64)) {
    pop <- suppressWarnings(copulapop::generate_population(
      dataset, n = min(n * k, 200000), age_range = age_range, sex = sex,
      seed = seed))
    if (is.null(pop) || nrow(pop) == 0) next
    if (nrow(pop) >= n) {
      set.seed(as.integer((as.numeric(seed) + 1) %% 2147483647))
      return(pop[sample(nrow(pop), n), , drop = FALSE])
    }
    if (is.null(best) || nrow(pop) > nrow(best)) best <- pop
  }
  if (is.null(best)) stop("copulapop kon geen populatie genereren voor deze instellingen.")
  warning(sprintf("Slechts %d van %d gevraagde kinderen gegenereerd in de leeftijdsrange %g-%g jaar.",
                  nrow(best), n, age_range[1], age_range[2]))
  best
}
