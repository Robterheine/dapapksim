#!/usr/bin/env Rscript
# =============================================================================
# Genereert het 1-kind voorbeeldbestand (nonmem/example_child_hf.csv) en de
# exacte R-referentietabel (nonmem/expected_r_reference.csv) vanuit dapapksim
# (R/pk_model.R), én verifieert automatisch de NONMEM-uitvoertabel
# (nonmem/sdtab_dapa_ped_hf.tab) zodra NONMEM lokaal is gedraaid.
#
# Gebruik vanuit de projectroot (of vanuit de map nonmem/):
#   Rscript nonmem/generate_and_verify.R
# =============================================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else if (dir.exists("nonmem")) {
  normalizePath("nonmem")
} else {
  normalizePath(".")
}
repo_root <- if (file.exists(file.path(script_dir, "..", "R", "pk_model.R"))) {
  normalizePath(file.path(script_dir, ".."))
} else {
  normalizePath(".")
}

source(file.path(repo_root, "R", "pk_model.R"))

# -----------------------------------------------------------------------------
# 1. Klinisch voorbeeldprofiel: zuigeling van 9 maanden (0.75 jaar) met
#    stabiel tot matig pediatrisch hartfalen:
#      - Actueel gewicht BWT = 8.40 kg bij +5% overvulling => drooggewicht 8.00 kg
#      - Lengte HT = 0.71 m, jongen (SEX = 1)
#      - Cardiorenaal GFR-behoud RFFAC = 0.80, leverstuwing HEPFAC = 0.90
#      - Hypoalbuminemie ALBIN = 32.0 g/L, enterale F1REL = 1.00
#      - UGT1A9-ontogenie: MATOPT = 1 (Farhan et al. 2024 in-vivo)
#      - Dosering: 0.15 mg/kg op drooggewicht = 1.20 mg 1x per 24 uur (7 giften)
# -----------------------------------------------------------------------------
child_pop <- data.frame(
  ID        = 1L,
  AGE       = 0.75,
  WEIGHT    = 8.40,
  HEIGHT    = 0.71,
  SEXISMALE = 1L
)

edema_pct     <- 5.0
rf_factor     <- 0.80
cl_hep_factor <- 0.90
alb_val       <- 32.0
f1_rel        <- 1.00
mat_opt_code  <- 1L
mat_opt_name  <- "farhan2024"
dose_mgkg     <- 0.15
bwt_dry       <- child_pop$WEIGHT / (1 + edema_pct / 100) # exact 8.00 kg
dose_mg       <- round(dose_mgkg * bwt_dry, 6)            # exact 1.20 mg
tau_h         <- 24.0
n_doses       <- 7L

tad_grid <- c(0.00, 0.25, 0.50, 0.60, 0.75, 1.00, 1.50, 2.00,
              3.00, 4.00, 6.00, 8.00, 12.00, 16.00, 20.00, 24.00)

# Doseringsrecord op TIME = 0 met ADDL = 6, II = 24 (7 giften op t = 0, 24, ..., 144 h)
dose_row <- data.frame(
  ID     = 1L,
  TIME   = 0.00,
  TAD    = 0.00,
  DAY    = 1L,
  EVID   = 1L,
  AMT    = dose_mg,
  ADDL   = n_doses - 1L,
  II     = tau_h,
  CMT    = 1L,
  DV     = 0,
  MDV    = 1L,
  AGE    = child_pop$AGE,
  BWT    = child_pop$WEIGHT,
  HT     = child_pop$HEIGHT,
  SEX    = child_pop$SEXISMALE,
  CREAT  = -99,
  EGFRIN = -99,
  ALBIN  = alb_val,
  EDEMA  = edema_pct,
  RFFAC  = rf_factor,
  HEPFAC = cl_hep_factor,
  F1REL  = f1_rel,
  MATOPT = mat_opt_code,
  DOSEMG = dose_mg,
  TAU    = tau_h
)

make_obs_rows <- function(day_num, t_offset) {
  data.frame(
    ID     = 1L,
    TIME   = t_offset + tad_grid,
    TAD    = tad_grid,
    DAY    = as.integer(day_num),
    EVID   = 0L,
    AMT    = 0,
    ADDL   = 0L,
    II     = 0,
    CMT    = 2L,
    DV     = 0,
    MDV    = 1L,
    AGE    = child_pop$AGE,
    BWT    = child_pop$WEIGHT,
    HT     = child_pop$HEIGHT,
    SEX    = child_pop$SEXISMALE,
    CREAT  = -99,
    EGFRIN = -99,
    ALBIN  = alb_val,
    EDEMA  = edema_pct,
    RFFAC  = rf_factor,
    HEPFAC = cl_hep_factor,
    F1REL  = f1_rel,
    MATOPT = mat_opt_code,
    DOSEMG = dose_mg,
    TAU    = tau_h
  )
}

nm_data <- rbind(
  dose_row,
  make_obs_rows(day_num = 1L, t_offset = 0.0),
  make_obs_rows(day_num = 7L, t_offset = 144.0)
)

csv_path <- file.path(script_dir, "example_child_hf.csv")
utils::write.csv(nm_data, csv_path, row.names = FALSE, quote = FALSE)
cat("Geschreven NONMEM dataset:", csv_path, sprintf("(%d rijen)\n", nrow(nm_data)))

# -----------------------------------------------------------------------------
# 2. Exacte R-referentie vanuit dapapksim (simulate_population_pk op 0.05 h grid)
# -----------------------------------------------------------------------------
sim_r <- simulate_population_pk(
  pop               = child_pop,
  dose_mg           = dose_mg,
  tau               = tau_h,
  n_doses           = n_doses,
  iiv               = FALSE,
  egfr_iiv          = FALSE,
  steady_state      = FALSE,
  mat_option        = mat_opt_name,
  out_by            = 0.05,
  rf_factor         = rf_factor,
  cl_hep_factor     = cl_hep_factor,
  edema_pct         = edema_pct,
  dose_weight_basis = "dry",
  alb_override      = alb_val,
  f1_rel            = f1_rel
)

pk_par <- individual_params(
  pop           = child_pop,
  th            = nonmem_theta(),
  egfr          = ngfr_rhodin(child_pop$AGE, bwt_dry),
  mat_option    = mat_opt_name,
  rf_factor     = rf_factor,
  cl_hep_factor = cl_hep_factor,
  edema_pct     = edema_pct,
  alb_override  = alb_val
)

idx_match <- vapply(tad_grid, function(tt) which.min(abs(sim_r$times - tt)), integer(1))
stopifnot(max(abs(sim_r$times[idx_match] - tad_grid)) < 1e-12)

ctot_d1 <- as.numeric(sim_r$curve_day1[1L, idx_match])
ctot_d7 <- as.numeric(sim_r$curve_ss[1L, idx_match])

obs_nm  <- nm_data[nm_data$EVID == 0L, , drop = FALSE]
ctot_v  <- c(ctot_d1, ctot_d7)
cu_v    <- pk_par$FU * ctot_v
sglt2_v <- 100.0 * cu_v / (0.20 + cu_v)
aucu24  <- pk_par$FU * ((dose_mg * 1000.0 * f1_rel) / pk_par$CL) * (24.0 / tau_h)
etgf_v  <- 100.0 * aucu24 / (5.2 + aucu24)
ehct_v  <- 100.0 * aucu24 / (20.7 + aucu24)

ref_df <- data.frame(
  ID     = obs_nm$ID,
  TIME   = obs_nm$TIME,
  TAD    = obs_nm$TAD,
  DAY    = obs_nm$DAY,
  WDRY   = round(pk_par$BWT_DRY, 6),
  FFM    = round(pk_par$FFM, 6),
  ALBU   = round(pk_par$ALB, 6),
  FU     = round(pk_par$FU, 6),
  NGFR   = round(pk_par$NGFR, 6),
  EGFR   = round(pk_par$EGFR, 6),
  RF     = round(pk_par$RF, 6),
  MAT    = round(pk_par$MAT, 6),
  CL     = round(pk_par$CL, 6),
  CLHEP  = round(pk_par$CL_hep, 6),
  CLRUGT = round(pk_par$CL_ren_ugt, 6),
  CLFILT = round(pk_par$CL_filt, 6),
  V2     = round(pk_par$V2, 6),
  V3     = round(pk_par$V3, 6),
  Q      = round(pk_par$Q, 6),
  KA     = round(pk_par$KA, 6),
  CTOT   = round(ctot_v, 6),
  CU     = round(cu_v, 6),
  SGLT2  = round(sglt2_v, 6),
  AUCU24 = round(aucu24, 6),
  ETGF   = round(etgf_v, 6),
  EHCT   = round(ehct_v, 6)
)

ref_path <- file.path(script_dir, "expected_r_reference.csv")
utils::write.csv(ref_df, ref_path, row.names = FALSE, quote = FALSE)
cat("Geschreven R-referentietabel:", ref_path, "\n\n")

cat("=== Samenvatting individuele PK/PD-parameters voor ID = 1 ===\n")
print(ref_df[1L, c("WDRY", "FFM", "ALBU", "FU", "NGFR", "EGFR", "RF", "MAT",
                   "CL", "CLHEP", "CLRUGT", "CLFILT", "V2", "V3", "Q", "KA",
                   "AUCU24", "ETGF", "EHCT")], row.names = FALSE)

cat("\n=== Verwachte concentraties & SGLT2-bezetting (selectie Dag 1 vs. Dag 7) ===\n")
sel <- ref_df$TAD %in% c(0, 0.5, 1, 2, 6, 12, 24)
print(ref_df[sel, c("DAY", "TIME", "TAD", "CTOT", "CU", "SGLT2")], row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Optionele automatische verificatie indien NONMEM sdtab aanwezig is
# -----------------------------------------------------------------------------
sdtab_path <- file.path(script_dir, "sdtab_dapa_ped_hf.tab")
if (file.exists(sdtab_path)) {
  cat("\n=== Verificatie van gevonden NONMEM-output (sdtab_dapa_ped_hf.tab) ===\n")
  nm_out <- utils::read.table(sdtab_path, header = TRUE, skip = 1)
  nm_obs <- nm_out[nm_out$EVID == 0, , drop = FALSE]
  stopifnot(nrow(nm_obs) == nrow(ref_df))
  check_cols <- c("CL", "V2", "V3", "Q", "KA", "FU", "EGFR", "CTOT", "CU", "SGLT2", "AUCU24")
  for (cc in check_cols) {
    max_abs <- max(abs(as.numeric(nm_obs[[cc]]) - as.numeric(ref_df[[cc]])))
    cat(sprintf("  %-8s : max |NONMEM - R| = %.3e\n", cc, max_abs))
  }
}
