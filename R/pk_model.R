# =====================================================================
# Fysiologisch geëxtrapoleerd populatie-PK model voor dapagliflozine
# Basis: NONMEM popPK-model van Penland et al. (J Clin Pharmacol 2023;
# 63:551-559, doi:10.1002/jcph.2196), uitgebreid voor fysiologisch
# plausibele extrapolatie over de gehele levensloop (preterm/neonaat tot
# oudere volwassene, obesitas en chronische nierschade):
#
#   1. Lichaamssamenstelling-allometrie via vetvrije massa (FFM; Al-Sallami
#      et al. 2015, doi:10.1007/s40262-015-0285-z; Janmahasatian et al. 2005)
#      en normal fat mass (NFM; Anderson & Holford 2014) voor CL, Q en V2.
#   2. Ontogenie van plasma-eiwitbinding aan albumine (McNamara & Alcorn
#      2002, doi:10.1208/ps040104; volwassen fu = 0.086) op CL/F, V3/F en
#      rapportage van ongebonden blootstelling (AUCu,ss).
#   3. Preterm-tot-volwassen nierfunctiematuratie en leeftijdsafname (Rhodin
#      et al. 2009, doi:10.1007/s00467-008-0997-5) + pediatrische CKiD U25
#      eGFR-formule (Pierce et al. 2021, doi:10.1016/j.kint.2020.10.047) voor
#      1-18 jaar en CKD-EPI 2021 voor >= 18 jaar.
#   4. Drie-paden klaringsdecompositie:
#        - Hepatische UGT1A9-glucuronidering: 5.44 * MAT * Fsize_CL
#        - Renale tubulaire UGT1A9-glucuronidering: 15.78 * MAT * Fsize_CL * RF
#          (met RF = eGFR / nGFR(AGE, BWT), zodat renale en hepatische UGT1A9
#          ontogenetisch matureren zonder dubbeltelling van GFR-rijping bij
#          gezonde zuigelingen, maar met behoud van CKD-afname bij lage eGFR)
#        - Onveranderde glomerulaire filtratie/excretie: 0.011 * eGFR * (BSA / 2)
#      Volledig geschaald met fu_ratio = fu(AGE) / fu_adult (laag-extractie).
#   5. Drie gevalideerde UGT1A9-ontogeniefuncties: Farhan et al. 2024 (in-vivo,
#      standaard), Badee et al. 2019 (Jo et al. 2021 PBPK scenario 2) en
#      Simcyp standaardontogenie (Jo et al. 2021 PBPK scenario 1).
#   6. Ontwikkelingsfysiologische GI-absorptie (neonatale maaglediging op KA;
#      gekalibreerde klinische Ka-BSV omega^2 = 0.36 als standaard naast het
#      bronmodel 8.52 FIX) + analytische ADVAN4 toestandspropagator.
#
# Structuur (ADVAN4 TRANS4): 3 compartimenten — depot (1) met lagtijd
# ALAG1, centraal (2, volume V2, eliminatie CL), perifeer (3, volume V3,
# uitwisseling Q). Eerste-orde absorptie KA van depot naar centraal.
#
# Eenheden: tijd in uur, volumes in L, BWT in kg, HEIGHT in m, AGE in jaar,
# EGFR in mL/min/1.73m2, concentraties in ng/mL (= microg/L).
# =====================================================================

nonmem_theta <- function() {
  list(
    ka         = 2.833,   # THETA(1)  volwassen KA (h^-1)
    cl_hep     = 5.44,    # hepatische UGT1A9 CL/F (L/h) bij ref-NFM, MAT=1, fu_ratio=1
    clnr       = 5.44,    # alias (niet-renale hepatische CL/F)
    cl_ren_ugt = 15.78,   # renale tubulaire UGT1A9 CL/F (L/h) bij ref-NFM, MAT=1, RF=1
    cl_filt_s  = 0.011,   # onveranderde renale excretie-slope (L/h per mL/min/1.73m2 bij BSA=2 m2)
    v2         = 73.47,   # THETA(3)  Vc/F (L) bij NFM_V2_ref = 61.68 kg
    v3         = 168.5,   # THETA(4)  Vp/F (L) bij BWT = 84 kg, fu_ratio = 1
    q          = 9.137,   # THETA(5)  Q/F (L/h) bij NFM_CL_ref = 63.08 kg
    bsa_ref    = 2.0,     # referentie-BSA (m2) bronpopulatie (84 kg, 1.71 m)
    bwt_ref    = 84.0,    # referentie-totaalgewicht (kg) Penland 2023
    ffm_ref    = 56.1,    # referentie vetvrije massa (kg) Penland 2023
    nfm_cl_ref = 63.08,   # FFM_ref + 0.25 * (84 - FFM_ref)
    nfm_v2_ref = 61.68,   # FFM_ref + 0.20 * (84 - FFM_ref)
    fu_adult   = 0.086,   # ongebonden fractie in volwassen plasma (91.4% gebonden)
    alb_adult  = 44.0,    # referentie plasma-albumine volwassene (g/L)
    sex        = -0.1441, # THETA(7) - niet toegepast (gedekt door geslachtsspecifieke FFM/NFM)
    alag       = 0.5,     # THETA(9)  ALAG1 (h), fixed
    gi_f0      = 0.45,    # relatieve maagledigingssnelheid bij geboorte t.o.v. volwassene
    gi_t50     = 0.10,    # leeftijd bij 50% GI-absorptiematuratie (jaar; ~5 weken)
    tm50            = 0.09,   # THETA(10) UGT1A9 maturatie TM50 (jaar), Farhan 2024 Table 1
    hill            = 1.9,    # THETA(11) UGT1A9 maturatie Hill-coeff, Farhan 2024
    fbirth          = 0.05,   # THETA(12) UGT1A9 Fbirth, Farhan 2024
    ic50_sglt2_ngml = 0.45,   # humaan SGLT2 in-vitro IC50 = 1.1 nM (MW 408.87 => 0.45 ng/mL; FDA NDA 202293)
    pt_conc_factor  = 2.25,   # proximale tubulus S1/S2 waterresorptie-concentratiefactor (Yakovleva 2014)
    cu50_sglt2_ngml = 0.20,   # effectieve vrije plasmaconcentratie voor 50% luminale SGLT2-remming (0.45 / 2.25)
    eauc50_u_tgf    = 5.2,    # vrije AUCu,50 voor acute GFR/TGF-respons (Kroonen 2020: 61 * 0.086 ng*h/mL)
    eauc50_u_hct    = 20.7    # vrije AUCu,50 voor hematocriet/volume-respons (Kroonen 2020: 241 * 0.086 ng*h/mL)
  )
}

# Maturatiefuncties voor hepatische + renale tubulaire UGT1A9, te kiezen in de app.
# Alle drie zijn geschaald op volwassen referentie = 1.0.
maturation_functions <- function(lang = "en") {
  is_nl <- identical(lang, "nl")
  list(
    farhan2024 = list(
      label = if (is_nl) {
        "Farhan et al. 2024 - doi:10.1002/jcph.2484 (in-vivo UGT1A9, standaard)"
      } else {
        "Farhan et al. 2024 - doi:10.1002/jcph.2484 (in-vivo UGT1A9, default)"
      },
      short = "Farhan 2024 (in-vivo)",
      fun = function(age, th = nonmem_theta()) {
        a <- pmax(as.numeric(age), 0)
        th$fbirth + (1 - th$fbirth) * (a^th$hill) / (th$tm50^th$hill + a^th$hill)
      }
    ),
    badee2019 = list(
      label = "Badee et al. 2019 via Jo et al. 2021 scenario 2 - doi:10.1002/psp4.12577",
      short = "Badee 2019 (Jo 2021 s2)",
      fun = function(age, th = NULL) {
        a <- pmax(as.numeric(age), 0)
        sig  <- 1.305 * (a^1.1) / (0.086^1.1 + a^1.1)
        expf <- pmax(1.0, 1.243 - 0.833 * exp(0.2 * (a - 25)))
        ifelse(a <= 0, 0, ifelse(a <= 1, sig, ifelse(a <= 20, expf, 1.0)))
      }
    ),
    simcyp_s1 = list(
      label = if (is_nl) {
        "Simcyp standaardontogenie via Jo et al. 2021 scenario 1 - doi:10.1002/psp4.12577"
      } else {
        "Simcyp default ontogeny via Jo et al. 2021 scenario 1 - doi:10.1002/psp4.12577"
      },
      short = "Simcyp s1 (Jo 2021 s1)",
      fun = function(age, th = NULL) {
        a <- pmax(as.numeric(age), 0)
        hill_s1 <- 0.086 + (1.041 - 0.086) * (a^1.2) / (1.16^1.2 + a^1.2)
        ifelse(a <= 15, hill_s1, 1.0)
      }
    )
  )
}

nonmem_omega <- function(ka_mode = c("calibrated", "penland_raw", "none")) {
  ka_mode <- match.arg(ka_mode)
  ka_var <- switch(
    ka_mode,
    calibrated  = 0.36, # SD(ETA1) = 0.60 (~66% CV); reproduceert klinische Cmax-CV (34-41%, Tirucherai 2016)
    penland_raw = 8.52, # Penland 2023 ruwe waarde (FIX in sparse dalspiegel-dataset)
    none        = 0.0
  )
  list(
    ka = ka_var,
    cl = 0.1136, # BSV op CL/F (ETA(2))
    v3 = 0.8229  # BSV op Vp/F (ETA(3)); geen BSV op V2 en Q in bronmodel
  )
}

# Fallback lengteschatting (m) alleen wanneer HEIGHT ontbreekt in synthetische dataframes
estimate_height_m <- function(bwt, age, sex_male = 1, bsa = NULL) {
  n <- length(bwt)
  if (!is.null(bsa) && length(bsa) == n && all(!is.na(bsa) & bsa > 0)) {
    hm <- (bsa^2 * 3600 / pmax(bwt, 0.5)) / 100
    return(pmin(pmax(hm, 0.35), 2.15))
  }
  a <- pmax(as.numeric(age), 0)
  h_adult <- ifelse(sex_male == 1, 1.76, 1.63)
  # Fysiologische groeicurve van 0.50 m bij geboorte naar volwassen lengte
  hm <- ifelse(
    a >= 18,
    h_adult,
    0.50 + (h_adult - 0.50) * ((a / 16)^0.78) / (1 + (a / 16)^0.78) * 1.88
  )
  hm <- pmin(hm, h_adult)
  # Zorg dat BMI bij extreme synthetische gewichten fysisch plausibel blijft (11-55 kg/m2)
  bmi_est <- bwt / (hm^2)
  hm <- ifelse(bmi_est > 55, sqrt(bwt / 45), ifelse(bmi_est < 11, sqrt(bwt / 13), hm))
  pmin(pmax(hm, 0.32), 2.15)
}

# Vetvrije massa (FFM, kg): Al-Sallami et al. 2015 (<18 j) + Janmahasatian et al. 2005 (>=18 j)
calc_ffm <- function(bwt, height_m, age, sex_male) {
  hm <- ifelse(height_m > 3, height_m / 100, height_m)
  bmi <- bwt / pmax(hm^2, 0.09)
  a <- pmax(as.numeric(age), 0)
  is_m <- (sex_male == 1)

  ffm_adult <- ifelse(
    is_m,
    9270 * bwt / (6680 + 216 * bmi),
    9270 * bwt / (8780 + 244 * bmi)
  )
  mat_m <- ifelse(a <= 0, 0.88, 0.88 + 0.12 / (1 + (a / 13.4)^(-12.7)))
  mat_f <- ifelse(a <= 0, 1.11, 1.11 - 0.11 / (1 + (a / 7.1)^(-1.1)))
  mat_ffm <- ifelse(a >= 18, 1.0, ifelse(is_m, mat_m, mat_f))

  pmin(ffm_adult * mat_ffm, 0.97 * bwt)
}

# Ontogenie plasma-albumine en ongebonden fractie fu (McNamara & Alcorn 2002)
fu_mcnamara <- function(age, alb_gL = NULL, fu_adult = 0.086, alb_adult = 44.0) {
  a <- pmax(as.numeric(age), 0)
  if (is.null(alb_gL) || all(is.na(alb_gL))) {
    days <- pmax(1.0, a * 365.25)
    alb <- ifelse(a >= 18, alb_adult, pmin(alb_adult, 1.1287 * log(days) + 33.746))
  } else {
    days <- pmax(1.0, a * 365.25)
    alb_def <- ifelse(a >= 18, alb_adult, pmin(alb_adult, 1.1287 * log(days) + 33.746))
    alb <- ifelse(is.na(alb_gL), alb_def, as.numeric(alb_gL))
  }
  ka_bind <- (1 - fu_adult) / fu_adult
  fu <- 1 / (1 + (alb / alb_adult) * ka_bind)
  list(fu = fu, fu_ratio = fu / fu_adult, alb = alb)
}

# Normale verwachte GFR (mL/min/1.73m2) voor leeftijd/gewicht:
# Rhodin et al. 2009 (PMA Hill-model) + fysiologische verouderingsafname boven 35 jaar
ngfr_rhodin <- function(age, bwt) {
  a <- pmax(as.numeric(age), 0)
  w <- pmax(as.numeric(bwt), 0.4)
  ga <- ifelse(a < 0.1 & w < 2.5, pmin(40, pmax(24, 22 + 7.2 * w)), 40)
  pma <- ga + 52.1775 * a
  115.0 * ((pma^3.4) / (47.7^3.4 + pma^3.4)) * (0.9925^pmax(0, a - 35))
}

# Pediatrische CKiD U25 eGFR (mL/min/1.73m2; Pierce et al. 2021, Kidney Int 99:948-956)
egfr_ckid_u25 <- function(creat_umol, age, sex_male, height_m) {
  scr_mgdl <- pmax(as.numeric(creat_umol) / 88.4, 0.05)
  hm <- ifelse(height_m > 3, height_m / 100, height_m)
  a <- as.numeric(age)
  is_m <- (sex_male == 1)
  k <- ifelse(
    is_m,
    ifelse(a < 12, 39.0 * (1.008^(a - 12)), 39.0 * (1.045^(a - 12))),
    ifelse(a < 12, 36.1 * (1.008^(a - 12)), 36.1 * (1.023^(a - 12)))
  )
  k * hm / scr_mgdl
}

# Volwassen CKD-EPI 2021 race-free eGFR (mL/min/1.73m2; Inker et al. 2021)
egfr_ckdepi_2021 <- function(creat_umol, age, sex_male) {
  scr_mgdl <- pmax(as.numeric(creat_umol) / 88.4, 0.05)
  is_m <- (sex_male == 1)
  kappa <- ifelse(is_m, 0.9, 0.7)
  alpha <- ifelse(is_m, -0.302, -0.241)
  sex_fac <- ifelse(is_m, 1.0, 1.012)
  cr_k <- scr_mgdl / kappa
  142 * (pmin(cr_k, 1)^alpha) * (pmax(cr_k, 1)^(-1.200)) * (0.9938^age) * sex_fac
}

# Individuele eGFR in mL/min/1.73m2 uit een (copulapop) populatie:
# - Als CREAT aanwezig is (bv. RUMC):
#     * >= 18 jaar : CKD-EPI 2021 (race-free)
#     * 1 - <18 jr : Pierce et al. 2021 CKiD U25 (leeftijds- en geslachtsafhankelijk)
#     * < 1 jaar   : Rhodin et al. 2009 PMA-GFR maturatie (cross-sectionele copula-CREAT
#                    mist maternale creatinineoverdracht en onderschat neonatale SCr)
#   Begrensd op [5, 180] mL/min/1.73m2 (Wang et al. 2020) tegen extreme spieratrofie-outliers.
# - Als CREAT ontbreekt maar EGFR expliciet is meegegeven (bv. deterministische validatie):
#     gebruik pop$EGFR direct.
# - Als zowel CREAT als EGFR ontbreken (bv. NHANES, HSE):
#     Rhodin et al. 2009 PMA-GFR maturatie + fysiologische BSV (SD = 0.215 op log-schaal).
egfr_from_copula <- function(pop, egfr_iiv = TRUE, egfr_cv_sd = 0.215, lang = "en") {
  is_nl <- identical(lang, "nl")
  n <- nrow(pop)
  age <- as.numeric(pop$AGE)
  bwt <- as.numeric(pop$WEIGHT)
  sex <- if ("SEXISMALE" %in% names(pop)) as.numeric(pop$SEXISMALE) else rep(1, n)
  bsa_in <- if ("BSA" %in% names(pop)) as.numeric(pop$BSA) else NULL
  hm <- if ("HEIGHT" %in% names(pop) && !all(is.na(pop$HEIGHT))) {
    h_raw <- as.numeric(pop$HEIGHT)
    ifelse(h_raw > 3, h_raw / 100, h_raw)
  } else {
    estimate_height_m(bwt, age, sex, bsa = bsa_in)
  }

  ngfr <- ngfr_rhodin(age, bwt)
  has_creat <- ("CREAT" %in% names(pop)) && any(!is.na(pop$CREAT))

  if (has_creat) {
    creat <- as.numeric(pop$CREAT)
    egfr <- rep(NA_real_, n)
    idx_adult <- which(!is.na(creat) & age >= 18)
    idx_ped   <- which(!is.na(creat) & age >= 1 & age < 18)
    idx_inf   <- which(is.na(creat) | age < 1)

    if (length(idx_adult) > 0) {
      egfr[idx_adult] <- egfr_ckdepi_2021(creat[idx_adult], age[idx_adult], sex[idx_adult])
    }
    if (length(idx_ped) > 0) {
      egfr[idx_ped] <- egfr_ckid_u25(creat[idx_ped], age[idx_ped], sex[idx_ped], hm[idx_ped])
    }
    if (length(idx_inf) > 0) {
      mult <- if (isTRUE(egfr_iiv)) exp(stats::rnorm(length(idx_inf), 0, egfr_cv_sd)) else 1.0
      egfr[idx_inf] <- ngfr[idx_inf] * mult
    }
    egfr <- pmin(pmax(egfr, 5), 180)
    parts_nl <- c(
      if (length(idx_inf) > 0) "Rhodin 2009 PMA (<1 j)" else NULL,
      if (length(idx_ped) > 0) "Pierce 2021 CKiD U25 (1-<18 j, uit CREAT)" else NULL,
      if (length(idx_adult) > 0) "CKD-EPI 2021 (>=18 j, uit CREAT)" else NULL
    )
    parts_en <- c(
      if (length(idx_inf) > 0) "Rhodin 2009 PMA (<1 y)" else NULL,
      if (length(idx_ped) > 0) "Pierce 2021 CKiD U25 (1-<18 y, from CREAT)" else NULL,
      if (length(idx_adult) > 0) "CKD-EPI 2021 (>=18 y, from CREAT)" else NULL
    )
    source_nl <- paste0("Dataset met CREAT: ", paste(parts_nl, collapse = " + "), " [begrensd 5-180 mL/min/1.73m2]")
    source_en <- paste0("Dataset with CREAT: ", paste(parts_en, collapse = " + "), " [bounded 5-180 mL/min/1.73m2]")
  } else if ("EGFR" %in% names(pop) && any(!is.na(pop$EGFR))) {
    egfr <- as.numeric(pop$EGFR)
    miss <- which(is.na(egfr))
    if (length(miss) > 0) {
      mult <- if (isTRUE(egfr_iiv)) exp(stats::rnorm(length(miss), 0, egfr_cv_sd)) else 1.0
      egfr[miss] <- pmin(pmax(ngfr[miss] * mult, 5), 180)
    }
    source_nl <- "Expliciet meegegeven EGFR-kolom (mL/min/1.73m2)"
    source_en <- "Explicitly supplied EGFR column (mL/min/1.73m2)"
  } else {
    mult <- if (isTRUE(egfr_iiv)) exp(stats::rnorm(n, 0, egfr_cv_sd)) else 1.0
    egfr <- pmin(pmax(ngfr * mult, 5), 180)
    source_nl <- paste0(
      "Rhodin et al. 2009 PMA-GFR maturatie + leeftijdsafname",
      if (isTRUE(egfr_iiv)) " incl. fysiologische GFR-BSV (22% CV)" else " (typische waarde)",
      " (dataset zonder creatinine)"
    )
    source_en <- paste0(
      "Rhodin et al. 2009 PMA-GFR maturation + age senescence",
      if (isTRUE(egfr_iiv)) " incl. physiological GFR-BSV (22% CV)" else " (typical value)",
      " (dataset without creatinine)"
    )
  }

  if (all(egfr >= 5 & egfr <= 180)) {
    unit_note_nl <- "sanity-check OK: alle waarden in fysiologisch bereik 5-180 mL/min/1.73m2"
    unit_note_en <- "sanity check OK: all values in physiological range 5-180 mL/min/1.73m2"
  } else if (all(egfr >= 1 & egfr <= 300)) {
    unit_note_nl <- "sanity-check OK: alle waarden in bereik 1-300 mL/min/1.73m2"
    unit_note_en <- "sanity check OK: all values in range 1-300 mL/min/1.73m2"
  } else {
    unit_note_nl <- "WAARSCHUWING: waarden buiten 1-300 - controleer de eenheid (mL/min/1.73m2 verwacht)"
    unit_note_en <- "WARNING: values outside 1-300 - check units (mL/min/1.73m2 expected)"
  }
  list(
    egfr         = egfr,
    ngfr         = ngfr,
    source       = if (is_nl) source_nl else source_en,
    source_en    = source_en,
    source_nl    = source_nl,
    unit_note    = if (is_nl) unit_note_nl else unit_note_en,
    unit_note_en = unit_note_en,
    unit_note_nl = unit_note_nl
  )
}

# Individuele PK-parameters uit de virtuele populatie (6-pijler fysiologisch model,
# uitbreidbaar met hartfalen-pathofysiologie: cardiorenaal RF, hepatische stuwing,
# extracellulair oedeem boven drooggewicht en hypoalbuminemie/PLE)
individual_params <- function(pop,
                              th = nonmem_theta(),
                              egfr = NULL,
                              mat_option = "farhan2024",
                              egfr_iiv = FALSE,
                              rf_factor = 1.0,
                              cl_hep_factor = 1.0,
                              edema_pct = 0.0,
                              alb_override = NULL) {
  mat_funs <- maturation_functions()
  if (!mat_option %in% names(mat_funs)) {
    stop("Onbekende maturatiefunctie: ", mat_option)
  }

  n   <- nrow(pop)
  bwt <- as.numeric(pop$WEIGHT)
  age <- pmax(as.numeric(pop$AGE), 0)
  sex <- if ("SEXISMALE" %in% names(pop)) as.numeric(pop$SEXISMALE) else rep(1, n)

  # Overvulling / extracellulair oedeem boven euvolemisch drooggewicht (PICU-definitie: W_obs = W_dry * (1 + FO/100))
  ed_pct  <- if (is.null(edema_pct) || is.na(edema_pct[1])) 0.0 else max(0.0, as.numeric(edema_pct[1]))
  bwt_dry <- bwt / (1 + ed_pct / 100)
  v_edema <- pmax(0.0, bwt - bwt_dry)

  bsa_in <- if ("BSA" %in% names(pop) && !all(is.na(pop$BSA))) as.numeric(pop$BSA) else NULL
  hm <- if ("HEIGHT" %in% names(pop) && !all(is.na(pop$HEIGHT))) {
    h_raw <- as.numeric(pop$HEIGHT)
    ifelse(h_raw > 3, h_raw / 100, h_raw)
  } else {
    estimate_height_m(bwt_dry, age, sex, bsa = bsa_in)
  }

  bsa <- if (!is.null(bsa_in)) bsa_in else sqrt((hm * 100 * bwt_dry) / 3600)
  bmi <- bwt_dry / pmax(hm^2, 0.09)

  # 1. Lichaamssamenstelling (FFM & NFM op euvolemisch drooggewicht: oedeemvocht heeft geen lever-/niermassa)
  ffm       <- calc_ffm(bwt_dry, hm, age, sex)
  nfm_cl    <- ffm + 0.25 * (bwt_dry - ffm)
  nfm_v2    <- ffm + 0.20 * (bwt_dry - ffm)
  f_size_cl <- (nfm_cl / th$nfm_cl_ref)^0.75
  f_size_v2 <- nfm_v2 / th$nfm_v2_ref

  # 2. Plasma-albumine binding ontogenie (fu) + optionele hypoalbuminemie/PLE override
  alb_in <- if (!is.null(alb_override) && length(alb_override) > 0 && !all(is.na(alb_override))) {
    rep(as.numeric(alb_override), length.out = n)
  } else if ("ALBUMIN" %in% names(pop)) {
    pop$ALBUMIN
  } else {
    NULL
  }
  fu_info  <- fu_mcnamara(age, alb_gL = alb_in, fu_adult = th$fu_adult, alb_adult = th$alb_adult)
  fu       <- fu_info$fu
  fu_ratio <- fu_info$fu_ratio

  # 3. Nierfunctie en relatieve renale functie t.o.v. gezonde leeftijdsnorm (+ cardiorenaal rf_factor)
  ngfr <- ngfr_rhodin(age, bwt_dry)
  if (is.null(egfr)) {
    egfr <- egfr_from_copula(pop, egfr_iiv = egfr_iiv)$egfr
  } else {
    egfr <- as.numeric(egfr)
  }
  rf_fac <- if (is.null(rf_factor) || is.na(rf_factor[1])) 1.0 else max(0.05, as.numeric(rf_factor[1]))
  if (abs(rf_fac - 1.0) > 1e-12) {
    egfr <- pmin(pmax(egfr * rf_fac, 5), 180)
  }
  rf <- pmin(pmax(egfr / pmax(ngfr, 1.0), 0.0), 1.75)

  # 4 & 5. Drie-paden klaring met hepatische + renale UGT1A9 maturatie (+ hepatische stuwing cl_hep_factor)
  hep_fac    <- if (is.null(cl_hep_factor) || is.na(cl_hep_factor[1])) 1.0 else max(0.05, as.numeric(cl_hep_factor[1]))
  mat        <- mat_funs[[mat_option]]$fun(age, th)
  cl_hep     <- fu_ratio * th$cl_hep * mat * f_size_cl * hep_fac
  cl_ren_ugt <- fu_ratio * th$cl_ren_ugt * mat * f_size_cl * rf
  cl_filt    <- fu_ratio * th$cl_filt_s * egfr * (bsa / th$bsa_ref)
  cl         <- cl_hep + cl_ren_ugt + cl_filt

  # 6. Verdelingsvolumes (+ extracellulair oedeemvolume in V2), Q en ontwikkelings-GI-absorptie
  v2   <- th$v2 * f_size_v2 + v_edema
  v3   <- th$v3 * (bwt_dry / th$bwt_ref) * fu_ratio
  q    <- th$q * f_size_cl
  f_gi <- th$gi_f0 + (1 - th$gi_f0) * (age / (th$gi_t50 + age))
  ka   <- th$ka * f_gi

  list(
    CL = cl, CL_hep = cl_hep, CL_ren_ugt = cl_ren_ugt, CL_filt = cl_filt,
    KA = ka, V2 = v2, V3 = v3, Q = q, ALAG = th$alag,
    BWT = bwt, BWT_DRY = bwt_dry, V_EDEMA = v_edema,
    HEIGHT = hm, BMI = bmi, FFM = ffm, NFM_CL = nfm_cl, NFM_V2 = nfm_v2,
    AGE = age, EGFR = egfr, NGFR = ngfr, RF = rf, BSA = bsa, SEX = sex,
    MAT = mat, FU = fu, FU_RATIO = fu_ratio, ALB = fu_info$alb,
    RF_FACTOR = rf_fac, CL_HEP_FACTOR = hep_fac, EDEMA_PCT = ed_pct
  )
}

# Exacte gevectoriseerde toestandspropagator Phi(dt) voor het lineaire
# 2-compartimentele orale absorptiemodel (ADVAN4):
#   dA1/dt = -ka * A1
#   dA2/dt =  ka * A1 - (ke + k12) * A2 + k21 * A3
#   dA3/dt =  k12 * A2 - k21 * A3
# Geeft de 7 niet-nul matrixcoëfficiënten terug voor tijdstap dt (scalair of vector >= 0).
.advan4_phi <- function(ka, alpha, beta, k12, k21, c_a, c_1, c_2, dt) {
  d_ab <- alpha - beta
  ea <- exp(-ka * dt)
  e1 <- exp(-alpha * dt)
  e2 <- exp(-beta * dt)

  p11 <- ea
  p21 <- ka * ((k21 - ka) * c_a * ea + (k21 - alpha) * c_1 * e1 + (k21 - beta) * c_2 * e2)
  p31 <- ka * k12 * (c_a * ea + c_1 * e1 + c_2 * e2)
  p22 <- ((alpha - k21) * e1 + (k21 - beta) * e2) / d_ab
  p32 <- k12 * (e2 - e1) / d_ab
  p23 <- k21 * (e2 - e1) / d_ab
  p33 <- ((k21 - beta) * e1 + (alpha - k21) * e2) / d_ab

  list(p11 = p11, p21 = p21, p31 = p31, p22 = p22, p32 = p32, p23 = p23, p33 = p33)
}

# Gevectoriseerde evaluatie van C2(t) = A2(t) / V2 voor een vector van tijdstippen `t_vec`
# binnen één doseringsinterval [0, tau], gegeven begintoestand (a1_0, a2_0, a3_0) op t = 0
# en een gelagde orale bolus `dose_in` die het depot binnenkomt op t = alag.
.eval_advan4_curve <- function(t_vec, a1_0, a2_0, a3_0, dose_in, alag,
                               v2, ka, alpha, beta, k12, k21, c_a, c_1, c_2) {
  n_subj <- length(v2)
  n_t    <- length(t_vec)
  d_ab   <- alpha - beta

  # Buitenproduct exp(-rate %o% t_vec): dimensie n_subj x n_t
  ea <- exp(-outer(ka, t_vec))
  e1 <- exp(-outer(alpha, t_vec))
  e2 <- exp(-outer(beta, t_vec))

  p21 <- (ka * (k21 - ka) * c_a) * ea +
         (ka * (k21 - alpha) * c_1) * e1 +
         (ka * (k21 - beta) * c_2) * e2
  p22 <- ((alpha - k21) / d_ab) * e1 + ((k21 - beta) / d_ab) * e2
  p23 <- (k21 / d_ab) * (e2 - e1)

  a2_mat <- a2_0 * p22 + a3_0 * p23 + a1_0 * p21

  idx_post <- which(t_vec > alag)
  if (length(idx_post) > 0) {
    tl <- t_vec[idx_post] - alag
    ea_l <- exp(-outer(ka, tl))
    e1_l <- exp(-outer(alpha, tl))
    e2_l <- exp(-outer(beta, tl))
    q21 <- (ka * (k21 - ka) * c_a) * ea_l +
           (ka * (k21 - alpha) * c_1) * e1_l +
           (ka * (k21 - beta) * c_2) * e2_l
    a2_mat[, idx_post] <- a2_mat[, idx_post] + dose_in * q21
  }
  a2_mat / v2
}

# Monte Carlo simulatie: n_doses doses met interval tau, per virtueel individu.
# Ondersteunt zowel gewichtgebaseerde dosering (dose_mgkg, zonder dosisplafond) als
# vaste mg-dosering (dose_mg), exacte analytische ADVAN4-oplossing over elk interval,
# optionele directe analytische steady state (steady_state = TRUE bij IOV = 0),
# en deterministische typische-patiëntsimulaties (iiv = FALSE, egfr_iiv = FALSE).
simulate_population_pk <- function(pop,
                                   dose_mgkg = NULL,
                                   dose_mg = NULL,
                                   tau = 24,
                                   n_doses = 7,
                                   iov_f1_cv = 0,
                                   ka_bsv = TRUE,
                                   ka_mode = "calibrated",
                                   mat_option = "farhan2024",
                                   iiv = TRUE,
                                   egfr_iiv = TRUE,
                                   steady_state = FALSE,
                                   seed = 1,
                                   th = nonmem_theta(),
                                   om = NULL,
                                   dt = NULL,
                                   out_by = 0.1,
                                   rf_factor = 1.0,
                                   cl_hep_factor = 1.0,
                                   edema_pct = 0.0,
                                   dose_weight_basis = c("dry", "actual"),
                                   alb_override = NULL,
                                   f1_rel = 1.0) {
  if (tau <= 0) stop("Dosisinterval (tau) moet groter dan 0 uur zijn.")
  if (is.null(dose_mgkg) && is.null(dose_mg)) {
    stop("Geef ofwel 'dose_mgkg' (mg/kg) ofwel 'dose_mg' (mg) op.")
  }
  if (!is.null(dose_mgkg) && !is.null(dose_mg)) {
    stop("Geef slechts één van 'dose_mgkg' of 'dose_mg' op, niet beide.")
  }
  dose_weight_basis <- match.arg(dose_weight_basis)

  set.seed(seed)
  n <- nrow(pop)

  if (is.null(om)) {
    eff_ka_mode <- if (!isTRUE(ka_bsv)) "none" else ka_mode
    om <- nonmem_omega(eff_ka_mode)
  } else if (!isTRUE(ka_bsv)) {
    om$ka <- 0
  }

  egfr_info <- egfr_from_copula(pop, egfr_iiv = isTRUE(egfr_iiv))
  p <- individual_params(
    pop,
    th            = th,
    egfr          = egfr_info$egfr,
    mat_option    = mat_option,
    rf_factor     = rf_factor,
    cl_hep_factor = cl_hep_factor,
    edema_pct     = edema_pct,
    alb_override  = alb_override
  )

  # Inter-individuele variabiliteit (log-normaal, exp(ETA))
  if (isTRUE(iiv)) {
    if (om$cl > 0) p$CL <- p$CL * exp(stats::rnorm(n, 0, sqrt(om$cl)))
    if (om$v3 > 0) p$V3 <- p$V3 * exp(stats::rnorm(n, 0, sqrt(om$v3)))
    if (om$ka > 0) {
      p$KA <- p$KA * exp(stats::rnorm(n, 0, sqrt(om$ka)))
      p$KA <- pmin(pmax(p$KA, 0.02), 50)
    }
  }

  wt_for_dose <- if (dose_weight_basis == "actual") p$BWT else p$BWT_DRY
  dose_mg_vec <- if (!is.null(dose_mg)) {
    rep(as.numeric(dose_mg), length.out = n)
  } else {
    as.numeric(dose_mgkg) * wt_for_dose
  }
  if (any(is.na(dose_mg_vec) | dose_mg_vec <= 0)) {
    stop("Alle individuele doses moeten positief zijn.")
  }
  f1_rel_val <- if (is.null(f1_rel) || is.na(f1_rel[1])) 1.0 else max(0.05, as.numeric(f1_rel[1]))
  dose_ug <- dose_mg_vec * 1000 * f1_rel_val  # mg -> microg (zodat microg / L = ng/mL)

  # Mean-gecentreerde log-normale IOV op F1: E[F1] = 1.0 bij elke CV%
  sigma_iov <- sqrt(log(1 + (iov_f1_cv / 100)^2))

  # Hybride snelheidsconstanten per individu
  ke  <- p$CL / p$V2
  k12 <- p$Q / p$V2
  k21 <- p$Q / p$V3
  s_k <- ke + k12 + k21
  disc <- sqrt(pmax(0, s_k^2 - 4 * ke * k21))
  alpha <- 0.5 * (s_k + disc)
  beta  <- 0.5 * (s_k - disc)

  # Bescherm tegen numerieke singulariteit wanneer ka exact samenvalt met alpha of beta
  ka <- p$KA
  near_a <- abs(ka - alpha) < 1e-5
  if (any(near_a)) ka[near_a] <- ka[near_a] + 2e-5
  near_b <- abs(ka - beta) < 1e-5
  if (any(near_b)) ka[near_b] <- ka[near_b] + 2e-5
  p$KA <- ka

  c_a <- 1 / ((alpha - ka) * (beta - ka))
  c_1 <- 1 / ((ka - alpha) * (beta - alpha))
  c_2 <- 1 / ((ka - beta) * (alpha - beta))

  alag <- th$alag
  if (tau <= alag) {
    stop(sprintf("Dosisinterval tau (%g h) moet groter zijn dan de absorptielagtijd ALAG1 (%g h).", tau, alag))
  }

  times <- seq(0, tau, by = out_by)
  if (tail(times, 1) < tau - 1e-12) times <- c(times, tau)
  n_out <- length(times)

  phi_tau  <- .advan4_phi(ka, alpha, beta, k12, k21, c_a, c_1, c_2, tau)
  phi_post <- .advan4_phi(ka, alpha, beta, k12, k21, c_a, c_1, c_2, tau - alag)

  # Dag-1 curve (eerste dosis vanaf lege compartimenten)
  f1_first <- if (sigma_iov > 0) exp(stats::rnorm(n, -0.5 * sigma_iov^2, sigma_iov)) else rep(1, n)
  d_first  <- dose_ug * f1_first
  curve_first <- .eval_advan4_curve(
    times, numeric(n), numeric(n), numeric(n), d_first, alag,
    p$V2, ka, alpha, beta, k12, k21, c_a, c_1, c_2
  )

  if (isTRUE(steady_state) && sigma_iov == 0) {
    # Exacte gesloten-vorm periodieke steady state: A_ss(0) = (I - Phi(tau))^-1 * (dose_ug * phi_1(tau - alag))
    b1 <- dose_ug * phi_post$p11
    b2 <- dose_ug * phi_post$p21
    b3 <- dose_ug * phi_post$p31

    a1_0 <- b1 / (1 - phi_tau$p11)
    r2 <- b2 + phi_tau$p21 * a1_0
    r3 <- b3 + phi_tau$p31 * a1_0
    m11 <- 1 - phi_tau$p22
    m12 <- -phi_tau$p23
    m21 <- -phi_tau$p32
    m22 <- 1 - phi_tau$p33
    det_m <- m11 * m22 - m12 * m21
    a2_0 <- (m22 * r2 - m12 * r3) / det_m
    a3_0 <- (-m21 * r2 + m11 * r3) / det_m
    d_last <- dose_ug
  } else {
    n_d <- max(1L, as.integer(round(n_doses)))
    A1 <- numeric(n)
    A2 <- numeric(n)
    A3 <- numeric(n)
    if (n_d == 1L) {
      a1_0 <- A1; a2_0 <- A2; a3_0 <- A3
      d_last <- d_first
    } else {
      # Propageer dosis 1 over [0, tau]
      A1 <- d_first * phi_post$p11
      A2 <- d_first * phi_post$p21
      A3 <- d_first * phi_post$p31
      if (n_d > 2L) {
        for (d in 2:(n_d - 1L)) {
          f1_d <- if (sigma_iov > 0) exp(stats::rnorm(n, -0.5 * sigma_iov^2, sigma_iov)) else 1.0
          dd   <- dose_ug * f1_d
          nA1 <- phi_tau$p11 * A1 + dd * phi_post$p11
          nA2 <- phi_tau$p21 * A1 + phi_tau$p22 * A2 + phi_tau$p23 * A3 + dd * phi_post$p21
          nA3 <- phi_tau$p31 * A1 + phi_tau$p32 * A2 + phi_tau$p33 * A3 + dd * phi_post$p31
          A1 <- nA1; A2 <- nA2; A3 <- nA3
        }
      }
      a1_0 <- A1; a2_0 <- A2; a3_0 <- A3
      f1_last <- if (sigma_iov > 0) exp(stats::rnorm(n, -0.5 * sigma_iov^2, sigma_iov)) else rep(1, n)
      d_last  <- dose_ug * f1_last
    }
  }

  # Evalueer laatste doseringsinterval op het uitvoerrooster `times`
  curve_last <- .eval_advan4_curve(
    times, a1_0, a2_0, a3_0, d_last, alag,
    p$V2, ka, alpha, beta, k12, k21, c_a, c_1, c_2
  )

  # Eindtoestand op t = tau voor exacte analytische massabalans-AUCtau
  a1_end <- phi_tau$p11 * a1_0 + d_last * phi_post$p11
  a2_end <- phi_tau$p21 * a1_0 + phi_tau$p22 * a2_0 + phi_tau$p23 * a3_0 + d_last * phi_post$p21
  a3_end <- phi_tau$p31 * a1_0 + phi_tau$p32 * a2_0 + phi_tau$p33 * a3_0 + d_last * phi_post$p31
  auc_tau <- ((a1_0 + a2_0 + a3_0 + d_last) - (a1_end + a2_end + a3_end)) / p$CL

  # Verfijn Cmax en Tmax lokaal (+/- out_by rond het rooster-maximum op 0.01 h resolutie; 40x sneller)
  imax_g <- max.col(curve_last, ties.method = "first")
  imax_g[is.na(imax_g)] <- 1L
  if (out_by > 0.01) {
    n_ref   <- max(3L, as.integer(round(2 * out_by / 0.01)) + 1L)
    offsets <- seq(-out_by, out_by, length.out = n_ref)
    t_loc   <- pmin(pmax(outer(times[imax_g], offsets, "+"), alag), tau)
    d_ab    <- alpha - beta

    ea <- exp(-ka * t_loc)
    e1 <- exp(-alpha * t_loc)
    e2 <- exp(-beta * t_loc)
    p21 <- (ka * (k21 - ka) * c_a) * ea +
           (ka * (k21 - alpha) * c_1) * e1 +
           (ka * (k21 - beta) * c_2) * e2
    p22 <- ((alpha - k21) / d_ab) * e1 + ((k21 - beta) / d_ab) * e2
    p23 <- (k21 / d_ab) * (e2 - e1)
    a2_loc <- a2_0 * p22 + a3_0 * p23 + a1_0 * p21

    tl <- pmax(t_loc - alag, 0)
    post_mask <- (t_loc > alag)
    ea_l <- exp(-ka * tl)
    e1_l <- exp(-alpha * tl)
    e2_l <- exp(-beta * tl)
    q21 <- (ka * (k21 - ka) * c_a) * ea_l +
           (ka * (k21 - alpha) * c_1) * e1_l +
           (ka * (k21 - beta) * c_2) * e2_l
    c_loc <- (a2_loc + post_mask * (d_last * q21)) / p$V2

    imax_l   <- max.col(c_loc, ties.method = "first")
    idx_pair <- cbind(seq_len(n), imax_l)
    cmax_val <- c_loc[idx_pair]
    tmax_val <- t_loc[idx_pair]
  } else {
    cmax_val <- curve_last[cbind(seq_len(n), imax_g)]
    tmax_val <- times[imax_g]
  }

  t_half_beta <- log(2) / beta
  auc_inf_sd  <- dose_ug / p$CL

  # Ongebonden concentratieprofiel Cu(t) = fu * C(t) & semi-mechanistische SGLT2/cardiorenaal PD
  curve_ss_u   <- p$FU * curve_last
  curve_day1_u <- p$FU * curve_first
  cu50         <- if (!is.null(th$cu50_sglt2_ngml)) th$cu50_sglt2_ngml else 0.20
  occ_mat      <- 100 * curve_ss_u / (cu50 + curve_ss_u)
  if (n_out >= 2L) {
    dt_vec    <- diff(times)
    sglt2_avg <- as.vector(((occ_mat[, -n_out, drop = FALSE] + occ_mat[, -1L, drop = FALSE]) %*% dt_vec) / (2 * tau))
  } else {
    sglt2_avg <- occ_mat[, 1L]
  }
  sglt2_min    <- occ_mat[, n_out]
  aucu_tau     <- p$FU * auc_tau
  aucu_24h     <- aucu_tau * (24 / tau)
  renal_pd_eff <- sglt2_avg * pmin(1.5, p$EGFR / 100.0)
  eauc50_tgf   <- if (!is.null(th$eauc50_u_tgf)) th$eauc50_u_tgf else 5.2
  eauc50_hct   <- if (!is.null(th$eauc50_u_hct)) th$eauc50_u_hct else 20.7
  e_tgf_pct    <- 100 * aucu_24h / (eauc50_tgf + aucu_24h)
  e_hct_pct    <- 100 * aucu_24h / (eauc50_hct + aucu_24h)

  metrics <- data.frame(
    ID            = if ("ID" %in% names(pop)) pop$ID else seq_len(n),
    AGE           = p$AGE,
    BWT           = p$BWT,
    BWT_DRY       = p$BWT_DRY,
    HEIGHT        = p$HEIGHT,
    BMI           = p$BMI,
    FFM           = p$FFM,
    EGFR          = p$EGFR,
    NGFR          = p$NGFR,
    RF            = p$RF,
    BSA           = p$BSA,
    SEXISMALE     = p$SEX,
    MAT           = p$MAT,
    ALB           = p$ALB,
    FU            = p$FU,
    CL            = p$CL,
    KA            = p$KA,
    V2            = p$V2,
    V3            = p$V3,
    Q             = p$Q,
    t_half        = t_half_beta,
    DOSE_MG       = dose_mg_vec,
    Cmax          = cmax_val,
    Tmax          = tmax_val,
    Cmin          = curve_last[, n_out],
    Cu_max        = p$FU * cmax_val,
    Cu_min        = curve_ss_u[, n_out],
    AUCtau        = auc_tau,
    AUC0_inf      = auc_inf_sd,
    Cavg          = auc_tau / tau,
    AUCu          = aucu_tau,
    AUCu_24h      = aucu_24h,
    Cavgu         = aucu_tau / tau,
    SGLT2_inh_avg = sglt2_avg,
    SGLT2_inh_min = sglt2_min,
    Renal_PD_eff  = renal_pd_eff,
    E_TGF_pct     = e_tgf_pct,
    E_HCT_pct     = e_hct_pct
  )

  ss_attainment <- stats::median(auc_tau / auc_inf_sd, na.rm = TRUE)

  egfr_src_en <- egfr_info$source_en
  egfr_src_nl <- egfr_info$source_nl
  if (abs(p$RF_FACTOR - 1.0) > 1e-12) {
    egfr_src_en <- sprintf("%s [cardiorenally scaled x %.2f]", egfr_src_en, p$RF_FACTOR)
    egfr_src_nl <- sprintf("%s [cardiorenaal geschaald x %.2f]", egfr_src_nl, p$RF_FACTOR)
  }

  list(
    pop          = pop,
    params       = p,
    times        = times,
    curve_ss     = curve_last,
    curve_day1   = curve_first,
    curve_ss_u   = curve_ss_u,
    curve_day1_u = curve_day1_u,
    occ_sglt2    = occ_mat,
    cs_ss        = curve_summary(curve_last, times),
    cs_day1      = curve_summary(curve_first, times),
    cs_ss_u      = curve_summary(curve_ss_u, times),
    cs_day1_u    = curve_summary(curve_day1_u, times),
    cs_sglt2     = curve_summary(occ_mat, times),
    metrics      = metrics,
    settings     = list(
      dose_mgkg         = if (!is.null(dose_mgkg)) as.numeric(dose_mgkg) else NA_real_,
      dose_mg           = if (!is.null(dose_mg)) as.numeric(dose_mg)[1] else NA_real_,
      dose_mode         = if (!is.null(dose_mg)) "mg" else "mgkg",
      dose_weight_basis = dose_weight_basis,
      tau               = tau,
      n_doses           = n_doses,
      steady_state      = steady_state,
      ss_attainment     = ss_attainment,
      iov_f1_cv         = iov_f1_cv,
      f1_rel            = f1_rel_val,
      rf_factor         = p$RF_FACTOR,
      cl_hep_factor     = p$CL_HEP_FACTOR,
      edema_pct         = p$EDEMA_PCT,
      alb_override      = if (!is.null(alb_override) && !all(is.na(alb_override))) as.numeric(alb_override)[1] else NA_real_,
      ka_bsv            = ka_bsv,
      ka_mode           = if (!isTRUE(ka_bsv)) "none" else ka_mode,
      seed              = seed,
      n                 = n,
      mat_option        = mat_option,
      mat_label         = maturation_functions()[[mat_option]]$short,
      egfr_source       = egfr_src_en,
      egfr_source_en    = egfr_src_en,
      egfr_source_nl    = egfr_src_nl,
      egfr_unit_note    = egfr_info$unit_note_en,
      egfr_unit_note_en = egfr_info$unit_note_en,
      egfr_unit_note_nl = egfr_info$unit_note_nl
    )
  )
}

# Populatie-quantielen per tijdspunt (exacte R type-7 quantielen via C-level quickselect)
curve_summary <- function(curve, times, probs = c(0.05, 0.5, 0.95)) {
  nr <- nrow(curve)
  if (nr == 0L) {
    return(data.frame(t = times, q05 = NA_real_, median = NA_real_, q95 = NA_real_))
  }
  if (nr == 1L || anyNA(curve)) {
    q <- apply(curve, 2, stats::quantile, probs = probs, names = FALSE, na.rm = TRUE)
    if (is.null(dim(q))) q <- matrix(q, nrow = length(probs))
    return(data.frame(t = times, q05 = q[1, ], median = q[2, ], q95 = q[3, ]))
  }
  idx <- 1 + (nr - 1) * probs
  lo  <- floor(idx)
  hi  <- ceiling(idx)
  h   <- idx - lo
  part_idx <- sort(unique(c(lo, hi)))
  q <- vapply(seq_len(ncol(curve)), function(j) {
    s <- sort.int(curve[, j], partial = part_idx)
    (1 - h) * s[lo] + h * s[hi]
  }, numeric(length(probs)))
  data.frame(t = times, q05 = q[1, ], median = q[2, ], q95 = q[3, ])
}

# Samenvatting steady-state metrieken: mediaan, geometrisch gemiddelde + P5/P95
summarize_metrics <- function(metrics, tau = 24, include_pd = FALSE, lang = "en") {
  is_nl <- identical(lang, "nl")
  cols <- c("Cmax", "Tmax", "Cmin", "AUCtau", "AUCu", "Cavg", "t_half")
  labels <- if (is_nl) {
    if (isTRUE(include_pd)) {
      c(
        "Cmax,ss totaal (ng/mL)",
        "Tmax (uur)",
        sprintf("Cmin,ss = C%.0f totaal (ng/mL)", tau),
        "AUCtau,ss totaal (ng*h/mL)",
        "AUCu,ss ongebonden (ng*h/mL)",
        "Cavg,ss totaal (ng/mL)",
        "Terminale t1/2,beta (uur)"
      )
    } else {
      c(
        "Cmax,ss (ng/mL)",
        "Tmax (uur)",
        sprintf("Cmin,ss = C%.0f (ng/mL)", tau),
        "AUCtau,ss totaal (ng*h/mL)",
        "AUCu,ss ongebonden (ng*h/mL)",
        "Cavg,ss totaal (ng/mL)",
        "Terminale t1/2,beta (uur)"
      )
    }
  } else {
    if (isTRUE(include_pd)) {
      c(
        "Cmax,ss total (ng/mL)",
        "Tmax (h)",
        sprintf("Cmin,ss = C%.0f total (ng/mL)", tau),
        "AUCtau,ss total (ng*h/mL)",
        "AUCu,ss unbound (ng*h/mL)",
        "Cavg,ss total (ng/mL)",
        "Terminal t1/2,beta (h)"
      )
    } else {
      c(
        "Cmax,ss (ng/mL)",
        "Tmax (h)",
        sprintf("Cmin,ss = C%.0f (ng/mL)", tau),
        "AUCtau,ss total (ng*h/mL)",
        "AUCu,ss unbound (ng*h/mL)",
        "Cavg,ss total (ng/mL)",
        "Terminal t1/2,beta (h)"
      )
    }
  }
  if (isTRUE(include_pd)) {
    cols <- c(cols, "Cu_max", "Cu_min", "AUCu_24h", "SGLT2_inh_avg", "SGLT2_inh_min", "Renal_PD_eff", "E_TGF_pct", "E_HCT_pct")
    labels <- c(
      labels,
      if (is_nl) {
        c(
          "Cu,max,ss ongebonden (ng/mL)",
          sprintf("Cu,min,ss = Cu%.0f ongebonden (ng/mL)", tau),
          "AUCu,24u-eq ongebonden (ng*h/mL)",
          "Luminale SGLT2-remming gem. over interval (%)",
          "Luminale SGLT2-remming op dalspiegel (%)",
          "GFR-gekoppelde renale PD-index (% van rijpe nefroncap.)",
          "Cardiorenaal TGF-effect (% van Emax; Kroonen 2020)",
          "Systemisch/HCT-effect (% van Emax; Kroonen 2020)"
        )
      } else {
        c(
          "Cu,max,ss unbound (ng/mL)",
          sprintf("Cu,min,ss = Cu%.0f unbound (ng/mL)", tau),
          "AUCu,24h-eq unbound (ng*h/mL)",
          "Luminal SGLT2 inhibition mean over interval (%)",
          "Luminal SGLT2 inhibition at trough (%)",
          "GFR-coupled renal PD index (% of mature nephron cap.)",
          "Cardiorenal TGF effect (% of Emax; Kroonen 2020)",
          "Systemic/HCT effect (% of Emax; Kroonen 2020)"
        )
      }
    )
  }
  present <- cols %in% names(metrics)
  cols <- cols[present]
  labels <- labels[present]

  geomean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) > 0 && all(x > 0)) exp(mean(log(x))) else NA_real_
  }
  out <- data.frame(
    Metric       = labels,
    Median       = sapply(metrics[cols], stats::median, na.rm = TRUE),
    `Geom. mean` = sapply(metrics[cols], geomean),
    P5           = sapply(metrics[cols], function(x) stats::quantile(x, 0.05, names = FALSE, na.rm = TRUE)),
    P95          = sapply(metrics[cols], function(x) stats::quantile(x, 0.95, names = FALSE, na.rm = TRUE)),
    check.names  = FALSE
  )
  if (is_nl) {
    names(out) <- c("Metriek", "Mediaan", "Geom. gem.", "P5", "P95")
  }
  out
}
