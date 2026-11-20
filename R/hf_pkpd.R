# =====================================================================
# Pediatrisch Hartfalen PK/PD & Dosisvoorspellingsmodule voor dapapksim
#
# Koppelt het fysiologische 6-pijler popPK-model (R/pk_model.R) aan:
#   1. Hartfalen-specifieke pathofysiologie presets (cardiorenaal syndroom,
#      veneuze leverstuwing / Fontan, extracellulair oedeem en hypoalbuminemie/PLE).
#   2. Semi-mechanistische duale PK-PD modellering:
#        - Systemische ongebonden blootstelling (AUCu,24h) gekoppeld aan
#          klinische cardiorenaal Emax-relaties (Kroonen et al. 2020,
#          Br J Clin Pharmacol 86:1905-1916, doi:10.1111/bcp.14318).
#        - Proximaal tubulaire luminale SGLT2-receptorbezetting (IC50 = 1.1 nM
#          = 0.45 ng/mL, proximale waterresorptie-factor 2.25x; Yakovleva et al.
#          2014, Front Pharmacol 5:218) én GFR-gekoppelde osmotische/natriuretische
#          tubulaire PD-capaciteit (gefiltreerde solute-vracht ~ eGFR / 100).
#   3. Analytische Probability of Target Attainment (PTA) dosis-optimalisatie
#      en UGT1A9-ontogenie gevoeligheidsanalyse (Farhan 2024 vs Badee 2019 vs
#      Simcyp s1) voor veilige pediatrische "start low, go slow" titratie.
# =====================================================================

hf_presets <- function() {
  list(
    standard = list(
      label         = "Standaard fysiologie (geen extra hartfalen-correctie)",
      rf_factor     = 1.00,
      cl_hep_factor = 1.00,
      edema_pct     = 0,
      alb_mode      = "auto",
      alb_val       = 40,
      f1_rel        = 1.00,
      note          = "Gezonde leeftijdsnorm voor nierdoorbloeding, leverklaring en serumalbumine."
    ),
    hf_stable = list(
      label         = "Stabiel chronisch pediatrisch HF (milde perfusiereductie, RF 80%)",
      rf_factor     = 0.80,
      cl_hep_factor = 0.90,
      edema_pct     = 0,
      alb_mode      = "auto",
      alb_val       = 38,
      f1_rel        = 1.00,
      note          = "Weerspiegelt stabiel behandeld HFrEF/DCM onder achtergrondtherapie (ACEi/ARNI + diuretica; ~20% hogere AUC conform volwassen DAPA-HF popPK)."
    ),
    hf_cardiorenal = list(
      label         = "Matig-ernstig cardiorenaal HF (RF 55%, leverstuwing 70%, +10% oedeem)",
      rf_factor     = 0.55,
      cl_hep_factor = 0.70,
      edema_pct     = 10,
      alb_mode      = "manual",
      alb_val       = 30,
      f1_rel        = 0.90,
      note          = "Low cardiac output + veneuze stuwing + 10% extracellulaire overvulling + milde hypoalbuminemie (30 g/L) en splanchnische malabsorptie (F1 = 90%)."
    ),
    hf_fontan_ple = list(
      label         = "Fontan / rechtszijdig HF met PLE (leverstuwing 65%, albumine 24 g/L)",
      rf_factor     = 0.75,
      cl_hep_factor = 0.65,
      edema_pct     = 5,
      alb_mode      = "manual",
      alb_val       = 24,
      f1_rel        = 0.85,
      note          = "Verhoogde centraal-veneuze druk (Fontan-leverstuwing) + uitgesproken hypoalbuminemie (24 g/L) door protein-losing enteropathy: totale spiegels dalen sterk terwijl vrije AUCu hoog blijft!"
    ),
    custom = list(
      label         = "Handmatig aangepast hartfalen- / toedieningsprofiel",
      rf_factor     = NULL,
      cl_hep_factor = NULL,
      edema_pct     = NULL,
      alb_mode      = NULL,
      alb_val       = NULL,
      f1_rel        = NULL,
      note          = "Vrij ingestelde pathofysiologische en farmaceutische parameters."
    )
  )
}

# Volwassen klinische referentiewaarden voor ongebonden 24-uurs blootstelling (AUCu,24h in ng*h/mL)
hf_reference_targets <- function(target_lo = 20, target_hi = 75, safety_max = 90) {
  list(
    adult_5mg_ref    = 19.2,  # 5 mg QD bij 84 kg, eGFR 90 (AUCtau ~224 * 0.086)
    adult_10mg_norm  = 38.5,  # 10 mg QD bij 84 kg, eGFR 90 (AUCtau ~447 * 0.086; Penland ref)
    adult_10mg_hf    = 48.2,  # 10 mg QD typische DAPA-HF patiënt (66 j, 80 kg, eGFR 66: ~560 * 0.086)
    adult_ckd_severe = 87.7,  # 10 mg QD ernstige CKD eGFR 15-29 mediaan (Penland Tabel S4: 1020 * 0.086)
    target_lo        = target_lo,
    target_hi        = target_hi,
    safety_max       = safety_max
  )
}

# Typische volwassen DAPA-HF referentiecurves (10 mg QD en 5 mg QD steady state) voor SGLT2- en Cu-vergelijking
adult_hf_reference_curves <- function(tau = 24, out_by = 0.1) {
  pop_ad <- rbind(
    data.frame(AGE = 66, WEIGHT = 82, HEIGHT = 1.75, SEXISMALE = 1, EGFR = 66),
    data.frame(AGE = 66, WEIGHT = 72, HEIGHT = 1.63, SEXISMALE = 0, EGFR = 66)
  )
  # Gebruik 24u QD als klinische volwassen standaard, geëvalueerd over [0, max(tau, 24)] en afgekapt op tau
  tau_eval <- max(24, tau)
  r10 <- simulate_population_pk(pop_ad, dose_mg = 10, tau = 24, steady_state = TRUE,
                                iiv = FALSE, egfr_iiv = FALSE, out_by = out_by)
  r05 <- simulate_population_pk(pop_ad, dose_mg = 5,  tau = 24, steady_state = TRUE,
                                iiv = FALSE, egfr_iiv = FALSE, out_by = out_by)
  df10 <- data.frame(
    t         = r10$times,
    Cu        = colMeans(r10$curve_ss_u),
    SGLT2_inh = colMeans(r10$occ_sglt2),
    ref       = "Volwassen DAPA-HF ref: 10 mg 1x/dag"
  )
  df05 <- data.frame(
    t         = r05$times,
    Cu        = colMeans(r05$curve_ss_u),
    SGLT2_inh = colMeans(r05$occ_sglt2),
    ref       = "Volwassen startdosis ref: 5 mg 1x/dag"
  )
  out <- rbind(df10, df05)
  out[out$t <= tau + 1e-9, , drop = FALSE]
}

# Analytische Probability of Target Attainment (PTA) scan over een dosisreeks op
# chronische steady state. Omdat voor onderhoudstherapie bij hartfalen de steady-state
# blootstelling bepalend is (ook als de gebruiker in de grafiek naar dosis 1 kijkt),
# gebruiken we exacte individuele steady-state AUCu,24h,ss = FU * AUC0_inf * (24 / tau).
compute_pta_scan <- function(res,
                             target_lo = 20,
                             target_hi = 75,
                             safety_max = 90,
                             n_grid = 150) {
  m  <- res$metrics
  st <- res$settings
  cur_dose    <- if (st$dose_mode == "mg") st$dose_mg else st$dose_mgkg
  aucu_ss_24h <- m$FU * m$AUC0_inf * (24 / st$tau)
  unit_aucu   <- aucu_ss_24h / cur_dose
  unit_aucu   <- unit_aucu[is.finite(unit_aucu) & unit_aucu > 0]
  if (length(unit_aucu) == 0L) return(NULL)

  med_u <- stats::median(unit_aucu)
  p95_u <- stats::quantile(unit_aucu, 0.95, names = FALSE)

  # Mijlpaaldoses op basis van de gesimuleerde populatieverdeling
  d_5mg_raw    <- 20.0 / med_u
  d_10mg_raw   <- 42.0 / med_u   # midden tussen gezonde 10 mg (38.5) en DAPA-HF (48.2)
  d_safety_cap <- safety_max / p95_u

  d_start_rec  <- min(d_5mg_raw, d_safety_cap)
  d_target_rec <- min(d_10mg_raw, d_safety_cap)

  # Adaptief dosisrooster dat zowel de huidige dosis als het optimale gebied bestrijkt
  max_grid <- if (st$dose_mode == "mg") {
    max(12.0, min(25.0, max(cur_dose * 1.6, d_10mg_raw * 2.0)))
  } else {
    max(0.35, min(1.0, max(cur_dose * 1.6, d_10mg_raw * 2.0)))
  }
  min_grid <- if (st$dose_mode == "mg") 0.05 else 0.005
  grid <- sort(unique(c(seq(min_grid, max_grid, length.out = n_grid), cur_dose, d_start_rec, d_target_rec)))

  # Gevectoriseerde evaluatie: matrix van (length(unit_aucu) x length(grid))
  mat_aucu <- outer(unit_aucu, grid)
  pct_below        <- 100 * colMeans(mat_aucu < target_lo)
  pct_in_target    <- 100 * colMeans(mat_aucu >= target_lo & mat_aucu <= target_hi)
  pct_above_target <- 100 * colMeans(mat_aucu > target_hi)
  pct_above_safety <- 100 * colMeans(mat_aucu > safety_max)
  med_vec          <- apply(mat_aucu, 2, stats::median)
  p05_vec          <- apply(mat_aucu, 2, stats::quantile, probs = 0.05, names = FALSE)
  p95_vec          <- apply(mat_aucu, 2, stats::quantile, probs = 0.95, names = FALSE)

  # Dosis met hoogste % binnen [target_lo, target_hi] onder voorwaarde P(AUCu > safety_max) <= 10%
  safe_idx <- which(pct_above_safety <= 10.0)
  opt_idx  <- if (length(safe_idx) > 0) safe_idx[which.max(pct_in_target[safe_idx])] else which.min(pct_above_safety)

  cur_idx <- which.min(abs(grid - cur_dose))

  list(
    grid_df = data.frame(
      dose             = grid,
      pct_below        = pct_below,
      pct_in_target    = pct_in_target,
      pct_above_target = pct_above_target,
      pct_above_safety = pct_above_safety,
      med_aucu         = med_vec,
      p05_aucu         = p05_vec,
      p95_aucu         = p95_vec
    ),
    cur_dose         = cur_dose,
    cur_pta          = pct_in_target[cur_idx],
    cur_below        = pct_below[cur_idx],
    cur_above_safety = pct_above_safety[cur_idx],
    dose_start_rec   = d_start_rec,
    dose_target_rec  = d_target_rec,
    dose_opt_pta     = grid[opt_idx],
    max_pta_pct      = pct_in_target[opt_idx],
    dose_safety_cap  = d_safety_cap
  )
}

# Vergelijk direct op exact dezelfde virtuele populatie (en dezelfde random BSV-trekkingen)
# wat de drie UGT1A9-ontogeniefuncties voorspellen voor blootstelling én benodigde dosis.
compare_ontogeny_dosing <- function(pop,
                                    st,
                                    target_lo = 20,
                                    target_hi = 75,
                                    safety_max = 90) {
  th <- nonmem_theta()
  eff_ka_mode <- if (!isTRUE(st$ka_bsv)) "none" else st$ka_mode
  om <- nonmem_omega(eff_ka_mode)
  n  <- nrow(pop)

  # Trek exact dezelfde individuele eGFR en CL-BSV als in simulate_population_pk()
  set.seed(st$seed)
  egfr_info <- egfr_from_copula(pop, egfr_iiv = TRUE)
  eta_cl    <- if (om$cl > 0) exp(stats::rnorm(n, 0, sqrt(om$cl))) else rep(1.0, n)

  alb_ov <- if (!is.na(st$alb_override)) st$alb_override else NULL
  cur_dose <- if (st$dose_mode == "mg") st$dose_mg else st$dose_mgkg

  mat_keys <- c("farhan2024", "badee2019", "simcyp_s1")
  mat_lbls <- c(
    farhan2024 = "Farhan et al. 2024 (in-vivo UGT1A9, standaard)",
    badee2019  = "Badee et al. 2019 (Jo 2021 PBPK s2)",
    simcyp_s1  = "Simcyp standaard (Jo 2021 PBPK s1, conservatief)"
  )
  mat_short <- c(
    farhan2024 = "Farhan 2024\n(standaard)",
    badee2019  = "Badee 2019\n(PBPK s2)",
    simcyp_s1  = "Simcyp s1\n(conservatief)"
  )

  rows <- list()
  dist_list <- list()
  quant_rows <- list()
  start_mgkg_vec  <- numeric(length(mat_keys))
  target_mgkg_vec <- numeric(length(mat_keys))

  for (k in seq_along(mat_keys)) {
    mk <- mat_keys[k]
    pk <- individual_params(
      pop,
      th            = th,
      egfr          = egfr_info$egfr,
      mat_option    = mk,
      rf_factor     = st$rf_factor,
      cl_hep_factor = st$cl_hep_factor,
      edema_pct     = st$edema_pct,
      alb_override  = alb_ov
    )
    cl_i <- pk$CL * eta_cl
    wt_dose <- if (identical(st$dose_weight_basis, "actual")) pk$BWT else pk$BWT_DRY
    med_wt  <- stats::median(wt_dose, na.rm = TRUE)

    # Steady-state AUCu,24h per 1 mg/kg en per 1 mg bij het ingestelde interval tau en f1_rel
    aucu24_per_mg   <- pk$FU * (1000 * st$f1_rel / cl_i) * (24 / st$tau)
    aucu24_per_mgkg <- aucu24_per_mg * wt_dose

    aucu_cur <- if (st$dose_mode == "mg") cur_dose * aucu24_per_mg else cur_dose * aucu24_per_mgkg

    med_u_mgkg <- stats::median(aucu24_per_mgkg, na.rm = TRUE)
    p95_u_mgkg <- stats::quantile(aucu24_per_mgkg, 0.95, names = FALSE, na.rm = TRUE)

    d_start_mgkg  <- min(20.0 / med_u_mgkg, safety_max / p95_u_mgkg, 5.0 / med_wt)
    d_target_mgkg <- min(42.0 / med_u_mgkg, safety_max / p95_u_mgkg, 10.0 / med_wt)
    start_mgkg_vec[k]  <- d_start_mgkg
    target_mgkg_vec[k] <- d_target_mgkg

    d_start_mg  <- min(5.0,  d_start_mgkg * med_wt)
    d_target_mg <- min(10.0, d_target_mgkg * med_wt)

    q_med <- stats::median(aucu_cur, na.rm = TRUE)
    q_05  <- stats::quantile(aucu_cur, 0.05, names = FALSE, na.rm = TRUE)
    q_95  <- stats::quantile(aucu_cur, 0.95, names = FALSE, na.rm = TRUE)
    p_in  <- 100 * mean(aucu_cur >= target_lo & aucu_cur <= target_hi)
    p_hi  <- 100 * mean(aucu_cur > safety_max)

    rows[[k]] <- data.frame(
      `UGT1A9-ontogeniemodel`         = mat_lbls[[mk]],
      `Mediane MAT`                   = sprintf("%.2f", stats::median(pk$MAT, na.rm = TRUE)),
      `AUCu,24h,ss huidige dosis [P5-P95]` = sprintf("%.1f [%.1f - %.1f]", q_med, q_05, q_95),
      `PTA streefvenster`             = sprintf("%.0f%%", p_in),
      `> Veiligheidsplafond`          = sprintf("%.0f%%", p_hi),
      `Startdosis (~5 mg-eq)`         = sprintf("%.2f mg/kg (%.2f mg)", d_start_mgkg, d_start_mg),
      `Streefdosis (~10 mg-eq)`       = sprintf("%.2f mg/kg (%.2f mg)", d_target_mgkg, d_target_mg),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    quant_rows[[k]] <- data.frame(
      key = mk, med = q_med, p05 = q_05, p95 = q_95, pct_in = p_in, pct_above_safety = p_hi,
      stringsAsFactors = FALSE
    )

    dist_list[[k]] <- data.frame(
      ontogenie = factor(mat_short[[mk]], levels = unname(mat_short)),
      AUCu_24h  = aucu_cur
    )
  }

  # Beknopte klinisch-farmacologische aandachtspunten op maat van het gesimuleerde cohort
  min_age    <- min(pop$AGE, na.rm = TRUE)
  max_age    <- max(pop$AGE, na.rm = TRUE)
  med_age    <- stats::median(pop$AGE, na.rm = TRUE)
  med_wt_dry <- stats::median(pop$WEIGHT / (1 + st$edema_pct / 100), na.rm = TRUE)
  ont_ratio  <- max(target_mgkg_vec) / max(min(target_mgkg_vec), 1e-6)

  alerts <- list()
  if (min_age < 0.8 && (max_age - min_age) > 1.2) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = "Heterogeen leeftijdsstratum over het eerste levensjaar:",
      body = sprintf("Tussen 0 en 12 maanden nemen zowel UGT1A9-activiteit als GFR niet-lineair toe. Stratificeer de huidige selectie (%.2f–%.2f jaar) bij voorkeur in smalle subcohorten (bv. 0.08–0.25, 0.25–0.5 en 0.5–1 jaar).", min_age, max_age)
    )
  }
  if (ont_ratio >= 1.5) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = sprintf("Divergentie in UGT1A9-ontogenie (factor %.1fx):", ont_ratio),
      body = sprintf("Onder 6–12 maanden wijken snelle in-vivo rijping (Farhan 2024) en conservatieve PBPK-ontogenie (Simcyp s1) substantieel uiteen. Initieer bij voorkeur op de conservatieve startdosis (%.2f mg/kg/dag) met stapsgewijze optitratie op geleide van klinische tolerantie.", min(start_mgkg_vec))
    )
  }
  if (min_age < 2.0) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = "Nauwere metabole veiligheidsmarge (< 2 jaar, 'sick-day rules'):",
      body = "Door beperkte hepatische glycogeenreserves kan SGLT2-geïnduceerde normoglykemische glucosurie bij intercurrent vasten, braken of gastro-enteritis leiden tot hypoglykemie en euglykemische ketoacidose; onderbreek toediening tijdelijk bij verminderde enterale intake."
    )
  }
  if (st$rf_factor < 0.85 || med_age < 0.25) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = "PK–PD dissociatie bij verlaagde GFR (geen opwaartse diuretische titratie):",
      body = "Bij onrijpe neonatale nefronen of cardiorenaal syndroom neemt AUCu,24h toe terwijl de natriuretische respons begrensd wordt door de verlaagde gefiltreerde solute-vracht (luminale SGLT2-bezetting ligt al > 80–90%). Bewaak op prerenale AKI bij gelijktijdig gebruik van lisdiuretica en RAS-remmers."
    )
  }
  if (!is.na(st$alb_override) && st$alb_override < 34) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = sprintf("Hypoalbuminemie ([Alb] = %g g/L) — interpretatie van TDM-spiegels:", st$alb_override),
      body = "Een verhoogde vrije fractie fu verhoogt voor dit laag-extractiefarmacon de totale CL/F en verlaagt totale plasmaconcentraties C(t), terwijl de farmacologisch actieve ongebonden blootstelling (AUCu,24h) ongewijzigd blijft."
    )
  }
  if (st$edema_pct > 0) {
    alerts[[length(alerts) + 1L]] <- list(
      lead = sprintf("Extracellulaire overvulling (+%g%%) — doseer op euvolemisch drooggewicht:", st$edema_pct),
      body = sprintf("Interstitieel oedeem vergroot het verdelingsvolume maar draagt niet bij aan hepatische of renale UGT1A9-capaciteit; baseer mg/kg-dosering op het geschatte drooggewicht (mediaan %.1f kg vs. %.1f kg actueel).", med_wt_dry, stats::median(pop$WEIGHT, na.rm = TRUE))
    )
  }

  list(
    summary_df      = do.call(rbind, rows),
    quant_df        = do.call(rbind, quant_rows),
    dist_df         = do.call(rbind, dist_list),
    start_rec_mgkg  = min(start_mgkg_vec),
    target_min_mgkg = min(target_mgkg_vec),
    target_max_mgkg = max(target_mgkg_vec),
    med_wt_dry      = med_wt_dry,
    ontogeny_ratio  = ont_ratio,
    alerts          = alerts
  )
}

# Gestratificeerde pediatrische hartfalen-referentietabel (0 mnd - 18 jaar)
# bij de actuele hartfalen-pathofysiologie en formulering-instellingen.
generate_pediatric_hf_regimen_table <- function(rf_factor = 1.0,
                                                cl_hep_factor = 1.0,
                                                edema_pct = 0.0,
                                                alb_override = NA_real_,
                                                f1_rel = 1.0,
                                                tau = 24,
                                                safety_max = 90) {
  strata <- data.frame(
    groep = c(
      "Pasgeborene (0 - 1 mnd)",
      "Jonge zuigeling (1 - 3 mnd)",
      "Zuigeling (3 - 6 mnd)",
      "Oudere zuigeling (6 - 12 mnd)",
      "Dreumes (1 - 2 jaar)",
      "Peuter / kleuter (2 - 6 jaar)",
      "Schoolkind (6 - 12 jaar)",
      "Adolescent (12 - 18 jaar)"
    ),
    age_mid = c(0.04, 0.17, 0.375, 0.75, 1.50, 4.00, 9.00, 15.00),
    wt_dry  = c(3.5,  5.0,  6.7,   8.8,  11.0, 16.5, 30.0, 56.0),
    ht_m    = c(0.51, 0.57, 0.63,  0.71, 0.81, 1.03, 1.34, 1.66),
    focus   = c(
      "Zeer lage GFR + onrijp UGT1A9; hoog hypoglykemie/AKI-risico -> uiterste terughoudendheid",
      "4x verschil Farhan vs Simcyp s1; start conservatief + strikte glucose/vochtmonitoring",
      "Snelle GFR- & UGT1A9-rijping; start laag, titreer op geleide van tolerantie/groei",
      "GFR nadert kinderwaarde; let op sick-day rules (vasten/gastro-enteritis) & lisdiuretica",
      "Hogere klaring per kg (allometrie): 0.1 -> 0.2 mg/kg geeft 5 -> 10 mg-eq AUCu",
      "Stabiele allometrische fase: start 0.1 mg/kg, streef 0.18-0.20 mg/kg (op drooggewicht)",
      "Overgang naar vaste tabletsterktes mogelijk (start 2.5 mg, streef 5 mg rond ~25-30 kg)",
      "Volwassen klaring bereikt: start 5 mg, streef 10 mg 1x/dag (plafond 10 mg)"
    ),
    stringsAsFactors = FALSE
  )

  th <- nonmem_theta()
  om <- nonmem_omega("calibrated")
  alb_ov <- if (!is.na(alb_override)) alb_override else NULL
  n_mc <- 400L

  out_rows <- lapply(seq_len(nrow(strata)), function(i) {
    s <- strata[i, ]
    set.seed(1000L + i)
    # Synthetisch gestandaardiseerd cohort rond het stratum-midden (50% jongen / 50% meisje + fysiologische spreiding)
    sex_v <- rep(c(1, 0), length.out = n_mc)
    wt_obs <- s$wt_dry * (1 + edema_pct / 100) * exp(stats::rnorm(n_mc, 0, 0.10))
    ht_v   <- s$ht_m * exp(stats::rnorm(n_mc, 0, 0.03))
    pop_s  <- data.frame(AGE = rep(s$age_mid, n_mc), WEIGHT = wt_obs, HEIGHT = ht_v, SEXISMALE = sex_v)

    egfr_s <- ngfr_rhodin(pop_s$AGE, s$wt_dry) * exp(stats::rnorm(n_mc, 0, 0.20))
    egfr_s <- pmin(pmax(egfr_s, 5), 180)
    eta_cl <- exp(stats::rnorm(n_mc, 0, sqrt(om$cl)))

    p_far <- individual_params(pop_s, th = th, egfr = egfr_s, mat_option = "farhan2024",
                               rf_factor = rf_factor, cl_hep_factor = cl_hep_factor,
                               edema_pct = edema_pct, alb_override = alb_ov)
    p_sim <- individual_params(pop_s, th = th, egfr = egfr_s, mat_option = "simcyp_s1",
                               rf_factor = rf_factor, cl_hep_factor = cl_hep_factor,
                               edema_pct = edema_pct, alb_override = alb_ov)

    aucu_mgkg_far <- p_far$FU * (1000 * f1_rel * p_far$BWT_DRY / (p_far$CL * eta_cl)) * (24 / tau)
    aucu_mgkg_sim <- p_sim$FU * (1000 * f1_rel * p_sim$BWT_DRY / (p_sim$CL * eta_cl)) * (24 / tau)

    d_start_sim <- min(20.0 / stats::median(aucu_mgkg_sim), safety_max / stats::quantile(aucu_mgkg_sim, 0.95, names = FALSE), 5.0 / s$wt_dry)
    d_start_far <- min(20.0 / stats::median(aucu_mgkg_far), safety_max / stats::quantile(aucu_mgkg_far, 0.95, names = FALSE), 5.0 / s$wt_dry)
    d_start_con <- min(d_start_sim, d_start_far)

    d_targ_sim  <- min(42.0 / stats::median(aucu_mgkg_sim), safety_max / stats::quantile(aucu_mgkg_sim, 0.95, names = FALSE), 10.0 / s$wt_dry)
    d_targ_far  <- min(42.0 / stats::median(aucu_mgkg_far), safety_max / stats::quantile(aucu_mgkg_far, 0.95, names = FALSE), 10.0 / s$wt_dry)
    d_lo        <- min(d_targ_sim, d_targ_far)
    d_hi        <- max(d_targ_sim, d_targ_far)

    # Typische PD bij de Farhan streefdosis
    typ_pair <- data.frame(AGE = c(s$age_mid, s$age_mid),
                           WEIGHT = c(s$wt_dry, s$wt_dry) * (1 + edema_pct / 100),
                           HEIGHT = c(s$ht_m, s$ht_m), SEXISMALE = c(1, 0))
    r_pd <- simulate_population_pk(
      typ_pair, dose_mgkg = d_targ_far, tau = tau, steady_state = TRUE,
      iiv = FALSE, egfr_iiv = FALSE, mat_option = "farhan2024",
      rf_factor = rf_factor, cl_hep_factor = cl_hep_factor,
      edema_pct = edema_pct, dose_weight_basis = "dry",
      alb_override = alb_ov, f1_rel = f1_rel
    )

    targ_txt <- if (abs(d_hi - d_lo) / max(d_lo, 1e-4) < 0.15) {
      sprintf("%.2f mg/kg (%.2f mg)", d_hi, min(10, d_hi * s$wt_dry))
    } else {
      sprintf("%.2f - %.2f mg/kg (%.2f - %.2f mg)", d_lo, d_hi, min(10, d_lo * s$wt_dry), min(10, d_hi * s$wt_dry))
    }

    data.frame(
      `Leeftijdsgroep`                 = s$groep,
      `Drooggewicht`                   = sprintf("%.1f kg", s$wt_dry),
      `Eff. eGFR`                      = sprintf("%.0f", stats::median(p_far$EGFR)),
      `MAT (Simcyp - Farhan)`          = sprintf("%.2f - %.2f", p_sim$MAT[1], p_far$MAT[1]),
      `Veilige startdosis (5 mg-eq)`   = sprintf("%.2f mg/kg (%.2f mg)", d_start_con, min(5, d_start_con * s$wt_dry)),
      `Streefdosis-bereik (10 mg-eq)`  = targ_txt,
      `SGLT2-remming / Renale PD`      = sprintf("%.0f%% / %.0f%%", mean(r_pd$metrics$SGLT2_inh_avg), mean(r_pd$metrics$Renal_PD_eff)),
      `Klinisch PK/PD-aandachtspunt`   = s$focus,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out_rows)
}

# Visuele liniaal van het therapeutisch venster: laat in één oogopslag zien waar de
# grenzen (20, 75 en 90 ng*h/mL) vandaan komen uit volwassen studies én waar het
# gesimuleerde kindercohort op deze schaal valt.
build_window_explainer_plot <- function(s, tg, ont) {
  st  <- s$res$settings
  qdf <- ont$quant_df

  chosen_key <- st$mat_option
  alt_key    <- if (identical(chosen_key, "simcyp_s1")) "farhan2024" else "simcyp_s1"
  short_map  <- c(
    farhan2024 = "Farhan 2024 (in-vivo UGT1A9)",
    badee2019  = "Badee 2019 (PBPK s2)",
    simcyp_s1  = "Simcyp s1 (conservatief)"
  )

  q_main <- qdf[qdf$key == chosen_key, , drop = FALSE]
  q_alt  <- qdf[qdf$key == alt_key, , drop = FALSE]

  bars <- data.frame(
    y = c(6, 5, 4, 3, 2, 1),
    groep = c(
      sprintf("Gesimuleerd cohort — %s", short_map[[chosen_key]]),
      sprintf("Gesimuleerd cohort — %s", short_map[[alt_key]]),
      "Volwassen 5 mg 1x/dag (referentie startdosis, eGFR 90)",
      "Volwassen 10 mg 1x/dag (normale nierfunctie, eGFR >= 90)",
      "Volwassen 10 mg 1x/dag (DAPA-HF / milde CKD, eGFR 60-89)",
      "Volwassen 10 mg 1x/dag (ernstige CKD, eGFR 15-29)"
    ),
    med = c(q_main$med, q_alt$med, 19.2, 38.5, 48.2, 87.7),
    p05 = c(q_main$p05, q_alt$p05, 11.8, 23.7, 29.2, 56.8),
    p95 = c(q_main$p95, q_alt$p95, 29.0, 57.9, 74.5, 119.5),
    type = c("ped_main", "ped_alt", "ad_5", "ad_10", "ad_hf", "ad_ckd"),
    stringsAsFactors = FALSE
  )

  x_max <- max(132, min(220, max(c(q_main$p95, tg$safety_max * 1.35), na.rm = TRUE) * 1.08))
  bars$p95_plot <- pmin(bars$p95, x_max * 0.985)
  bars$med_plot <- pmin(bars$med, x_max * 0.97)
  bars$lbl <- ifelse(
    bars$p95 > x_max * 0.985,
    sprintf("%.1f [%.0f - %.0f+]", bars$med, bars$p05, bars$p95),
    sprintf("%.1f [%.0f - %.0f]", bars$med, bars$p05, bars$p95)
  )

  zone_top <- 7.05
  z_mid_lo <- tg$target_lo / 2
  z_mid_ok <- (tg$target_lo + tg$target_hi) / 2
  z_mid_hi <- (tg$target_hi + tg$safety_max) / 2
  z_mid_mx <- (tg$safety_max + x_max) / 2

  ggplot(bars, aes(y = y)) +
    annotate("rect", xmin = 0, xmax = tg$target_lo, ymin = 0.4, ymax = 6.60, fill = "#eef3f8", alpha = 0.85) +
    annotate("rect", xmin = tg$target_lo, xmax = tg$target_hi, ymin = 0.4, ymax = 6.60, fill = "#dff2e6", alpha = 0.85) +
    annotate("rect", xmin = tg$target_hi, xmax = tg$safety_max, ymin = 0.4, ymax = 6.60, fill = "#fdf0d5", alpha = 0.85) +
    annotate("rect", xmin = tg$safety_max, xmax = x_max, ymin = 0.4, ymax = 6.60, fill = "#fbe2e2", alpha = 0.85) +
    geom_vline(xintercept = tg$target_lo, color = "#2b5c8f", linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = tg$target_hi, color = "#1a7f37", linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = tg$safety_max, color = "#b2182b", linetype = "dashed", linewidth = 0.85) +
    geom_hline(yintercept = 4.5, color = "#7b8a9b", linetype = "dotted", linewidth = 0.7) +
    annotate("text", x = z_mid_lo, y = zone_top, label = sprintf("Subtherap.\n(< %g)", tg$target_lo),
             size = 3.0, fontface = "bold", color = "#345474", lineheight = 0.95) +
    annotate("text", x = z_mid_ok, y = zone_top, label = sprintf("Streefvenster DAPA-HF\n(%g – %g ng*h/mL)", tg$target_lo, tg$target_hi),
             size = 3.3, fontface = "bold", color = "#145a2b", lineheight = 0.95) +
    annotate("text", x = z_mid_hi, y = zone_top, label = sprintf("Hoog\n(%g–%g)", tg$target_hi, tg$safety_max),
             size = 3.0, fontface = "bold", color = "#8a5300", lineheight = 0.95) +
    annotate("text", x = z_mid_mx, y = zone_top, label = sprintf("> Veiligheidsplafond\n(> %g)", tg$safety_max),
             size = 3.05, fontface = "bold", color = "#991b1b", lineheight = 0.95) +
    geom_errorbarh(aes(xmin = p05, xmax = p95_plot, color = type, linetype = type),
                   height = 0.22, linewidth = 1.15) +
    geom_point(aes(x = med_plot, color = type, shape = type), size = 3.8) +
    geom_text(aes(x = med_plot, y = y + 0.28, label = lbl, color = type),
              size = 3.2, fontface = "bold") +
    scale_y_continuous(breaks = bars$y, labels = bars$groep) +
    scale_x_continuous(expand = c(0.01, 0)) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0.4, 7.55), clip = "off") +
    scale_color_manual(values = c(
      ped_main = "#08306b",
      ped_alt  = "#35689a",
      ad_5     = "#2b6cb0",
      ad_10    = "#1a7f37",
      ad_hf    = "#145a2b",
      ad_ckd   = "#b2182b"
    )) +
    scale_linetype_manual(values = c(
      ped_main = "solid",
      ped_alt  = "dashed",
      ad_5     = "solid",
      ad_10    = "solid",
      ad_hf    = "solid",
      ad_ckd   = "solid"
    )) +
    scale_shape_manual(values = c(
      ped_main = 18,
      ped_alt  = 16,
      ad_5     = 15,
      ad_10    = 16,
      ad_hf    = 17,
      ad_ckd   = 15
    )) +
    labs(
      title = "Afleiding therapeutisch streefvenster: volwassen referenties vs. gesimuleerd cohort",
      subtitle = "Mediaan en 90%-predictie-interval [P5 – P95] van steady-state AUCu,24h (Penland et al. 2023 Tabel S4; DAPA-HF)",
      x = "Ongebonden blootstelling AUCu,24h,ss (ng*h/mL)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", color = "#08306b", size = 12.5),
      plot.subtitle = element_text(color = "#3d4f63", size = 10.5),
      plot.margin = margin(8, 20, 8, 8),
      axis.text.y = element_text(face = "bold", color = "#1f2d3d", size = 10.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

