# =====================================================================
# Populatie-PK model voor dapagliflozine — R-implementatie van het
# NONMEM model van Penland et al. (J Clin Pharmacol 2023;63:551-559,
# doi:10.1002/jcph.2196) met de gevraagde aanpassingen:
#   - CL gesplitst in niet-renale klaring (gematureerd, UGT1A9) en
#     renale klaring (THETA(6)*EGFR)
#   - allometrische schaling van alle parameters t.o.v. (BWT/84)
#   - maturatie van niet-renale klaring volgens Farhan et al. 2024
#     (doi:10.1002/jcph.2484): MAT = Fbirth + (1-Fbirth)*AGE^Hill/(TM50^Hill+AGE^Hill)
#
# Structuur (ADVAN4 TRANS4): 3 compartimenten — depot (1) met lagtijd
# ALAG1, centraal (2, volume V2, eliminatie CL), perifeer (3, volume V3,
# uitwisseling Q). Eerste-orde absorptie KA van depot naar centraal.
#
# Eenheden: tijd in uur, volumes in L, BWT in kg, AGE in jaar,
# EGFR in mL/min/1.73m2, concentraties in ng/mL (= microg/L).
# =====================================================================

nonmem_theta <- function() {
  list(
    ka     = 2.833,   # THETA(1)  KA (h-1) bij BWT=84 kg
    clnr   = 5.44,    # THETA(2)  niet-renale CL/F, volwassene (L/h) bij BWT=84, MAT=1
    v2     = 73.47,   # THETA(3)  Vc/F (L) bij BWT=84
    v3     = 168.5,   # THETA(4)  Vp/F (L) bij BWT=84
    q      = 9.137,   # THETA(5)  Q/F (L/h) bij BWT=84
    egfr_s = 0.177,   # THETA(6)  renale slope (L/h per mL/min/1.73m2)
    sex    = -0.1441, # THETA(7)  CL/F factor = (1 + sex*(1-SEX)); SEX=1 man, 0 vrouw
    alag   = 0.5,     # THETA(9)  ALAG1 (h), fixed
    tm50   = 0.09,    # THETA(10) UGT1A9 maturatie TM50 (jaar), Farhan 2024 Table 1
    hill   = 1.9,     # THETA(11) UGT1A9 maturatie Hill-coeff, Farhan 2024
    fbirth = 0.05     # THETA(12) UGT1A9 Fbirth, Farhan 2024
  )
}

nonmem_omega <- function() {
  list(
    ka = 8.52,   # BSV op KA, FIX in publicatie (ETA(1))
    cl = 0.1136, # BSV op CL/F (ETA(2))
    v3 = 0.8229  # BSV op Vp/F (ETA(3)); geen BSV op V2 en Q in dit model
  )
}

# EGFR in mL/min/1.73m2 ophalen uit de copulapop-populatie.
# copulapop berekent EGFR zelf, maar alleen als de dataset CREAT bevat
# (in de app: NHANES, RUMC):
#   <18 jaar : bedside Schwartz  eGFR = 0.413 * lengte(cm) / Scr(mg/dL)
#   >=18 jaar: CKD-EPI 2021 (race-free)
# Beide formules leveren mL/min/1.73m2 — precies de eenheid waarin het
# NONMEM covariaat EGFR is gedefinieerd (THETA(6)*EGFR). Er is dus
# GEEN omrekening met BSA nodig (dat zou absoluut GFR in mL/min geven).
# Fallback bij datasets zonder creatinine (in de app: HSE): geschatte
# normale eGFR per 1.73m2 als functie van leeftijd (GFR-maturatie, Hill).
egfr_from_copula <- function(pop) {
  n <- nrow(pop)
  egfr <- rep(NA_real_, n)
  source <- "fallback: geschatte normale eGFR per 1.73m2 per leeftijd (dataset zonder creatinine)"
  if ("EGFR" %in% names(pop) && !all(is.na(pop$EGFR))) {
    egfr <- pop$EGFR
    source <- "copulapop: bedside Schwartz (<18 j) / CKD-EPI 2021 (>=18 j), in mL/min/1.73m2"
  }
  miss <- is.na(egfr)
  if (any(miss)) {
    age <- pop$AGE[miss]
    egfr[miss] <- pmax(30, 120 * age^1.7 / (0.75^1.7 + age^1.7))
  }
  list(egfr = egfr, source = source)
}

# Individuele parameters uit de virtuele populatie (covariaten + allometrie)
individual_params <- function(pop, th = nonmem_theta(), egfr = NULL) {
  n <- nrow(pop)
  bwt <- pop$WEIGHT          # BWT in kg
  age <- pop$AGE             # AGE in jaar
  sex <- if ("SEXISMALE" %in% names(pop)) pop$SEXISMALE else rep(1, n)  # 1=man, 0=vrouw
  if (is.null(egfr)) egfr <- egfr_from_copula(pop)$egfr  # mL/min/1.73m2

  w   <- bwt / 84
  mat <- th$fbirth + (1 - th$fbirth) * age^th$hill / (th$tm50^th$hill + age^th$hill)

  clnr <- th$clnr * mat * w^0.75
  cl   <- (clnr + th$egfr_s * egfr) * (1 + th$sex * (1 - sex))
  ka   <- th$ka * w^(-0.25)
  v2   <- th$v2 * w
  v3   <- th$v3 * w
  q    <- th$q * w^0.75

  list(CL = cl, KA = ka, V2 = v2, V3 = v3, Q = q, ALAG = th$alag,
       BWT = bwt, AGE = age, EGFR = egfr, SEX = sex, MAT = mat)
}

# Monte Carlo simulatie: n_doses doses met interval tau, per kind.
# IIV: ETA(2) op CL, ETA(3) op V3, optioneel ETA(1) op KA (omega^2 = 8.52 FIX).
# IOV: per dosing occasion een eigen log-normale multiplier op F1 (toegediende dosis).
# ODEs worden met RK4 vector-gewijs over alle kinderen gentegreerd.
simulate_population_pk <- function(pop,
                                   dose_mgkg,
                                   tau = 24,
                                   n_doses = 7,
                                   iov_f1_cv = 0,
                                   ka_bsv = TRUE,
                                   seed = 1,
                                   th = nonmem_theta(),
                                   om = nonmem_omega(),
                                   dt = 0.05,
                                   out_by = 0.1) {
  set.seed(seed)
  n <- nrow(pop)
  egfr_info <- egfr_from_copula(pop)
  p <- individual_params(pop, th, egfr = egfr_info$egfr)

  # inter-individuele variabiliteit (log-normaal, exp(ETA))
  p$CL <- p$CL * exp(rnorm(n, 0, sqrt(om$cl)))
  p$V3 <- p$V3 * exp(rnorm(n, 0, sqrt(om$v3)))
  if (ka_bsv) {
    p$KA <- p$KA * exp(rnorm(n, 0, sqrt(om$ka)))
    # numerieke stabilisatie RK4: ka clampen; treft alleen de extreme staart
    # van ETA(1), waar absorptie nagenoeg momentaan is
    p$KA <- pmin(pmax(p$KA, 0.05), 25)
  }

  dose_ug <- dose_mgkg * p$BWT * 1000  # mg -> microg

  sigma_iov <- sqrt(log(1 + (iov_f1_cv / 100)^2))

  times <- seq(0, tau, by = out_by)
  n_out <- length(times)
  steps_per_out <- max(1, round(out_by / dt))
  n_steps <- round(tau / dt)

  A1 <- numeric(n)  # depot (gut)
  A2 <- numeric(n)  # centraal compartiment
  A3 <- numeric(n)  # perifeer compartiment
  curve_first <- NULL
  curve_last <- NULL

  for (d in seq_len(n_doses)) {
    f1 <- if (sigma_iov > 0) exp(rnorm(n, 0, sigma_iov)) else rep(1, n)
    A1 <- A1 + dose_ug * f1
    ke  <- p$CL / p$V2
    k12 <- p$Q / p$V2
    k21 <- p$Q / p$V3
    ka  <- p$KA
    alag <- p$ALAG

    deriv <- function(a1, a2, a3, tt) {
      ka_eff <- ifelse(tt >= alag, ka, 0)  # lagtijd ALAG1: tot dan geen absorptie
      list(-ka_eff * a1,
           ka_eff * a1 - (ke + k12) * a2 + k21 * a3,
           k12 * a2 - k21 * a3)
    }

    curve <- matrix(NA_real_, nrow = n, ncol = n_out)
    curve[, 1] <- A2 / p$V2
    oi <- 2
    t <- 0
    for (s in seq_len(n_steps)) {
      k1 <- deriv(A1, A2, A3, t)
      k2 <- deriv(A1 + 0.5 * dt * k1[[1]], A2 + 0.5 * dt * k1[[2]], A3 + 0.5 * dt * k1[[3]], t + 0.5 * dt)
      k3 <- deriv(A1 + 0.5 * dt * k2[[1]], A2 + 0.5 * dt * k2[[2]], A3 + 0.5 * dt * k2[[3]], t + 0.5 * dt)
      k4 <- deriv(A1 + dt * k3[[1]], A2 + dt * k3[[2]], A3 + dt * k3[[3]], t + dt)
      A1 <- A1 + (dt / 6) * (k1[[1]] + 2 * k2[[1]] + 2 * k3[[1]] + k4[[1]])
      A2 <- A2 + (dt / 6) * (k1[[2]] + 2 * k2[[2]] + 2 * k3[[2]] + k4[[2]])
      A3 <- A3 + (dt / 6) * (k1[[3]] + 2 * k2[[3]] + 2 * k3[[3]] + k4[[3]])
      t <- t + dt
      if (s %% steps_per_out == 0 && oi <= n_out) {
        curve[, oi] <- A2 / p$V2
        oi <- oi + 1
      }
    }
    if (d == 1) curve_first <- curve
    if (d == n_doses) curve_last <- curve
  }

  metrics <- data.frame(
    ID = if ("ID" %in% names(pop)) pop$ID else seq_len(n),
    AGE = p$AGE, BWT = p$BWT, EGFR = p$EGFR, SEXISMALE = p$SEX, MAT = p$MAT,
    CL = p$CL, KA = p$KA, V2 = p$V2, V3 = p$V3, Q = p$Q,
    DOSE_MG = dose_mgkg * p$BWT
  )
  metrics$Cmax <- apply(curve_last, 1, max)
  metrics$Tmax <- times[apply(curve_last, 1, which.max)]
  metrics$Cmin <- curve_last[, n_out]  # trough = C(tau)
  metrics$AUCtau <- apply(curve_last, 1, function(x) {
    sum(diff(times) * (head(x, -1) + tail(x, -1)) / 2)
  })
  metrics$Cavg <- metrics$AUCtau / tau

  list(
    pop = pop,
    params = p,
    times = times,
    curve_ss = curve_last,
    curve_day1 = curve_first,
    metrics = metrics,
    settings = list(dose_mgkg = dose_mgkg, tau = tau, n_doses = n_doses,
                    iov_f1_cv = iov_f1_cv, ka_bsv = ka_bsv, seed = seed, n = n,
                    egfr_source = egfr_info$source)
  )
}

# Populatie-quantielen per tijdspunt (voor mediaan + 90% predictie-interval)
curve_summary <- function(curve, times, probs = c(0.05, 0.5, 0.95)) {
  q <- apply(curve, 2, quantile, probs = probs, names = FALSE)
  data.frame(t = times, q05 = q[1, ], median = q[2, ], q95 = q[3, ])
}

# Samenvatting steady-state metrieken: mediaan + P5/P95
summarize_metrics <- function(metrics, tau = 24) {
  cols <- c("Cmax", "Tmax", "Cmin", "AUCtau", "Cavg")
  labels <- c("Cmax,ss (ng/mL)", "Tmax (uur)", sprintf("Cmin,ss = C%.0f (ng/mL)", tau),
              "AUCtau,ss (ng*h/mL)", "Cavg,ss (ng/mL)")
  data.frame(
    Metriek = labels,
    Mediaan = sapply(metrics[cols], median),
    P5  = sapply(metrics[cols], function(x) quantile(x, 0.05, names = FALSE)),
    P95 = sapply(metrics[cols], function(x) quantile(x, 0.95, names = FALSE))
  )
}
