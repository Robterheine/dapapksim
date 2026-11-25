$SIZES PD=50
$PROBLEM Dapagliflozine 6-pijler fysiologisch popPK/PD model + pediatrisch hartfalen (dapapksim)
; =============================================================================
; Volledige NONMEM-implementatie van het dapapksim-model (R/pk_model.R & R/hf_pkpd.R)
; Basis: Penland et al. J Clin Pharmacol 2023;63:551-559 (doi:10.1002/jcph.2196),
; uitgebreid met:
;   1. FFM/NFM allometrie op euvolemisch drooggewicht (Al-Sallami 2015; Janmahasatian 2005)
;   2. Plasma-albumine ontogenie & vrije fractie fu (McNamara & Alcorn 2002; fu_volw = 0.086)
;   3. Nierfunctie: Rhodin 2009 PMA-nGFR (<1 j / geen CREAT), Pierce 2021 CKiD U25 (1-<18 j),
;      en CKD-EPI 2021 (>=18 j), begrensd op [5, 180] mL/min/1.73m2
;   4. Drie-paden klaring: hepatisch UGT1A9 + renaal tubulair UGT1A9 + passieve filtratie
;   5. UGT1A9-ontogenie keuze via MATOPT: 1 = Farhan 2024 (standaard), 2 = Badee 2019, 3 = Simcyp s1
;   6. Ontwikkelings-GI-absorptie op KA + hartfalen-pathofysiologie (RFFAC, HEPFAC, EDEMA, ALBIN, F1REL)
;   7. Semi-mechanistische duale PD in $ERROR: Cu(t), luminale SGLT2-bezetting (%) en TGF/HCT Emax (%)
;
; Eenheden:
;   TIME, TAD, II : uur (h)
;   AMT, DOSEMG   : milligram (mg) -> via S2 = V2/1000 geeft F = A(2)/S2 exact ng/mL (= ug/L)
;   AGE           : jaar
;   BWT           : actueel lichaamsgewicht (kg; bij EDEMA > 0 geldt WDRY = BWT / (1 + EDEMA/100))
;   HT            : lengte in meter (m; waarden > 3 worden automatisch als cm -> m omgerekend)
;   SEX           : 1 = man/jongen, 0 = vrouw/meisje
;   CREAT         : serumcreatinine in umol/L (zet op -99 indien niet beschikbaar)
;   EGFRIN        : directe eGFR in mL/min/1.73m2 (zet op -99 om uit Rhodin/CREAT te berekenen)
;   ALBIN         : serumalbumine in g/L (zet op -99 voor fysiologische McNamara leeftijdsformule)
;   EDEMA         : percentage overvulling boven drooggewicht (%; 0 = euvolemisch)
;   RFFAC         : cardiorenaal GFR-behoud (fractie; 1.0 = normaal, bv. 0.80 bij stabiel HF)
;   HEPFAC        : hepatische UGT1A9-functie bij veneuze stuwing (fractie; 1.0 = normaal)
;   F1REL         : relatieve enterale biologische beschikbaarheid (fractie; 1.0 = 100%)
;   MATOPT        : UGT1A9-ontogenie (1 = Farhan 2024, 2 = Badee 2019, 3 = Simcyp s1)
;   TAU           : doseringsinterval in uur (bv. 24) voor berekening van AUCU24
; =============================================================================

$INPUT ID TIME TAD DAY EVID AMT ADDL II CMT DV MDV AGE BWT HT SEX CREAT EGFRIN ALBIN EDEMA RFFAC HEPFAC F1REL MATOPT DOSEMG TAU

$DATA example_child_hf.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

$PK
  ; --- 0. Euvolemisch drooggewicht, lengte, BMI en Mosteller-BSA ---
  EDPCT = EDEMA
  IF (EDPCT.LT.0.0) EDPCT = 0.0
  WDRY  = BWT / (1.0 + EDPCT / 100.0)
  VEDEM = BWT - WDRY
  IF (VEDEM.LT.0.0) VEDEM = 0.0

  IF (HT.GT.3.0) THEN
    HM = HT / 100.0
  ELSE
    HM = HT
  ENDIF
  H2 = HM * HM
  IF (H2.LT.0.09) H2 = 0.09
  BMI = WDRY / H2
  BSA = SQRT((HM * 100.0 * WDRY) / 3600.0)

  ; --- 1. Lichaamssamenstelling: FFM (Al-Sallami 2015 / Janmahasatian 2005) & NFM ---
  IF (SEX.EQ.1) THEN
    FFMAD = 9270.0 * WDRY / (6680.0 + 216.0 * BMI)
  ELSE
    FFMAD = 9270.0 * WDRY / (8780.0 + 244.0 * BMI)
  ENDIF

  IF (AGE.GE.18.0) THEN
    MFFM = 1.0
  ELSE
    IF (SEX.EQ.1) THEN
      IF (AGE.LE.0.0) THEN
        MFFM = 0.88
      ELSE
        MFFM = 0.88 + 0.12 / (1.0 + (AGE / 13.4)**(-12.7))
      ENDIF
    ELSE
      IF (AGE.LE.0.0) THEN
        MFFM = 1.11
      ELSE
        MFFM = 1.11 - 0.11 / (1.0 + (AGE / 7.1)**(-1.1))
      ENDIF
    ENDIF
  ENDIF

  FFM = FFMAD * MFFM
  IF (FFM.GT.0.97*WDRY) FFM = 0.97 * WDRY

  NFMCL = FFM + 0.25 * (WDRY - FFM)
  NFMV2 = FFM + 0.20 * (WDRY - FFM)
  FSZCL = (NFMCL / 63.08)**0.75
  FSZV2 = NFMV2 / 61.68

  ; --- 2. Plasma-albumine ontogenie (McNamara & Alcorn 2002) & eiwitbinding fu ---
  FUAD   = 0.086
  ALBAD  = 44.0
  KABIND = (1.0 - FUAD) / FUAD
  DAYS   = AGE * 365.25
  IF (DAYS.LT.1.0) DAYS = 1.0

  IF (AGE.GE.18.0) THEN
    ALBMOD = ALBAD
  ELSE
    ALBMOD = 1.1287 * LOG(DAYS) + 33.746
    IF (ALBMOD.GT.ALBAD) ALBMOD = ALBAD
  ENDIF

  IF (ALBIN.GT.0.0) THEN
    ALBU = ALBIN
  ELSE
    ALBU = ALBMOD
  ENDIF

  FU    = 1.0 / (1.0 + (ALBU / ALBAD) * KABIND)
  FURAT = FU / FUAD

  ; --- 3. Verwachte gezonde GFR (Rhodin 2009) & individuele eGFR ---
  WGFR = WDRY
  IF (WGFR.LT.0.4) WGFR = 0.4
  GA = 40.0
  IF (AGE.LT.0.1.AND.WGFR.LT.2.5) THEN
    GA = 22.0 + 7.2 * WGFR
    IF (GA.LT.24.0) GA = 24.0
    IF (GA.GT.40.0) GA = 40.0
  ENDIF
  PMA   = GA + 52.1775 * AGE
  AGE35 = AGE - 35.0
  IF (AGE35.LT.0.0) AGE35 = 0.0
  NGFR  = 115.0 * ((PMA**3.4) / (47.7**3.4 + PMA**3.4)) * (0.9925**AGE35)

  IF (EGFRIN.GT.0.0) THEN
    EGFR0 = EGFRIN
  ELSE
    IF (CREAT.GT.0.0.AND.AGE.GE.18.0) THEN
      ; Volwassen CKD-EPI 2021 (race-free)
      SCR = CREAT / 88.4
      IF (SCR.LT.0.05) SCR = 0.05
      IF (SEX.EQ.1) THEN
        KAP  = 0.9
        ALP  = -0.302
        SFAC = 1.0
      ELSE
        KAP  = 0.7
        ALP  = -0.241
        SFAC = 1.012
      ENDIF
      CRK   = SCR / KAP
      MINCR = CRK
      IF (MINCR.GT.1.0) MINCR = 1.0
      MAXCR = CRK
      IF (MAXCR.LT.1.0) MAXCR = 1.0
      EGFR0 = 142.0 * (MINCR**ALP) * (MAXCR**(-1.200)) * (0.9938**AGE) * SFAC
      IF (EGFR0.LT.5.0) EGFR0 = 5.0
      IF (EGFR0.GT.180.0) EGFR0 = 180.0
    ELSE
      IF (CREAT.GT.0.0.AND.AGE.GE.1.0) THEN
        ; Pediatrische Pierce 2021 CKiD U25 (1 - <18 jaar)
        SCR = CREAT / 88.4
        IF (SCR.LT.0.05) SCR = 0.05
        IF (SEX.EQ.1) THEN
          IF (AGE.LT.12.0) THEN
            KCKID = 39.0 * (1.008**(AGE - 12.0))
          ELSE
            KCKID = 39.0 * (1.045**(AGE - 12.0))
          ENDIF
        ELSE
          IF (AGE.LT.12.0) THEN
            KCKID = 36.1 * (1.008**(AGE - 12.0))
          ELSE
            KCKID = 36.1 * (1.023**(AGE - 12.0))
          ENDIF
        ENDIF
        EGFR0 = KCKID * HM / SCR
        IF (EGFR0.LT.5.0) EGFR0 = 5.0
        IF (EGFR0.GT.180.0) EGFR0 = 180.0
      ELSE
        ; Zuigelingen (<1 jaar) of ontbrekend creatinine: Rhodin PMA-maturatie
        EGFR0 = NGFR
        IF (EGFR0.LT.5.0) EGFR0 = 5.0
        IF (EGFR0.GT.180.0) EGFR0 = 180.0
      ENDIF
    ENDIF
  ENDIF

  ; Cardiorenaal GFR-behoud (RFFAC)
  EGFR = EGFR0
  IF (RFFAC.LT.0.999999.OR.RFFAC.GT.1.000001) THEN
    EGFR = EGFR0 * RFFAC
    IF (EGFR.LT.5.0) EGFR = 5.0
    IF (EGFR.GT.180.0) EGFR = 180.0
  ENDIF

  NGFRD = NGFR
  IF (NGFRD.LT.1.0) NGFRD = 1.0
  RF = EGFR / NGFRD
  IF (RF.LT.0.0) RF = 0.0
  IF (RF.GT.1.75) RF = 1.75

  ; --- 4 & 5. UGT1A9-ontogenie (MAT) & drie-paden klaringsdecompositie ---
  TM50   = THETA(11)
  HILL   = THETA(12)
  FBIRTH = THETA(13)

  IF (MATOPT.EQ.2) THEN
    ; Badee et al. 2019 / Jo et al. 2021 scenario 2
    IF (AGE.LE.0.0) THEN
      MAT = 0.0
    ELSE
      IF (AGE.LE.1.0) THEN
        MAT = 1.305 * (AGE**1.1) / (0.086**1.1 + AGE**1.1)
      ELSE
        IF (AGE.LE.20.0) THEN
          MAT = 1.243 - 0.833 * EXP(0.2 * (AGE - 25.0))
          IF (MAT.LT.1.0) MAT = 1.0
        ELSE
          MAT = 1.0
        ENDIF
      ENDIF
    ENDIF
  ELSE
    IF (MATOPT.EQ.3) THEN
      ; Simcyp s1 / Jo et al. 2021 scenario 1
      IF (AGE.LE.15.0) THEN
        MAT = 0.086 + (1.041 - 0.086) * (AGE**1.2) / (1.16**1.2 + AGE**1.2)
      ELSE
        MAT = 1.0
      ENDIF
    ELSE
      ; Standaard (MATOPT = 1): Farhan et al. 2024 in-vivo UGT1A9
      IF (AGE.LE.0.0) THEN
        MAT = FBIRTH
      ELSE
        MAT = FBIRTH + (1.0 - FBIRTH) * (AGE**HILL) / (TM50**HILL + AGE**HILL)
      ENDIF
    ENDIF
  ENDIF

  CLHEP  = FURAT * THETA(2) * MAT * FSZCL * HEPFAC
  CLRUGT = FURAT * THETA(3) * MAT * FSZCL * RF
  CLFILT = FURAT * THETA(4) * EGFR * (BSA / 2.0)
  TVCL   = CLHEP + CLRUGT + CLFILT

  ; --- 6. Volumes (incl. extracellulair oedeem in V2), Q en GI-absorptierijping ---
  TVV2 = THETA(5) * FSZV2 + VEDEM
  TVV3 = THETA(6) * (WDRY / 84.0) * FURAT
  TVQ  = THETA(7) * FSZCL
  FGI  = THETA(9) + (1.0 - THETA(9)) * (AGE / (THETA(10) + AGE))
  TVKA = THETA(1) * FGI

  ; --- Individuele PK-parameters (ADVAN4 TRANS4) ---
  CL    = TVCL * EXP(ETA(2))
  V2    = TVV2
  V3    = TVV3 * EXP(ETA(3))
  Q     = TVQ
  KARAW = TVKA * EXP(ETA(1))
  IF (KARAW.LT.0.02) THEN
    KA = 0.02
  ELSE
    IF (KARAW.GT.50.0) THEN
      KA = 50.0
    ELSE
      KA = KARAW
    ENDIF
  ENDIF

  ALAG1 = THETA(8)
  F1    = F1REL
  S2    = V2 / 1000.0 ; AMT in mg, V2 in L => F = A(2)/S2 in ng/mL (ug/L)

  ; Analytische steady-state 24-uurs vrije blootstelling & systemische Emax-markers
  AUCINF = (DOSEMG * 1000.0 * F1REL) / CL
  AUCU24 = FU * AUCINF * (24.0 / TAU)
  ETGF   = 100.0 * AUCU24 / (5.2 + AUCU24)
  EHCT   = 100.0 * AUCU24 / (20.7 + AUCU24)

$ERROR
  CTOT  = F
  CU    = FU * CTOT
  SGLT2 = 100.0 * CU / (0.20 + CU)
  Y     = CTOT + EPS(1)

$THETA
  2.833 FIX  ; 1  KA_adult (1/h) - Penland 2023 THETA(1)
  5.44  FIX  ; 2  CL_hep_UGT1A9 (L/h) bij ref-NFM, MAT=1, fu=0.086
  15.78 FIX  ; 3  CL_ren_UGT1A9 (L/h) bij ref-NFM, MAT=1, RF=1, fu=0.086
  0.011 FIX  ; 4  CL_filt_slope (L/h per mL/min/1.73m2 bij BSA=2.0 m2)
  73.47 FIX  ; 5  V2_ref (L bij NFM_V2 = 61.68 kg) - Penland 2023 THETA(3)
  168.5 FIX  ; 6  V3_ref (L bij BWT = 84 kg, fu = 0.086) - Penland 2023 THETA(4)
  9.137 FIX  ; 7  Q_ref (L/h bij NFM_CL = 63.08 kg) - Penland 2023 THETA(5)
  0.5   FIX  ; 8  ALAG1 absorptielagtijd (h)
  0.45  FIX  ; 9  GI_F0 relatieve neonatale maagledigingssnelheid bij geboorte
  0.10  FIX  ; 10 GI_T50 rijpingshalfwaardetijd maaglediging (jaar)
  0.09  FIX  ; 11 TM50 UGT1A9 maturatie (jaar) - Farhan et al. 2024
  1.9   FIX  ; 12 HILL UGT1A9 maturatie - Farhan et al. 2024
  0.05  FIX  ; 13 FBIRTH UGT1A9 fractie bij geboorte - Farhan et al. 2024

; Standaard op 0 FIX voor exacte 1-op-1 deterministische reproductie van het voorbeeldkind.
; Voor populatie Monte Carlo simulatie met inter-individuele variabiliteit (IIV):
; vervang onderstaande drie regels door:
;   0.36   FIX  ; ETA(1) BSV op KA (gekalibreerd klinisch; of 8.52 FIX uit Penland 2023)
;   0.1136 FIX  ; ETA(2) BSV op CL/F
;   0.8229 FIX  ; ETA(3) BSV op V3/F
$OMEGA
  0 FIX  ; ETA(1) BSV KA (deterministische modus: ETA = 0)
  0 FIX  ; ETA(2) BSV CL
  0 FIX  ; ETA(3) BSV V3

$SIGMA
  0 FIX  ; Geen residuele fout (zuivere PK/PD-voorspelling)

$SIMULATION (20261125) ONLYSIM SUBPROBLEMS=1

$TABLE ID TIME TAD DAY EVID AMT WDRY FFM ALBU FU NGFR EGFR RF MAT
       CL CLHEP CLRUGT CLFILT V2 V3 Q KA CTOT CU SGLT2 AUCU24 ETGF EHCT
       NOAPPEND NOPRINT ONEHEADER FORMAT=s1PE14.6 FILE=sdtab_dapa_ped_hf.tab
