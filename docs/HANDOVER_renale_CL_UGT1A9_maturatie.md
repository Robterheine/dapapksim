================================================================================
HANDOVER — Bepaling renale/niet-renale klaring en UGT1A9-maturatie
Dapagliflozin popPK model (basis: Penland et al., J Clin Pharmacol 2023;63:551-559)
Bestand: plain text, bedoeld voor hergebruik in een andere sessie / artikel
================================================================================

0. BRONNEN
-----------
- Paper: Penland RC et al. "Evaluation of the PK of Dapagliflozin in Patients With
  CKD With or Without T2DM". J Clin Pharmacol 2023;63(5):551-559 (DOI 10.1002/jcph.2196).
  Fig. 2b = model-predicted AUC (10 mg genormaliseerd) vs eGFR; Tabel S4 = mediaan/
  5e-95e percentiel AUC per eGFR-stratum (de brondata achter fig. 2).
- Gepubliceerde model (Tabel S3): 2-compartment, ADVAN4, Ka 2.833 h-1, CL/F 21.56
  L/h, Vc/F 73.47 L, Vp/F 168.5 L, Q/F 9.137 L/h, eGFR~CL 0.6374 (macht, eq. 1:
  CL/F = CLREF/F*(eGFR/eGFRREF)^0.6374 * (SEXM + SEX~CL), SEX~CL = -0.1441),
  BWT~Vc 0.5424 (eq. 2), BSV op Ka (FIX 8.52), CL/F (0.1136), Vp/F (0.8229).
- UGT1A9 maturatie literatuur:
  * Miyagi SJ et al. Drug Metab Dispos 2012;40:1321-7 (zdd1321): one-phase exponential
    association, eq. 1: Y = Y0 + (Ymax - Y0)*(1 - exp(-K*X)); gerapporteerd: Y0 = 0
    (geboorte), plateau Ymax = 27.9 +/- 1.5 nmol/min/mg, 90% van plateau op ~4.0
    maanden (= 1/3 jaar).
  * Farhan et al. J Clin Pharmacol 2024;64:1224 (jcph.2484): meta-analyse, eq. 1
    (Hill): fraction of adult = Fbirth + (AdultMax - Fbirth)*AGE^n/(Age50^n + AGE^n).
    Tabel 1, UGT1A9: IN VIVO (propofol CL): Fbirth=0.05, Age50=0.09 jaar, n=1.9,
    AdultMax=1.1. In vitro activity (Miyagi refit): Fbirth=0.12, Age50=0.23, n=6.7.
    In vitro abundance (LC-MS/MS): Fbirth=0.00, Age50=8.20, n=0.5 (outlier).
    Conclusies paper: in-vivo activity-based ontogeny functions beste PBPK-
    voorspelling (MSE 0.05 vs 0.104/0.129 in vitro).

1. BEPALING RENALE EN NIET-RENALE COMPONENT VAN CL/F
--------------------------------------------------------
1.1 Aannames
- Lineaire PK: AUC(tau) = Dosis/(F*CL)  =>  AUC evenredig met 1/CL.
- Parallelle eliminatie sommeren: CLtotaal = CLrenaal + CLniet-renaal (additief).
- Anker: CL(90) = 21.6 L/h (paper: CL/F = 21.56 L/h bij eGFRREF = 90 mL/min/1.73m2).
- eGFR is de CKD-EPI 2009 waarde genormaliseerd per 1.73 m2 (mL/min/1.73m2).

1.2 Stap-voor-stap afleiding
Stap 1. Tabel S4 strata (midpoint eGFR, mediaan AUC ng*h/mL, 5e, 95e):
  15-29: 22, 1020, 661, 1390
  30-44: 37,  819, 552, 1200
  45-59: 52,  666, 455, 1060
  60-89: 75,  528, 339,  866
  >=90: 100,  420, 275,  673
Stap 2. AUC(90) interpoleren m.b.v. de machtsvorm van de paper (AUC ~ eGFR^-0.6374):
  AUC(90) = 420 * (100/90)^0.6374 = 448.9 ng*h/mL.
Stap 3. Converteer per stratum naar absolute CL (L/h):
  CL(x) = 21.6 * AUC(90)/AUC(x)
  => (22, 9.51), (37, 11.85), (52, 14.57), (75, 18.38), (100, 23.10), plus anker (90, 21.60).
Stap 4. Least-squares fit van CL = CLnr + S*eGFR op de 6 punten:
  S    = cov(eGFR,CL)/var(eGFR) = 0.1765 L/h per mL/min/1.73m2
  CLnr = mean(CL) - S*mean(eGFR) = 5.44 L/h
  R^2  = 0.9985
  Resultaat:  CL/F = 5.44 + 0.1765*eGFR   (L/h)
  Renale fractie bij eGFR 90 = 0.1765*90/21.3 = ~74%.
  Dimensionezele renale slope: kr = S/0.06 = 2.94 mL/min per mL/min/1.73m2
  (0.06 = conversie mL/min -> L/h).
Stap 5. Controle vs anchors:
  CL(90) = 5.44 + 0.1765*90 = 21.33 ~ 21.6 L/h (anker) OK.
  AUC(0)/AUC(120) = CL(120)/CL(0) = (5.44+21.18)/5.44 = 4.89 (lineaire extrapolatie;
  fig. 2 loopt maar tot eGFR ~15-20, dus eGFR=0 is extrapolatie).

1.3 Alternatief dat is overwogen (twee-anker oplossing, NIET gekozen)
Ankers: CL(90)=21.6 en AUC(0)/AUC(120)=3.25 (losse lezing fig. 2 bij eGFR 0):
  (CLnr + 7.2*kr)/CLnr = 3.25  ->  kr = 0.3125*CLnr   (7.2 = 0.06*120)
  CLnr*(1 + 5.4*0.3125) = 21.6 -> CLnr = 8.04 L/h, kr = 2.51.
Verworpen omdat het maar 1 extrapolatiepunt forceert en de tussenliggende strata
slechter fit (bij eGFR 22: AUC-ratio 1.9 vs 2.4 geobserveerd).

1.4 Hybride variant (eveneens gefit, ter vergelijking)
CL = CLnr + a*eGFR^theta  (scipy curve_fit, bounds CLnr>=0, a>=0, 0<theta<=3):
  CL = 6.47 + 0.0914*eGFR^1.132, R^2 = 0.9991, AUC(0)/AUC(120) = 4.18.
Niet weerhouden: praktisch gelijk aan lineair, extra parameter.

1.5 Reproduceerbare code (python)
  import numpy as np
  egfr = np.array([22.,37.,52.,75.,100.]); auc = np.array([1020.,819.,666.,528.,420.])
  auc90 = 420*(100/90)**0.6374
  CL = 21.6*auc90/auc; CL = np.append(CL, 21.6); egfr = np.append(egfr, 90.)
  A = np.vstack([np.ones_like(egfr), egfr]).T
  (CLnr, S), *_ = np.linalg.lstsq(A, CL, rcond=None)
  # -> CLnr = 5.44, S = 0.1765

2. BEPALING UGT1A9-MATURATIE (leeftijd, relatief t.o.v. adult)
---------------------------------------------------------------
2.1 Gekozen functie (Farhan 2024, IN VIVO, UGT1A9, AdultMax genormaliseerd naar 1):
  MAT(AGE) = Fbirth + (1 - Fbirth)*AGE^n / (Age50^n + AGE^n)
           = 0.05 + 0.95*AGE^1.9 / (0.09^1.9 + AGE^1.9)      (AGE in jaren)
  AdultMax=1.1 (propofol-specifiek: kinderen > adult per kg) is bewust op 1 gezet
  omdat in dit model adult-CLnr het referentiepunt is.
2.2 Waarden: AGE=0 -> 0.05; 0.09 j (~1.1 mnd, Age50) -> 0.525; 0.23 j -> 0.863;
  0.33 j -> 0.926; 0.5 j -> 0.965; 1 j -> 0.990; >=2 j -> ~1.0; 11-17 j -> 1.0.
2.3 Vergelijking met Miyagi-exponentieel (K = -ln(0.10)/(4/12) = 6.9078 j^-1,
  MAT = 1 - exp(-6.9078*AGE)): 0 j -> 0; 0.09 -> 0.46; 0.33 -> 0.90; 1 j -> 0.999.
  => beide bronnen komen overeen: UGT1A9 practisch volwassen binnen ~1 jaar.
  Afweging: in-vivo curve heeft beste PBPK-voorspelaccuracy (Farhan MSE 0.05);
  abundance-curve (Age50 8.2 j) is outlier en NIET bruikbaar voor klaring.
2.4 Identificeerbaarheid: pediatrische studie NCT01525238 = 11-17 jaar, dus MAT = 1
  voor ALLE subjecten in de dataset. Age50/n/Fbirth zijn DAAROM GEFIXEERD
  (THETA 10-12 FIX), niet schatbaar; effect alleen bij extrapolatie < 2 jaar.
2.5 Toepassing: ALLEEN op de niet-renale component, als pure leeftijdsfactor,
  los van allometrie:
  CLnr = THETA(2) * MAT(AGE) * (BWT/84)^0.75
  Renale component schaalt via eGFR (bij kinderen al leeftijd/size-gecorrigeerd).

3. IMPLEMENTATIE IN NONMEM (dapa_popPK_run1.ctl)
-------------------------------------------------
$INPUT ID TIME AMT DV EVID MDV CMT EXCL EGFR BWT SEX PED CKD AGE   ; AGE in jaren
$PK
 MAT  = THETA(12) + (1.0D0-THETA(12))*AGE**THETA(11)/(THETA(10)**THETA(11)+AGE**THETA(11))
 CLNR = THETA(2)*MAT*((BWT/84.0D0)**0.75)
 TVCL = CLNR + THETA(6)*EGFR
 CL   = TVCL*(1.0D0 + THETA(7)*(1.0D0-SEX))*EXP(ETA(2))
 ...
$THETA
 (0.1 5.44 20.0)      ; 2 CLnr adult (L/h)  - startwaarde uit par. 1.2
 (0.01 0.177 1.0)     ; 6 renale slope S (L/h per mL/min/1.73m2)
 (-1.0 -0.1441 0.0)   ; 7 SEX~CL/F (vrouw)
 (0.01 0.09) FIX      ; 10 Age50 (jaar)
 (0.1 1.9) FIX        ; 11 Hill n
 (0.0 0.05 1.0) FIX   ; 12 Fbirth
Rest ongewijzigd t.o.v. paper: KA 2.833, Vc 73.47, Vp 168.5, Q 9.137, ALAG1 0.5 FIX
(placeholder, niet gepubliceerd), BSV Ka 8.52 FIX / CL 0.1136 / Vp 0.8229,
$SIGMA: prop adult 0.2845^2, prop peds 0.2519^2, add adult 0.1849^2, CKD log-add 0.5491.
Draaien: nmfe73 dapa_popPK_run1.ctl dapa_run1.lst (NONMEM 7.3, gfortran).

4. VALIDATIE (figuren)
----------------------
Fig 1 (fig1_auc_vs_egfr.png): AUC vs eGFR - Tabel S4 mediaan/5e-95e (zwart) vs
  paper-macht (blauw) vs additief CL=5.44+0.1765*eGFR (rood). Additief loopt door
  alle mediaanpunten: eGFR 22/37/52/75/100 -> geobs 1020/819/666/528/420 vs
  additief 1040/810/663/519/420 (macht: 1102/791/637/504/420).
Fig 2 (fig2_maturation.png): maturatiecurves vs leeftijd; paarse lijn = 11 j
  (pediatrische studie, MAT=1).
Fig 3 (fig3_cl_decomposition.png): CL/F = 5.44 + 0.1765*eGFR; eGFR 90 -> 21.3 L/h,
  eGFR 0 -> 5.44 L/h (alleen niet-renaal, eindige AUC).
Beperking: validatie op gepubliceerde stratum-samenvatting; striktere check =
NONMEM-run op individuele data (9715 obs) + GOF (DV vs IPRED, CWRES vs TIME) + pcVPC.

5. AANNAMES / OPEN PUNTEN
-------------------------
- eGFR per 1.73 m2 (CKD-EPI 2009). Absolute GFR = eGFR*BSA/1.73 (optioneel, kolom BSA).
- BWT-referentie 84 kg (mediaan dataset); anker CL(90)=21.6 L/h.
- ALAG1 = 0.5 h FIX is placeholder (lag time niet gepubliceerd).
- THETA 8 (BWT^0.5424 op Vc) ongebruikt omdat V2 lineair (BWT/84) schaalt.
- Data-cleaning upstream (BLQ, postdose > Cmax/2, outliers); in ctl via EXCL/MDV.
- Maturatieparameters gefixeerd uit literatuur (niet identificeerbaar in data 11-17 j).
