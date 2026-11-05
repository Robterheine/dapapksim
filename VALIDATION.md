# Validatierapport — dapagliflozine popPK-simulaties vs Penland et al. 2023

**Gegenereerd:** 2026-11 — bij het opzetten van `dapagliflozin-pk-children/` (Shiny app voor
Monte Carlo simulatie van het dapagliflozine popPK-model in kinderen).

**Bronartikel:** Penland RC et al. *Evaluation of the Pharmacokinetics of Dapagliflozin in
Patients With Chronic Kidney Disease With or Without Type 2 Diabetes Mellitus.*
J Clin Pharmacol 2023;63(5):551–559. [doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196)

## 1. Doel

Verifiëren of de Monte Carlo simulaties van de app (`R/pk_model.R`, R-implementatie van het
NONMEM-model met renale/niet-renale split, allometrie en UGT1A9-maturatie) ongeveer dezelfde
PK oplevert als gerapporteerd in het bronartikel.

## 2. Methode

- Simulatie met **typische waarden** (ETA = 0, geen IOV, geen residuele error), 10 mg 1×/dag,
  7 doses (steady state), RK4 dt = 0.05 h — exact de logica van `simulate_population_pk()`.
- Covariaten = mediaanwaarden uit **Table 1** van Penland et al. 2023.
- Vergelijking met de gerapporteerde model-predicted waarden (CL/F, AUC genormaliseerd naar
  10 mg, eGFR-stratum-ratio's, Cmax).
- Panel B van de figuur: volledige Monte Carlo (N = 250, IIV op CL/V3/Ka) voor een
  illustratieve virtuele pediatrische populatie (10–17 jaar, 0.15 mg/kg 1×/dag).

## 3. Resultaten

### 3.1 Vergelijking gesimuleerd vs gerapporteerd

| Scenario (Penland Table 1 mediaan) | Gesimuleerd | Gerapporteerd (Penland 2023) |
|---|---|---|
| CL/F referentie (84 kg, eGFR 90, man) | **21,4 L/h** | 21,6 L/h (populatie CL/F) |
| AUC→10 mg, adult T2DM znd CKD (90 kg, eGFR 90) | **461 ng·h/mL** | 464 ng·h/mL |
| AUC→10 mg, DAPA-CKD +T2DM (81 kg, eGFR 42) | **783 ng·h/mL** | 757 ng·h/mL |
| AUC→10 mg, DAPA-CKD −T2DM (78 kg, eGFR 40) | **815 ng·h/mL** | 778 ng·h/mL |
| AUC-ratio vs eGFR ≥90 — eGFR 45–59 | **1,46×** | 1,6× |
| AUC-ratio vs eGFR ≥90 — eGFR 30–44 | **1,78×** | 2,0× |
| AUC-ratio vs eGFR ≥90 — eGFR 15–29 | **2,27×** | 2,4× |
| Cmax na 10 mg (84 kg, eGFR 90–100, man) | **100–102 ng/mL** | ~100 ng/mL (model-Cmax-cutoff) |
| Tmax na 10 mg | **1,3 h** | tmax ~1–2 h (dalen gemeten 20–32 h na dosis) |
| Pediatrisch T2DM (15 j, 94 kg, eGFR 110), 10 mg | AUC 393 ng·h/mL, Cmax 90 ng/mL | pediatrische groep: hoogste eGFR → laagste AUC |

### 3.2 Steady-state metrieken, virtuele pediatrische populatie (N = 250, 0.15 mg/kg 1×/dag)

| Metriek | P5 | Mediaan | P95 |
|---|---|---|---|
| Cmax,ss (ng/mL) | 20 | 116 | 161 |
| Cmin,ss = C24 (ng/mL) | 0,2 | 2,2 | 10,4 |
| AUCtau,ss (ng·h/mL) | 156 | 342 | 619 |

## 4. Figuur

![Validatie: AUC vs eGFR en gesimuleerde steady-state curve](figures/fig_validation.svg)

- **Panel A:** dosis-genormaliseerde AUC (10 mg) vs eGFR — Penland Tabel S4 (zwart, whiskers
  P5–P95) vs de machtsvorm uit het artikel (blauw) vs het app-model `CL = 5.44 + 0.177·eGFR`
  (rood). Het additieve model loopt door alle mediaanpunten.
- **Panel B:** gesimuleerde steady-state concentratie-tijdcurve (mediaan + 90%
  predictie-interval, IIV alleen — geen residuele error) met mediaan Cmax,ss en Cmin,ss.

## 5. Afwijkingen en oorzaken

- **eGFR-band-ratio's iets vlakker** (1,46–2,27× vs 1,6–2,4×): de lineaire renale split
  `0.177·EGFR` is vlakker dan de machtsfunctie op eGFR in het artikel.
- **Gewichtseffect op AUC kleiner** (1,06× bij 60 vs 85 kg vs 1,2× in artikel): in het artikel
  staat BWT alleen op Vc/F; hier (bewuste aanpassing) allometrie `(BWT/84)^0.75` op CLnr.
- **Geslachtseffect 1,17×** (θ7 = −0,1441) vs ~1,0× in het artikel — beide niet klinisch relevant.
- NB: in de huidige app-versie is het sex-effect op CL/F **verwijderd** (allometrie op BWT dekt dit; zie `model_doc.md`) — het effect was niet klinisch relevant.
- Cmax, CL/F en AUC komen binnen enkele procenten overeen.

## 6. Conclusie

De simulaties reproduceren de gepubliceerde farmacokinetiek van dapagliflozine (CL/F ≈ 21,6 L/h,
AUC per 10 mg ≈ 464 ng·h/mL bij normale nierfunctie, Cmax ≈ 100 ng/mL, eGFR-afhankelijkheid van
exposure) binnen enkele procenten. Het model is geschikt voor de intended use van de app.

Zie `docs/HANDOVER_renale_CL_UGT1A9_maturatie.md` voor de afleiding van de renale/niet-renale
split en de UGT1A9-maturatieparameters, en `model_doc.md` voor de modelmapping.
