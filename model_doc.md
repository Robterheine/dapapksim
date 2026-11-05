# Model & assumpties

De app simuleert het **NONMEM populatiemodel van Penland et al.**
(J Clin Pharmacol 2023;63:551–559, [doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196))
met de volgende aanpassingen:

1. **Klaring gesplitst** in niet-renale klaring (gematureerd) en renale klaring: `CL = CLnr + THETA(6)*EGFR`.
2. **Allometrische schaling van alle parameters** t.o.v. `(BWT/84)` (84 kg = referentiewicht in de bronanalyse).
3. **Maturatie van de niet-renale klaring (UGT1A9)** volgens Farhan et al. 2024
   ([doi:10.1002/jcph.2484](https://doi.org/10.1002/jcph.2484), Table 1, propofol):
   `MAT = 0.05 + 0.95 * AGE^1.9 / (0.09^1.9 + AGE^1.9)`.

## Structuur (ADVAN4 TRANS4)

3 compartimenten: **depot (1)** met lagtijd `ALAG1 = 0.5 h` → eerste-orde absorptie
(`KA`) naar **centraal (2)** (volume `V2`, eliminatie `CL`) ↔ **perifeer (3)** (volume `V3`,
uitwisseling `Q`). Concentratie = A2/V2.

## Parameterwaarden (THETA)

| Parameter | THETA | Waarde | Betekenis |
|---|---|---|---|
| KA | 1 | 2.833 h⁻¹ | absorptie bij BWT=84; geschaald met (BWT/84)^−0.25 |
| CLnr | 2 | 5.44 L/h | niet-renale CL/F, volwassene bij BWT=84, MAT=1 |
| V2 | 3 | 73.47 L | centraal volume bij BWT=84; lineair in BWT/84 |
| V3 | 4 | 168.5 L | perifeer volume bij BWT=84; lineair in BWT/84 |
| Q | 5 | 9.137 L/h | intercompartimentale klaring; (BWT/84)^0.75 |
| Renale slope | 6 | 0.177 L/h per mL/min/1.73m² | renale component = slope × EGFR |
| Sex-effect | 7 | −0.1441 | **niet toegepast** — verwijderd omdat allometrie op BWT dit effect dekt (zie "Verschillen t.o.v. bronmodel") |
| ALAG1 | 9 | 0.5 h (fixed) | lagtijd depot; placeholder — niet gerapporteerd in de publicatie |
| TM50 UGT1A9 | 10 | 0.09 jaar | maturatie (Farhan 2024) |
| Hill UGT1A9 | 11 | 1.9 | maturatie (Farhan 2024) |
| Fbirth UGT1A9 | 12 | 0.05 | maturatie (Farhan 2024) |

Individuele parameters:

```
MAT  = 0.05 + 0.95 * AGE^1.9 / (0.09^1.9 + AGE^1.9)
CLnr = 5.44 * MAT * (BWT/84)^0.75
CL   = CLnr + 0.177*EGFR
KA   = 2.833 * (BWT/84)^(-0.25)
V2   = 73.47 * (BWT/84);  V3 = 168.5 * (BWT/84);  Q = 9.137 * (BWT/84)^0.75
```

## Variabiliteit

- **IIV** (log-normaal, exp(ETA)), conform `$OMEGA BLOCK(3)`:
  - ETA(1) op KA: ω² = **8.52 FIX** (vastgezet in de publicatie; altijd meegenomen).
  - ETA(2) op CL/F: ω² = 0.1136 (≈34% CV)
  - ETA(3) op Vp/F: ω² = 0.8229 (≈105% CV)
  - geen BSV op V2 en Q in dit model.
- **IOV** (schuifje, extensie t.o.v. het NONMEM model): per dosing occasion een eigen
  log-normale multiplier op de **relatieve biologische beschikbaarheid F1** (de dosis die
  het depot in gaat), met CV = ingestelde %. Gemodelleerd vanwege variabiliteit in de
  toegediende dosis bij kinderen (tablet opgelost, via sonde toegediend).
- **Residuele fout ($SIGMA / EPS) wordt expres NIET gesimuleerd** — alleen
  inter-individuele variabiliteit (en optionele IOV). De gesimuleerde curve per kind is
  de deterministische verwachte concentratie (A2/V2 gegeven de ETA's); de 90% band is
  daarmee een predictie-interval over individuen (P5–P95 van de gesimuleerde kinderen),
  niet over individuele observaties (geen VPC-style band).
- Numeriek: RK4 met dt = 0.05 h. Om numerieke instabiliteit door de zeer grote Ka-BSV
  (ω² = 8.52) te voorkomen, wordt individuele KA geclamped op [0.05, 25] h⁻¹; dit treft
  alleen de extreme staart waar absorptie nagenoeg momentaan is.

## Nierfunctie (EGFR) uit de copulapop-populatie — in mL/min/1.73m²

Het NONMEM covariaat is **eGFR in mL/min/1.73m²** (renale slope in L/h per mL/min/1.73m²),
dus er is **geen omrekening met BSA** nodig (BSA × eGFR/1.73 zou absoluut GFR in mL/min geven).

`copulapop::generate_population()` berekent zelf een kolom `EGFR` in exact deze eenheid,
maar **alleen voor datasets die creatinine bevatten** (in de app: NHANES, RUMC):

- **< 18 jaar**: bedside Schwartz — `eGFR = 0.413 × lengte(cm) / Scr(mg/dL)`
- **≥ 18 jaar**: CKD-EPI 2021 (race-free)

Bij datasets **zonder creatinine** (in de app: HSE) ontbreekt de EGFR-kolom; de app
valt dan terug op een geschatte **normale eGFR per 1.73m² per leeftijd**
(GFR-maturatie: `pmax(30, 120·AGE^1.7/(0.75^1.7 + AGE^1.7))`). De gebruikte bron staat in het
tabblad "Virtuele populate".

**Verificatie eenheid (mL/min/1.73m²):** copulapop berekent EGFR op exact één plek —
`calculate_egfr()` in `R/generate.R`, aangeroepen door zowel `generate_population()` als
`generate_population_from_copula()`: CKD-EPI 2021 race-free (≥18 j) en bedside Schwartz
`0.413 × lengte(cm) / Scr(mg/dL)` (<18 j). Beide formules zijn **per definitie genormaliseerd
per 1.73 m²** (mL/min/1.73m²) — dit staat ook in de copulapop-README ("EGFR … mL/min/1.73m²").
De app gebruikt `pop$EGFR` **direct** in `THETA(6)·EGFR` (0.177 L/h per mL/min/1.73m²) en doet
**nooit** een BSA-conversie naar absoluut GFR (mL/min). De fallback is eveneens per 1.73m².
Als permanente guard checkt `egfr_from_copula()` dat alle EGFR-waarden in het plausibele
bereik 5–300 mL/min/1.73m² liggen; het resultaat staat in het tabblad "Virtuele populate".

## Covariaten uit de virtuele populatie

- **BWT** = `WEIGHT` (kg) — allometrische referent 84 kg.
- **AGE** = `AGE` (jaar) — maturatie van niet-renale CL.
- **SEX** = `SEXISMALE` (1 = jongen, 0 = meisje) — alleen gebruikt in de eGFR-formules (via copulapop); niet meer in CL.
- **EGFR** = `EGFR` (mL/min/1.73m²) — zie boven.

## Steady state & metrieken

Dosis per gift = mg/kg × BWT (dit is de **keerdosis per toediening**, geen dagdosis); dosisinterval τ instelbaar in de app (0–24 h, stap 1 h, standaard 24 h; τ = 0 wordt afgevangen met een foutmelding). Random seed instelbaar in de app (standaard 42).
Na `n_doses` doses geldt het laatste interval als steady state. Per kind: **Cmax,ss**,
**Tmax**, **Cmin,ss (= C(τ))**, **AUCtau,ss** (trapeziumregel 0–τ), **Cavg,ss = AUCtau/τ**. In de app-tabel per metriek ook het **geometrisch gemiddelde** (naast mediaan en P5–P95).
In het tabblad "AUC & Cmax vs nierfunctie / leeftijd" worden AUCtau,ss en Cmax,ss
geplot vs eGFR en vs leeftijd met een LOESS-trendlijn.

## Verschillen t.o.v. het bronmodel (Penland 2023) — en waarom

- **eGFR-vergelijking**: het artikel gebruikte CKD-EPI **2009**; copulapop gebruikt CKD-EPI
  **2021** (≥18 j) en bedside Schwartz (<18 j). Het artikel rapporteert dat de overstap
  2009→2021 "minimal impact" had op de parameterschattingen.
- **eGFR-effect op CL**: in het artikel een machtsfunctie `(eGFR/eGFR_REF)^0.6374`; hier de
  lineaire split `CL = 5.44 + 0.177·eGFR` (afgeleid uit Tabel S4, R² = 0.9985 — zie
  `docs/HANDOVER_renale_CL_UGT1A9_maturatie.md`). Gevalideerd in `VALIDATION.md`
  (AUC-ratio's per eGFR-stratum iets vlakker dan in het artikel).
- **Maturatiefunctie (UGT1A9)** — welke, hoe, waarom:
  - *Functie*: Hill-functie `MAT(AGE) = Fbirth + (1 − Fbirth)·AGE^n / (Age50^n + AGE^n)`
    = `0.05 + 0.95·AGE^1.9 / (0.09^1.9 + AGE^1.9)` (AGE in jaren; THETA 10–12).
    Waarden: 0 j → 0,05; 0,09 j (Age50) → 0,525; 0,5 j → 0,965; ≥ ~2 j → ≈ 1.
  - *Implementatie*: alleen op de **niet-renale** klaring, als pure leeftijdsfactor los van
    allometrie — `CLnr = 5.44 · MAT(AGE) · (BWT/84)^0.75` (in de app: `individual_params()` in
    `R/pk_model.R`). De renale component schaalt via eGFR, die bij kinderen al
    leeftijd/size-gecorrigeerd is (bedside Schwartz).
  - *Waarom*: Farhan et al. 2024 (doi:10.1002/jcph.2484, Table 1, **in-vivo** UGT1A9 via
    propofol-CL): Fbirth = 0.05, Age50 = 0.09 jaar, n = 1.9 (AdultMax = 1.1 is genormaliseerd
    naar 1 omdat adult-CLnr hier het referentiepunt is). De in-vivo activity-based
    ontogeny gaf de beste PBPK-voorspelling (MSE 0.05 vs 0.104/0.129 in-vitro). Miyagi 2012
    (exponentieel, 90% plateau op ~4 maanden) komt hiermee overeen: UGT1A9 is praktisch
    volwassen binnen ~1 jaar. De abundance-curve (Age50 8.2 j) is een outlier en niet bruikbaar
    voor klaring.
  - *Identificeerbaarheid*: de pediatrische brondata betroffen 11–17 jaar (MAT = 1 voor alle
    subjecten), dus THETA 10–12 zijn **gefixeerd uit de literatuur**, niet schatbaar uit de
    data; het maturatie-effect is alleen zichtbaar bij extrapolatie < ~2 jaar.
- **IOV op F1** is een extensie t.o.v. het NONMEM-model. De getoonde steady-state curve is het
  **laatste doseringsinterval** en reflecteert per kind één F1-trekking; de 90% band bevat
  daarmee óók de occasion-to-occasion spreiding van die laatste dosis.
- **Sex-effect op CL/F (THETA(7) = −0.1441) verwijderd**: in het bronmodel correleerde dit met lichaamsgrootte; nu alle parameters allometrisch geschaald worden op (BWT/84), is dat effect niet meer nodig. Het effect was klein (~1.0–1.17× AUC) en volgens het artikel niet klinisch relevant.
- **ALAG1 = 0.5 h** is een placeholder (lagtijd niet gepubliceerd in het artikel).

## Verificatie tegen Penland et al. 2023 (doi:10.1002/jcph.2196)

Simulatie met typische waarden (ETA = 0, geen IOV), 10 mg 1×/dag, 7 doses (steady state),
met de mediaan-covariaten uit Table 1 van het artikel. Vergelijking met de gerapporteerde
model-predicted waarden:

| Scenario | Gesimuleerd | Gerapporteerd (Penland 2023) |
|---|---|---|
| CL/F referentie (84 kg, eGFR 90, man) | 21,4 L/h | 21,6 L/h (populatie CL/F) |
| AUC gen. naar 10 mg — adult T2DM znd CKD (90 kg, eGFR 90) | 461 ng·h/mL | 464 ng·h/mL |
| AUC gen. naar 10 mg — DAPA-CKD +T2DM (81 kg, eGFR 42) | 783 ng·h/mL | 757 ng·h/mL |
| AUC gen. naar 10 mg — DAPA-CKD −T2DM (78 kg, eGFR 40) | 815 ng·h/mL | 778 ng·h/mL |
| AUC-ratio vs eGFR ≥90 — eGFR 45–59 | 1,46× | 1,6× |
| AUC-ratio vs eGFR ≥90 — eGFR 30–44 | 1,78× | 2,0× |
| AUC-ratio vs eGFR ≥90 — eGFR 15–29 | 2,27× | 2,4× |
| Cmax na 10 mg (84 kg, eGFR 90–100, man) | 100–102 ng/mL | ~100 ng/mL (model-Cmax-cutoff) |
| Tmax na 10 mg | 1,3 h | tmax ~1–2 h; dalen gemeten 20–32 h na dosis |
| Pediatrisch T2DM (15 j, 94 kg, eGFR 110), 10 mg | AUC 393 ng·h/mL, Cmax 90 ng/mL | eGFR ≥90-groep pediatrisch: hoogste eGFR, dus laagste AUC |

Conclusie: CL/F, AUC (per 10 mg) en Cmax komen binnen enkele procenten overeen met de
gepubliceerde model-predicted waarden. Kleine afwijkingen zijn te verwachten:

- de lineaire renale split (`0.177·EGFR`) geeft iets vlakkere eGFR-band-ratio's dan de
  power-functie op eGFR in het artikel;
- de allometrische schaling op CLnr (gewenste aanpassing) maakt het gewichtseffect op AUC
  kleiner (1,06× bij 60 vs 85 kg) dan in het artikel (1,2×; daar staat BWT alleen op Vc/F);
- het geslachts­effect is 1,17× (θ7 = −0,1441) vs ~1,0× in het artikel (niet klinisch relevant).

## Bronnen

- Penland RC et al. J Clin Pharmacol 2023;63:551–559 — doi:10.1002/jcph.2196 (basismodel, 2-compartment + lagtijd; CL/F ≈ 21.6 L/h in CKD-populatie).
- Farhan H et al. 2024 — doi:10.1002/jcph.2484 (UGT1A9-ontogenie: TM50 = 0.09 j, Hill = 1.9, Fbirth = 0.05).
- Tirucherai GS et al. Diabetes Obes Metab 2016;18:678–684 (PK 10–17 jaar vergelijkbaar met volwassenen).
- Zwep LHM et al. Clin Pharmacol Ther 2024;115:795–804 — rationale copulapop (doi:10.1002/cpt.3099).
- [copulapop](https://github.com/Robterheine/copulapop) — virtuele populaties (vine-copulas).

## Beperkingen

- Vereenvoudiging: IOV op F1 is een extensie (zit niet in het NONMEM model); residuele
  fout niet gesimuleerd; Ka-BSV geclamped voor numerieke stabiliteit.
- Virtuele populaties zijn algemene-bevolkings-/ziekenhuispopulaties, geen T2DM-patiënten.
- Educatief/visualisatiemodel — geen doseringsadvies.
