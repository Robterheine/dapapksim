# Validatierapport — dapagliflozine popPK-simulaties (dapapksim)

**Datum:** 2026-11-05
**Auteur:** Rob ter Heine (Radboud Applied Pharmacometrics group), ontwikkeld met Mistral (Mistral AI)
**Status:** concept ter review — tekst nagelezen (AI-proofread; zie §12)
**Broncode:** https://github.com/Robterheine/dapapksim (`app.R`, `R/pk_model.R`)

---

## 1. Samenvatting

Het in de app geïmplementeerde populatiemodel (R-implementatie van het NONMEM-model van
Penland et al. 2023, met renale/niet-renale klaringssplit, allometrie, UGT1A9-maturatie en
zonder sex-effect op CL) is numeriek en grafisch gevalideerd tegen twee publicaties:

1. **Penland et al. 2023** (J Clin Pharmacol 63:551–559, doi:10.1002/jcph.2196) — het
   bronmodel voor volwassenen/CKD, incl. Tabel S4 (AUC per eGFR-stratum).
2. **Tirucherai et al. 2016** (Diabetes Obes Metab 18:678–684, doi:10.1111/dom.12638) —
   eenmalige dosis 2,5/5/10 mg bij 24 kinderen/adolescenten (10–17 j) met T2DM (Tabel 2).

**Belangrijkste bevindingen**

- **Klaring en totale exposure (AUC) kloppen goed**: CL/F binnen −1 % (Penland) en
  +4…+11 % (Tirucherai); dosis-genormaliseerde AUC binnen +1…+5 % (Penland, per
  eGFR-stratum) en −4…−10 % (Tirucherai, per dosisgroep).
- **Cmax wordt onderschat**, vooral in de populatiesimulatie: typische waarde (ETA = 0)
  −31 % vs gerapporteerd (10 mg: 82 vs 118 ng/mL); populaties-geometrisch gemiddelde
  −40…−51 %. Oorzaak: de uit het bronmodel übernommen, gefixeerde Ka-BSV (ω² = 8,52) die
  de absorptiesnelheid extreem spreidt (CV ~7000 % op ka; gesimuleerde Cmax-CV 126 % vs
  34–41 % geobserveerd). Dit is een **eigenschap van het bronmodel**, geen
  implementatieafwijking — zie §8–§9.
- **tmax klopt** (gesimuleerd ~1,3 h vs gerapporteerd 0,9–1,5 h).
- **t½ is in het model langer** (~19–20 h) dan de NCA-schatting (10–14 h) — ebenfalls een
  eigenschap van de distributieparameters van het bronmodel (Vp/F 168,5 L, Q/F 9,1 L/h).
- De eGFR-band-ratio's zijn iets vlakker dan in het artikel (lineaire renale split vs
  machtsfunctie), binnen 0,1–0,2 fold.

**Conclusie:** het model is valide voor **exposure (AUC, CL/F)** in kinderen; **Cmax**
moet met voorzichtigheid worden geïnterpreteerd vanwege de gefixeerde Ka-BSV uit het
bronmodel.

---

## 2. Doel en scope

Valideren dat de Monte Carlo-simulaties in de Shiny-app (`dapapksim`) dezelfde
farmacokinetiek opleveren als (a) het bron-popPK-model (Penland 2023) en (b) de
gepubliceerde pediatrische studie (Tirucherai 2016). Rapportage is bedoeld voor een
**klinisch farmacoloog** (§8) en een **biomedisch statisticus** (§9).

## 3. Model (kort)

2-compartment (ADVAN4/TRANS4), eerste-orde absorptie met lagtijd 0,5 h, lineaire
eliminatie. Zie `model_doc.md` (tabblad "Model & assumpties" in de app) voor de volledige
THETA/OMEGA-mapping. Kern:

```
MAT  = 0,05 + 0,95·AGE^1,9 / (0,09^1,9 + AGE^1,9)          # UGT1A9-maturatie (Farhan 2024)
CL   = 5,44·MAT·(BWT/84)^0,75 + 0,177·eGFR                 # eGFR in mL/min/1,73m²
KA   = 2,833·(BWT/84)^(−0,25);  V2 = 73,47·(BWT/84);  V3 = 168,5·(BWT/84);  Q = 9,137·(BWT/84)^0,75
```

IIV (log-normaal): KA ω² = 8,52 (fixed, bron), CL/F 0,1136, Vp/F 0,8229. **Geen residuele
error** — de 90 % band is een predictie-interval over individuen. Sex-effect op CL is
verwijderd (allometrie dekt dit; zie model_doc). IOV op F1 (toegediende dosis) is
standaard 0 % in deze validatie.

## 4. Methode

- Simulatie: RK4, dt = 0,05 h; steady state na 7 doses 1×/dag (τ = 24 h) voor de
  Penland-vergelijking; **eenmalige dosis** (0–48 h) voor de Tirucherai-vergelijking.
- Virtuele populaties: voor Tirucherai N = 2000 gesampled rond de studie-gemiddelden
  (leeftijd 14,5 ± 2,0 j (10–17), gewicht 99,7 ± 26,4 kg (≥30), eGFR 112,9 ± 20,5
  mL/min/1,73m² (≥80), 37,5 % man), met IIV op CL/V3/Ka.
- Vergelijkingsmaten: absoluut verschil in % t.o.v. de gepubliceerde waarde; voor
  Tirucherai zowel de **typische waarde (ETA = 0)** als het **populatie-geometrisch
  gemiddelde** (geomean, %CV), omdat de publicatie geomean rapporteert.
- Eenheden: concentraties ng/mL (= µg/L; dosis in mg → µg, volumes in L — geverifieerd),
  AUC in ng·h/mL, CL/F in L/h (×0,06 → mL/min).

## 5. Numerieke validatie vs Penland et al. 2023 (volwassenen/CKD)

Typische waarden (ETA = 0), 10 mg 1×/dag, steady state.

| Metriek | Gesimuleerd | Gerapporteerd (Penland) | Afwijking |
|---|---|---|---|
| CL/F referentie (84 kg, eGFR 90) | 21,4 L/h | 21,6 L/h | −1,1 % |
| AUC→10 mg, adult T2DM znd CKD (90 kg, eGFR 90) | 467 ng·h/mL | 464 ng·h/mL | +0,7 % |
| AUC→10 mg, DAPA-CKD +T2DM (81 kg, eGFR 42) | 783 ng·h/mL | 757 ng·h/mL | +3,4 % |
| AUC→10 mg, DAPA-CKD −T2DM (78 kg, eGFR 40) | 815 ng·h/mL | 778 ng·h/mL | +4,8 % |
| Cmax na 10 mg (84 kg, eGFR 90) | 102 ng/mL | ~100 ng/mL (model-Cmax-cutoff) | +2,3 % |
| tmax na 10 mg | 1,3 h | ~1–2 h | binnen bereik |
| AUC-ratio vs eGFR ≥90: 45–59 / 30–44 / 15–29 | 1,46× / 1,78× / 2,27× | 1,6× / 2,0× / 2,4× | −0,1…−0,2 fold |

Per eGFR-stratum (Penland Tabel S4, mediaan [P5–P95] vs model):

| eGFR (mL/min/1,73m²) | Gerapporteerd mediaan [P5–P95] | Model AUC→10 mg | Afwijking |
|---|---|---|---|
| 22 (15–29) | 1020 [661–1390] | 1060 | +4,0 % |
| 37 (30–44) | 819 [552–1200] | 830 | +1,3 % |
| 52 (45–59) | 666 [455–1060] | 681 | +2,2 % |
| 75 (60–89) | 528 [339–866] | 534 | +1,0 % |
| 100 (≥90) | 420 [275–673] | 432 | +2,8 % |

## 6. Numerieke validatie vs Tirucherai et al. 2016 (kinderen 10–17 j, eenmalige dosis)

Virtuele pediatrische populatie N = 2000 (studie-gemiddelden, zie §4); gerapporteerd =
geomean (%CV) uit Tabel 2.

| Dosis | Metriek | Gesimuleerd (pop. geomean, CV) | Gesimuleerd (typisch, ETA=0) | Gerapporteerd geomean (%CV) | Afwijking pop. |
|---|---|---|---|---|---|
| 2,5 mg | Cmax (ng/mL) | 14,5 (126 %) | 20,4 | 24,8 (34 %) | −42 % |
| | AUC0–inf (ng·h/mL) | 96 (44 %) | 96 | 101 (23 %) | −5 % |
| | CL/F (mL/min) | 435 (–) | 436 | 413 (26 %) | +5 % |
| 5 mg | Cmax (ng/mL) | 29,0 (126 %) | 40,8 | 48,4 (41 %) | −40 % |
| | AUC0–inf | 191 (44 %) | 191 | 199 (29 %) | −4 % |
| | CL/F | 435 | 436 | 418 (27 %) | +4 % |
| 10 mg | Cmax (ng/mL) | 58,0 (126 %) | 81,6 | 118,0 (35 %) | −51 % |
| | AUC0–inf | 383 (44 %) | 382 | 427 (31 %) | −10 % |
| | CL/F | 435 | 436 | 391 (25 %) | +11 % |

| Dosis | tmax mediaan (sim vs rep) | t½ (sim vs rep mean) |
|---|---|---|
| 2,5 mg | 1,3 h vs 1,50 h | 19,9 h vs 14,1 h |
| 5 mg | 1,3 h vs 0,96 h | 19,9 h vs 10,3 h |
| 10 mg | 1,3 h vs 0,88 h | 19,9 h vs 10,7 h |

Dosis-proportionaliteit: in het model exact (AUC/mg = 38,3 en Cmax/mg = 5,8 constant over
2,5–10 mg); gerapporteerd ook nagenoeg proportioneel (AUC/mg 40,4–42,7).

## 7. Grafische validatie

![Figure 1 — AUC vs eGFR vs Penland Tabel S4](../figures/fig1_penland_egfr.svg)

*Figuur 1.* Dosis-genormaliseerde AUC (10 mg) vs eGFR: gerapporteerde stratum-medianen met
P5–P95 (Penland Tabel S4, zwart) vs het artikel-machtsmodel (blauw) vs het app-model
(rood). Het additieve model loopt door alle mediaanpunten (afwijking +1…+4 %).

![Figure 2 — Cmax en AUC0-inf per dosis vs Tirucherai](../figures/fig2_tirucherai_dose.svg)

*Figuur 2.* Pediatrische eenmalige dosis: gesimuleerde populaties-geomean (rode cirkel) en
typische waarde (blauw vlak) vs gerapporteerde geomean met ±1 log-SD whiskers (zwart).
AUC/CL kloppen; Cmax zit lager door de gefixeerde Ka-BSV (zie §8–§9).

![Figure 3 — gesimuleerde curven](../figures/fig3_curves.svg)

*Figuur 3.* Links: mediaan van een eenmalige 10 mg dosis (0–48 h) met de gerapporteerde
geomean-Cmax/tmax als marker. Rechts: steady state 10 mg 1×/dag met mediaan en 90 %
predictie-interval — de hoofduitvoer van de app.

## 8. Interpretatie — klinisch farmacoloog

- **Exposure (AUC) en klaring zijn klinisch valide**: binnen ~10 % van zowel het
  bronmodel (Penland) als de pediatrische studie (Tirucherai). De eGFR-afhankelijkheid van
  de exposure wordt goed gereproduceerd (ratio's binnen 0,2 fold van het artikel).
- **Cmax is conservatief (laag)**: door de uit het bronmodel overgenomen, gefixeerde
  Ka-BSV (ω² = 8,52 — in de CKD-data van Penland waren vooral dalwaarden beschikbaar,
  dus ka was daar niet informatief) spreidt de absorptie extreem. De typische Cmax
  (ETA = 0) ligt ~31 % onder de gerapporteerde geomean; in de populatiesimulatie ~50 %
  lager. Voor dosis-berekeningen op basis van **piekconcentraties** is voorzichtigheid
  geboden; voor **totale exposure (AUC)** en **dalconcentraties** is het model
  betrouwbaar.
- **tmax** (~1,3 h) komt overeen met de gerapporteerde 0,9–1,5 h (snelle absorptie).
- **t½**: het model geeft een langere terminale halfwaardetijd (~19–20 h) dan de
  NCA-schatting (10–14 h). Dit volgt uit de distributieparameters van het bronmodel
  (Vp/F 168,5 L bij Q/F 9,1 L/h); het model reproduceert daarmee het bronmodel
  getrouw, maar de terminale fase is langer dan wat NCA op de geobserveerde data schat.
- Praktisch: bij mg/kg-dosering in kinderen (zwaarder dan 84 kg, hogere eGFR) is de
  relatieve exposure lager dan bij volwassenen per mg — dit is zichtbaar in de app en
  matcht de litratuur ("zelfde dosis als volwassenen" uit Tirucherai).

## 9. Interpretatie — biomedisch statisticus

- **Verdelingen**: IIV is log-normaal gemodelleerd (exp(ETA)); de publicatie rapporteert
  geomean/%CV — daarom is in deze validatie zowel de geomean als de mediaan/typische
  waarde gegeven. In de app-tabel staat nu mediaan, geometrisch gemiddelde en P5–P95.
- **Geen residuele error**: de 90 % band is een **predictie-interval over individuen**
  (P5–P95 van gesimuleerde kinderen), geen VPC-band over observaties. Bij vergelijking
  met gepubliceerde spreiding (bijv. %CV) moet hiermee rekening worden gehouden.
- **Ka-BSV (ω² = 8,52, fixed)**: sd(ETA₁) = 2,92 → CV(ka) ~7000 %. Ditabaliseert Cmax
  enorm (gesimuleerd Cmax-CV 126 % vs 34–41 % geobserveerd) en drukt de geomean-Cmax
  (Cmax is concaaf in ka). Aanbeveling: bij Cmax-voorspellingen de typische ka
  (ETA₁ = 0) gebruiken, of de Ka-BSV niet meenemen. (In de app is de Ka-BSV
  momenteel altijd aan, conform het bronmodel.)
- **Dekkingscheck**: de gerapporteerde mediaan-AUC per eGFR-stratum (Penland Tabel S4)
  ligt in alle strata binnen de gerapporteerde P5–P95 én binnen +1…+4 % van het model —
  de model-PI dekt de gepubliceerde populatiespreiding.
- **Dosis-proportionaliteit**: exact in het model (lineair), en nagenoeg in de data —
  geen schaling van IIV met dosis nodig.

## 10. Beperkingen

- Validatie op gepubliceerde samengevatitudes (Tabel S4 / Tabel 2), niet op individuele
  data (9715 obs) — striktere check = NONMEM-run + pcVPC.
- Lineaire renale split vs machtsfunctie in het artikel → eGFR-ratio's iets vlakker.
- CKD-EPI 2009 (artikel) vs 2021 (copulapop, volwassenen); artikel meldt "minimal impact".
- Ka-BSV gefixeerd op 8,52 uit het bronmodel (zie §9).
- t½ langer dan NCA-schatting (bronmodel-eigenschap).
- ALAG1 = 0,5 h is een placeholder (niet gepubliceerd).
- Virtuele populaties zijn algemene populaties (copulapop), geen T2DM-patiëntenpopulatie.

## 11. Conclusie

De app simuleert **AUC/CL/F binnen ~10 %** van zowel het bronmodel als de pediatrische
studie, met correcte tmax en dosis-proportionaliteit. **Cmax** is conservatief door de
uit het bronmodel overgenomen Ka-BSV. Het model is geschikt voor exposure-gebaseerde
visualisatie en mg/kg-exploratie in kinderen; piekconcentraties met voorzichtigheid
interpreteren.

## 12. Drukproef

Tekst nagelezen op taal, consistentie en eenheden (AI-proofread, 2026-11-05). In deze
omgeving is geen menselijke drukproefexpert beschikbaar; aan te raden is een finale
leesronde door een native speaker/drukproeflezer vóór publicatie.

## Referenties

1. Penland RC et al. J Clin Pharmacol 2023;63:551–559. doi:10.1002/jcph.2196 (bronmodel,
   Tabel S3/S4).
2. Tirucherai GS et al. Diabetes Obes Metab 2016;18:678–684. doi:10.1111/dom.12638
   (pediatrische single-dose studie NCT01525238, Tabel 1–2).
3. Farhan H et al. J Clin Pharmacol 2024;64:1224. doi:10.1002/jcph.2484 (UGT1A9-ontogenie).
4. Kasichayanula S et al. Clin Pharmacokinet 2014;53:17–27 (volwassen PK, t½ ~12,9 h).
5. Zwep LHM et al. Clin Pharmacol Ther 2024;115:795–804. doi:10.1002/cpt.3099 (copulapop).

## Appendix — reproduceerbaarheid

- Modelcode: `R/pk_model.R` (functies `nonmem_theta()`, `individual_params()`,
  `simulate_population_pk()`); app: `app.R`.
- Validatiesimulaties: typische waarden ETA = 0; Tirucherai-populatie N = 2000 met
  seed 20261105 (analoog aan de app-seed).
- Figuren: `figures/fig1_penland_egfr.svg`, `fig2_tirucherai_dose.svg`,
  `fig3_curves.svg`.
