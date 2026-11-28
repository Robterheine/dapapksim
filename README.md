# Dapagliflozine PK & fysiologische extrapolatie — interactieve Shiny app

Interactieve Shiny-app waarmee de farmacokinetiek van dapagliflozine van pasgeborenen
tot volwassenen (inclusief obesitas en chronische nierschade) wordt gesimuleerd en
gevisualiseerd met Monte Carlo simulaties van een fysiologisch geëxtrapoleerd
populatie-PK model. Virtuele patiëntenpopulaties worden gegenereerd met het
[copulapop](https://github.com/Robterheine/copulapop) package (vine-copulas).

## Functionaliteit

- **PK-curve (eenmalige gift of steady state)**: mediaan + 90% predictie-interval
  (P5–P95) over de virtuele populatie, optioneel individuele curves (deterministisch
  geselecteerd per seed), mediaan van dag 1, en logaritmische y-as.
- **Volledige doseringsvrijheid zonder plafond**: keuze tussen **op lichaamsgewicht
  (`mg/kg`)** of **vaste dosis (`mg`)** per gift, zonder kunstmatig dosisplafond —
  bedoeld voor vrije dosis- en leeftijds-extrapolatie.
- **Vrije leeftijdsrange (pasgeborene t/m oudere)**: fijne numerieke invoer voor
  minimum- en maximumleeftijd (bv. `0–0.08` jaar voor neonaten, `0.08–2` jaar voor
  zuigelingen/peuters, `10–17` jaar voor adolescenten, of `18–100` jaar voor volwassenen).
- **Totale én ongebonden expositiemetrieken + PK-curve**: `Cmax,ss`, `Tmax`, `Cmin,ss (= C(τ))`,
  `AUCtau,ss` (totaal), **`AUCu,ss` en `AUCu,24h` (ongebonden, gecorrigeerd voor albumine-ontogenie en hypoalbuminemie)**,
  `Cavg,ss`, terminale halfwaardetijd `t1/2,β`, én directe weergave van de **ongebonden plasmacurve `Cu(t)` met luminale SGLT2-drempels**.
- **Pediatrisch hartfalen: pathofysiologie, duale PK–PD & dosisvoorspelling (`R/hf_pkpd.R`)**:
  - Klinische hartfalenprofielen en fijnafstelling voor **cardiorenale hypoperfusie (`RF`)**, **hepatische stuwing** (bv. biventriculair falen / Fontan), **extracellulair oedeem boven drooggewicht**, **hypoalbuminemie / PLE** en **relatieve biologische beschikbaarheid `F1`** bij sonde-/suspensietoediening.
  - **Probability of Target Attainment (PTA)** scan op steady-state `AUCu,24h` + directe gevoeligheidsvergelijking over alle drie UGT1A9-ontogeniefuncties en een leeftijdsgestratificeerde referentietabel (`0 mnd – 18 jaar`).
  - **Semi-mechanistische PK–PD koppeling**: onderscheid tussen proximaal tubulaire luminale SGLT2-receptorbezetting (`%`), `GFR`-gekoppelde effectieve renale osmotische/natriuretische PD-drive (`%`) en systemische cardiorenaal `Emax`-responsen.
- **Drie gevalideerde UGT1A9-maturatiefuncties** (werkzaam op zowel hepatische als
  renale tubulaire UGT1A9-glucuronidering):
  1. **Farhan et al. 2024** (*in-vivo* UGT1A9, standaard)
  2. **Badee et al. 2019 / Jo et al. 2021 scenario 2**
  3. **Simcyp standaardontogenie / Jo et al. 2021 scenario 1**
- **Keuze absorptie-variabiliteit (`Ka`-BSV)**: gekalibreerde klinische variabiliteit
  (`ω² = 0.36`, standaard; reproduceert klinische `Cmax`-spreiding van 34–41%), het ruwe
  gefixeerde bronmodel (`ω² = 8.52 FIX` uit Penland 2023), of geen `Ka`-BSV (`ω² = 0`).
- **Gemeten PK-data als overlay**: upload een CSV/TSV of plak direct in het tekstveld
  (komma-, puntkomma- of tab-gescheiden; ondersteunt ook Nederlandse decimaalkomma's bij `;`).
- **Publication-ready figuur**: download de curve als PNG/PDF/SVG/TIFF met instelbare
  afmetingen en resolutie (met Engelstalige as- en titelteksten voor publicatie).
- **Inter-occasion variabiliteit (IOV) op `F1`**: mean-gecentreerd log-normaal schuifje
  (`0–100% CV`, zodat `E[F1] = 1.0`) voor variabiliteit in toegediende sonde-/dispersiedosis.
- **Tweetalig (Engels / Nederlands, standaard Engels)**: schakelaar rechtsboven in de app
  (`English | Nederlands`) waarmee de volledige interface, keuzelijsten, figuren, tabellen,
  waarschuwingen en in-app documentatie direct tussen Engels en Nederlands wisselen zonder
  herstart van de simulatie.
- **Snelle exacte analytische solver**: gesloten-vorm `ADVAN4` toestandspropagator
  (machineprecisie, ~50× sneller dan numerieke RK4-integratie, met exacte behandeling van
  absorptielagtijd `ALAG1` en optionele directe analytische steady state).

## Installatie

```r
install.packages(c("shiny", "ggplot2", "markdown"))
install.packages(c("VineCopula", "data.table", "jsonlite", "Matrix", "MASS"))
# copulapop installeren vanaf GitHub:
devtools::install_github("Robterheine/copulapop")
```

Vereist R ≥ 4.1 en ggplot2 ≥ 3.4.

## Runnen

```r
shiny::runApp("/pad/naar/dapapksim")
```

Als de copulapop JSON-data niet automatisch worden gevonden, stel dan vóór het starten
in R in: `copulapop::set_data_path(system.file("data", package = "copulapop"))`.

## Bestandsstructuur

```
dapapksim/
├── app.R                          # Shiny UI + server (tweetalig EN/NL; standaard Engels)
├── R/
│   ├── pk_model.R                 # Fysiologisch 6-pijler popPK-model + analytische ADVAN4-solver + PD
│   ├── hf_pkpd.R                  # Pediatrisch hartfalen: pathofysiologie, PTA & duale PK-PD module
│   └── population.R               # Adaptieve wrapper rond copulapop::generate_population()
├── scripts/
│   └── validate_model.R           # Reproduceerbare R-validatiesuite & SVG-figuurgenerator
├── hf_pkpd_rationale_en.md        # PK/PD rationale voor pediatrisch hartfalen (Engels, standaard)
├── hf_pkpd_rationale.md           # PK/PD rationale voor pediatrisch hartfalen (Nederlands)
├── model_rationale_en.md          # Stapsgewijze fysiologische afleiding & kalibratie (Engels)
├── model_rationale.md             # Stapsgewijze fysiologische afleiding & kalibratie (Nederlands)
├── model_doc_en.md                # Volledige wiskundige formules & R/NONMEM-code (Engels)
├── model_doc.md                   # Volledige wiskundige formules & R/NONMEM-code (Nederlands)
├── VALIDATION_en.md               # Samenvattende validatie (Engels)
├── VALIDATION.md                  # Samenvattende validatie (Nederlands)
├── about_en.md                    # Over deze app (Engels)
├── about.md                       # Over deze app (Nederlands)
├── docs/
│   ├── VALIDATION_REPORT.md       # Volledig validatierapport (Engels: Penland, Kasichayanula, Tirucherai, Jo PBPK)
│   ├── VALIDATION_REPORT.html     # Zelfstandig HTML-rapport met ingebedde SVG-figuren
│   ├── VALIDATION_REPORT.pdf      # PDF-versie van het validatierapport
│   ├── validation_results.json    # Exacte numerieke output van scripts/validate_model.R
│   ├── EVALUATION.md              # Expert-review (klinisch farmacoloog + UI-designer)
│   └── HANDOVER_renale_CL_UGT1A9_maturatie.md
├── figures/                       # Gegenereerde validatiefiguren (SVG)
├── www/style.css
└── README.md
```

## Fysiologisch extrapolatiemodel (kort)

Uitbreiding van het NONMEM popPK-model van Penland et al. 2023
([doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196), `ADVAN4/TRANS4`) voor
betrouwbare extrapolatie naar neonaten, zuigelingen, kinderen met obesitas en CKD:

1. **Lichaamssamenstelling (`FFM` & `NFM`)**: vetvrije massa volgens Al-Sallami et al. 2015
   (`< 18 j`) en Janmahasatian et al. 2005 (`≥ 18 j`), met *normal fat mass* voor `CL/Q`
   (`Ffat = 0.25`, macht `0.75`) en `V2` (`Ffat = 0.20`, lineair).
2. **Plasma-eiwitbinding (`fu`)**: ontogenie van serumalbumine volgens McNamara & Alcorn
   2002 (`fu_adult = 0.086`), toegepast op `CL/F`, lipofiel `V3/F` en `AUCu,ss`.
3. **Nierfunctie (`eGFR` & `nGFR`)**: `CKD-EPI 2021` (`≥ 18 j`), `Pierce et al. 2021 CKiD U25`
   (`1–<18 j` bij datasets met creatinine zoals `RUMC`), en `Rhodin et al. 2009` PMA-rijping
   bij `< 1 j` en datasets zonder creatinine (`NHANES`, `HSE`).
4. **Drie-paden klaring**: hepatische UGT1A9 (`5.44 L/h`), renale tubulaire UGT1A9
   (`15.78 L/h`, geschaald met zowel `MAT(AGE)` als relatieve nierfunctie `RF = eGFR / nGFR`),
   en onveranderde glomerulaire excretie (`0.011 × eGFR × BSA/2`).
5. **Ontwikkelings-GI-absorptie**: rijping van neonatale maaglediging op `KA` + gekalibreerde
   klinische `Ka`-BSV (`ω² = 0.36`).

Zie [`model_doc.md`](model_doc.md) en [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md)
voor alle vergelijkingen en kwantitatieve validatie tegen klinische volwassen en pediatrische PK-studies.

Let op: dit is een onderzoeks- en educatief simulatiemodel — geen gevalideerd klinisch doseringsadvies.

## Auteur

Rob ter Heine — Radboud Applied Pharmacometrics group, Radboudumc
([website](https://www.radboudumc.nl/en/research/research-groups/applied-pharmacometrics)),
r.terheine@radboudumc.nl.
Laatste versie: [github.com/Robterheine/dapapksim](https://github.com/Robterheine/dapapksim).
