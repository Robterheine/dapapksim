# Dapagliflozine PK in kinderen — interactieve Shiny app

Interactieve Shiny-app waarmee de farmacokinetiek van dapagliflozine in kinderen
wordt gevisualiseerd met Monte Carlo simulaties van een populatie-PK model.
De virtuele kinderen worden gegenereerd met het
[copulapop](https://github.com/Robterheine/copulapop) package (vine-copulas).

## Functionaliteit

- **PK-curve bij steady state**: mediaan + 90% predictie-interval (P5–P95) over de
  virtuele populatie, optioneel individuele curves, mediaan van dag 1, log-y-as.
- **Steady-state metrieken**: Cmax,ss, Tmax, Cmin,ss (= C24), AUCtau,ss en Cavg,ss —
  als populatiesamenvatting (mediaan [P5–P95]) én als verdelingen per kind.
- **Dosis in mg/kg** (schuifje), omgerekend per kind op basis van gewicht.
- **Leeftijdsrange** instelbaar; de virtuele populatie wordt Fermi-gefilterd gegenereerd
  via `copulapop::generate_population(..., age_range = ...)`.
- **Inter-occasion variabiliteit (IOV)**: schuifje (0–100% CV) dat per dosing occasion
  een eigen log-normale multiplier op de klaring (CL) legt.
- Keuze uit copulapop datasets (NHANES, HSE, RUMC),
  geslacht, aantal kinderen, aantal doses tot steady state, random seed.
- Download van curve-samenvatting en per-kind metrieken als CSV.

## Installatie

```r
install.packages(c("shiny", "ggplot2", "markdown"))
install.packages(c("VineCopula", "data.table", "jsonlite", "Matrix", "MASS"))
# copulapop zit nog niet op CRAN:
devtools::install_github("Robterheine/copulapop")
```

Vereist R ≥ 4.1 en ggplot2 ≥ 3.4.

## Runnen

```r
shiny::runApp("/pad/naar/dapagliflozin-pk-app")
```

Als de copulapop JSON-data niet automatisch worden gevonden, vul dan in de app het
veld "copulapop data-pad" in (bv. `system.file("data", package = "copulapop")`).

## Bestandsstructuur

```
dapagliflozin-pk-app/
├── app.R              # Shiny UI + server (Nederlands)
├── R/
│   ├── pk_model.R     # popPK model: parameters, IIV/IOV, RK4 Monte Carlo simulatie, metrieken
│   └── population.R   # wrapper rond copulapop::generate_population()
├── model_doc.md       # modelbeschrijving, assumpties en bronnen (tab in de app)
├── www/style.css
└── README.md
```

## Model (kort)

R-implementatie van het NONMEM popPK-model van Penland et al. 2023
([doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196)) — ADVAN4/TRANS4:
2-compartment met eerste-orde absorptie + lagtijd (0,5 h). Klaring gesplitst in
niet-renaal (`5.44·MAT·(BWT/84)^0.75`) en renaal (`0.177·EGFR`, EGFR in
mL/min/1.73m² — direct uit de copulapop-kolom `EGFR`, Schwartz <18 j / CKD-EPI 2021
≥18 j; fallback per leeftijd als de dataset geen creatinine heeft). Maturatie van
niet-renale CL (UGT1A9) volgens Farhan 2024 (doi:10.1002/jcph.2484). Allometrie op
alle parameters t.o.v. BWT/84. IIV op KA (ω²=8.52, fixed), CL/F en Vp/F; optionele
IOV op CL per dosis. Zie `model_doc.md` in de app voor de volledige mapping.

Let op: dit is een educatief/visualisatiemodel — geen gevalideerd doseringsadvies.
