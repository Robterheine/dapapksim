# Evaluatie van de app — klinisch farmacoloog & UI-designer

Datum: 2026-11-04. Review op basis van de code en documentatie in deze repo
(`app.R`, `R/`, `model_doc.md`, `VALIDATION.md`, `docs/HANDOVER_renale_CL_UGT1A9_maturatie.md`).

---

## 1. Review: klinisch farmacoloog

**Wat sterk is**
- Modelmapping is 1-op-1 met het NONMEM-model (ADVAN4/TRANS4, THETA-tabel in `model_doc.md`),
  inclusief de aanpassingen: renale/niet-renale split `CL = 5.44 + 0.177·eGFR`, allometrie
  op BWT/84, UGT1A9-maturatie (Farhan 2024), sexeffect −0.1441.
- **eGFR-eenheid is hard geverifieerd**: copulapop berekent EGFR op één plek
  (`calculate_egfr()` in `R/generate.R`) met CKD-EPI 2021 (≥18 j) en bedside Schwartz
  (<18 j) — beide per definitie mL/min/1.73m². De app gebruikt `pop$EGFR` direct in
  `THETA(6)·EGFR`, **zonder BSA-conversie** naar absoluut mL/min, plus een runtime
  sanity-check (bereik 5–300) die in het tabblad "Virtuele populatie" staat.
- Geen residuele error gesimuleerd (alleen IIV + optionele IOV op F1) — staat expliciet
  in de app én in `model_doc.md`; de 90% band is daarmee een PI over individuen.
- Validatie tegen het bronartikel is gedocumenteerd en navolgbaar (`VALIDATION.md` +
  figuur): CL/F 21,4 vs 21,6 L/h; AUC10 461 vs 464 ng·h/mL; Cmax ~100 ng/mL.
- IOV op F1 (toegediende dosis) is klinisch relevant voor sonde-toediening en uitgelegd.
- Bronnen zijn te achterhalen: DOI's in `model_doc.md`, afleiding in het handover-document,
  validatietabel in `VALIDATION.md`.

**Aandachtspunten (opgelost in deze commit)**
- ALAG1 = 0.5 h was alleen in het handover-document als "placeholder (niet gepubliceerd)"
  gemarkeerd → nu ook in de THETA-tabel van `model_doc.md`.
- Verschil in eGFR-vergelijking (artikel: CKD-EPI 2009; app: CKD-EPI 2021 via copulapop)
  en machts- vs lineair eGFR-effect waren alleen impliciet → nu een eigen sectie
  "Verschillen t.o.v. het bronmodel" in `model_doc.md`.
- Identificeerbaarheid van de maturatie (brondata 11–17 j ⇒ MAT=1; THETA 10–12 gefixeerd
  uit literatuur; effect alleen < ~2 j) → nu gedocumenteerd.
- Bij IOV op F1 is de getoonde steady-state curve het laatste interval (één F1-trekking
  per kind) → nu verduidelijkt in `model_doc.md`.

**Conclusie farmacoloog:** modeltechnisch verantwoord en transparant gedocumenteerd;
bruisbaar als educatief/visualisatiemiddel (geen doseringsadvies — staat er ook).

---

## 2. Review: UI-designer

**Wat sterk is**
- Duidelijke hiërarchie: prominente "Simuleer"-knop bovenin de sidebar; logische groepering
  (Simulatie-instellingen → Weergave → Gemeten data).
- Tabs volgen de workflow: curve → metrieken → exposure vs covariaten → populatie → model.
- Kleurblind-vriendelijke paletten (Okabe-Ito) voor geslacht en gemessen data.
- Foutafhandeling is zichtbaar (rode alert + notificatie), incl. duidelijke messages bij
  ontbrekende kolommen in de overlay-data.
- Publication-ready export zit op de logische plek (onder de curve) met zinnige defaults
  (7×5 inch, 600 dpi PNG).

**Aandachtspunten (opgelost in deze commit)**
- Er was geen introductietekst: een nieuwe gebruiker zag alleen een titel en tabbladen →
  nu een korte intro boven de tabs (wat doet de app, waar beginnen, waar is de onderbouwing).
- De curve-pagina was tamelijk vol (plot + helptext + exportpaneel + CSV-knoppen) —
  acceptabel; de helptext staat direct onder de plot waar hij hoort.

**Minor suggesties (niet doorgevoerd, optioneel)**
- Overweeg op termijn bslib/theming voor een modernere look; huidig plain-Shiny is
  bewust licht in dependencies.
- Bij hele smalle vensters klapt de sidebar boven de content — standaard Shiny-gedrag.

**Conclusie UI-designer:** intuïtief en bruikbaar; degegroepeerde sidebar en duidelijke
call-to-action werken; documentatie is vanaf de app zelf bereikbaar (tabblad
"Model & assumpties").

---

## 3. Samenvattend antwoord op de vier vragen

| Vraag | Antwoord |
|---|---|
| Bruikbaar? | Ja — één klik op "Simuleer" geeft direct curve + metrieken; fouten worden direct zichtbaar gemeld |
| Intuïtief? | Ja — logische sidebar-groepen, tabs in workflow-volgorde, intro toegevoegd |
| Onderbouwingen duidelijk? | Ja — `model_doc.md` in de app: THETA/OMEGA-mapping, assumpties, verschillen t.o.v. bronmodel, eGFR-eenheid-verificatie |
| Bronnen te achterhalen? | Ja — DOI's (Penland 2023 jcph.2196, Farhan 2024 jcph.2484, Zwep 2024 cpt.3099), afleiding in `docs/HANDOVER_...md`, validatie in `VALIDATION.md` + figuur |
