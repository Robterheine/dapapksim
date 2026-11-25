# NONMEM-reproductiepakket (`dapapksim`)

Deze map bevat een zelfstandige **NONMEM 7.x control stream** (`ADVAN4 TRANS4`) en een **1-kind voorbeelddataset** waarmee het volledige fysiologische 6-pijler populatie-PK/PD model inclusief de pediatrische hartfalen-uitbreidingen 1-op-1 in NONMEM kan worden gereproduceerd.

## Bestanden in deze map

| Bestand | Beschrijving |
|---|---|
| `dapa_ped_hf_sim.ctl` | Volledige NONMEM control stream (`$SUBROUTINES ADVAN4 TRANS4`) met exact dezelfde fysiologische en hartfalen-vergelijkingen als `R/pk_model.R` en `R/hf_pkpd.R`. |
| `example_child_hf.csv` | NONMEM-invoerbestand voor **1 voorbeeldkind (`ID = 1`)**: zuigeling van `9 maanden` (`0.75 j`) met pediatrisch hartfalen, `7` dagelijkse giften van `1.20 mg` (`0.15 mg/kg` op drooggewicht) en rijke bemonstering op **Dag 1** (`0–24 h`) én **Dag 7** (`144–168 h`). |
| `expected_r_reference.csv` | Exacte referentie-output berekend door `dapapksim` (`R/pk_model.R`) op alle 32 observatietijdstippen ter directe controle van de NONMEM-tabel `sdtab_dapa_ped_hf.tab`. |
| `generate_and_verify.R` | R-script dat `example_child_hf.csv` en `expected_r_reference.csv` genereert én — zodra NONMEM lokaal gedraaid heeft — automatisch `sdtab_dapa_ped_hf.tab` vergelijkt met R. |

---

## Klinisch profiel van het voorbeeldkind (`ID = 1`)

Het voorbeeldkind activeert gelijktijdig alle 6 fysiologische modelpijlers én de hartfalen-pathofysiologie zodat elke vergelijking toetsbaar is:

- **Demografie & overvulling:** jongen (`SEX = 1`), leeftijd `AGE = 0.75 jaar` (`9 mnd`), lengte `HT = 0.71 m`, actueel gewicht `BWT = 8.40 kg` bij `EDEMA = 5%` overvulling $\rightarrow$ euvolemisch drooggewicht `WDRY = 8.00 kg` (`VEDEM = 0.40 L` extra in `V2`).
- **Hartfalen-pathofysiologie:** cardiorenaal GFR-behoud `RFFAC = 0.80` (`80%`), hepatische stuwing `HEPFAC = 0.90` (`90%`), hypoalbuminemie `ALBIN = 32.0 g/L`, enterale biologische beschikbaarheid `F1REL = 1.00`.
- **UGT1A9-ontogenie:** `MATOPT = 1` (*Farhan et al. 2024* in-vivo standaard; zet op `2` voor *Badee 2019* of `3` voor *Simcyp s1*).
- **Dosering:** `0.15 mg/kg` op drooggewicht = **`1.20 mg` 1× per 24 uur** (`AMT = 1.2`, `ADDL = 6`, `II = 24`, `TAU = 24`).
  *(Opmerking: expliciete meervoudige dosering met `ADDL = 6, II = 24` vermijdt de bekende NONMEM-eigenschap waarbij `SS = 1` i.c.m. `ALAG1 = 0.5 h` op `0 ≤ TIME < 0.5 h` nog nul voorspelt vóór het lagtijdstip.)*

### Verwachte individuele PK/PD-parameters (`ETA = 0`)

| Parameter | Verwachte waarde | Eenheid |
|---|---|---|
| `WDRY` / `FFM` | `8.000000` / `6.456421` | `kg` |
| `ALBU` / `FU` | `32.000000` / `0.114556` | `g/L` / `–` |
| `NGFR` / `EGFR` / `RF` | `97.55082` / `78.04065` / `0.800000` | `mL/min/1.73m²` / `–` |
| `MAT` (Farhan 2024) | `0.983385` | `–` |
| `CL` (`CLHEP` + `CLRUGT` + `CLFILT`) | `4.564798` (`1.212178` + `3.125517` + `0.227104`) | `L/h` |
| `V2` / `V3` / `Q` / `KA` | `8.458278` / `21.37610` / `1.726980` / `2.649688` | `L` / `L` / `L/h` / `h⁻¹` |
| `AUCU24` (steady-state) | `30.11452` | `ng·h/mL` |
| `ETGF` / `EHCT` | `85.27518` / `59.26361` | `%` |

### Verwachte concentraties & luminale SGLT2-bezetting (selectie)

| Dag | `TIME (h)` | `TAD (h)` | `CTOT (ng/mL)` | `CU (ng/mL)` | `SGLT2 (%)` |
|---|---|---|---|---|---|
| 1 | `0.5` | `0.5` | `0.000000` | `0.000000` | `0.00000` |
| 1 | `1.0` | `1.0` | `83.623496` | `9.579538` | `97.95491` |
| 1 | `2.0` | `2.0` | `61.783627` | `7.077659` | `97.25186` |
| 1 | `12.0` | `12.0` | `2.575936` | `0.295088` | `59.60313` |
| 1 | `24.0` | `24.0` | `1.288932` | `0.147654` | `42.47159` |
| 7 | `144.0` | `0.0` | `1.731832` | `0.198391` | `49.79807` |
| 7 | `144.5` | `0.5` | `1.683354` | `0.192838` | `49.08838` |
| 7 | `145.0` | `1.0` | `85.259729` | `9.766978` | `97.99337` |
| 7 | `146.0` | `2.0` | `63.329540` | `7.254752` | `97.31715` |
| 7 | `156.0` | `12.0` | `3.452095` | `0.395457` | `66.41234` |
| 7 | `168.0` | `24.0` | `1.732194` | `0.198433` | `49.80330` |

---

## Uitvoeren in NONMEM en automatische verificatie

1. **Draai de control stream vanuit de map `nonmem/`:**
   ```bash
   cd nonmem
   nmfe75 dapa_ped_hf_sim.ctl dapa_ped_hf_sim.lst
   # of via PsN:
   execute dapa_ped_hf_sim.ctl
   ```
2. **Vergelijk `sdtab_dapa_ped_hf.tab` automatisch met R (`dapapksim`):**
   ```bash
   Rscript generate_and_verify.R
   ```
   Het script leest `sdtab_dapa_ped_hf.tab` in en rapporteert per parameter en concentratiekolom het maximale absolute verschil `max |NONMEM - R|` (`~10⁻⁶`, uitsluitend begrensd door de tabel-afronding `s1PE14.6`).

---

## Kolomdefinities (`$INPUT`) & populatie-Monte Carlo modus

- **Eenheidsschaling:** `AMT` en `DOSEMG` staan in **`mg`**; door `S2 = V2 / 1000.0` levert `F = A(2) / S2` direct concentraties in **`ng/mL` (`µg/L`)**.
- **Ontbrekende labwaarden (`-99`):**
  - `ALBIN = -99`: berekent leeftijdsafhankelijk serumalbumine volgens *McNamara & Alcorn (2002)*.
  - `EGFRIN = -99` en `CREAT = -99`: berekent gezonde leeftijds-`nGFR` volgens *Rhodin et al. (2009)*.
  - `EGFRIN = -99` en `CREAT > 0` (`µmol/L`): gebruikt *Pierce 2021 CKiD U25* (`1–<18 j`) of *CKD-EPI 2021* (`≥18 j`), begrensd op `[5, 180] mL/min/1.73m²`.
- **Monte Carlo populatiesimulatie met `IIV`:**
  In `dapa_ped_hf_sim.ctl` staat `$OMEGA` standaard op `0 FIX` voor exacte deterministische verificatie. Vervang de drie `0 FIX` regels onder `$OMEGA` door `0.36 FIX` (`KA`), `0.1136 FIX` (`CL`) en `0.8229 FIX` (`V3`) om inter-individuele variabiliteit in NONMEM mee te simuleren.
