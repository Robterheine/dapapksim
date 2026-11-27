# Volledige wiskundige, algoritmische & NONMEM-reproduceerbare modelspecificatie

Deze sub-pagina bevat de **exacte wiskundige vergelijkingen, parameterwaarden, randvoorwaarden, stochastische trekkingen en gesloten-vorm `ADVAN4`-solver** van **dapapksim** (`R/pk_model.R`, `R/population.R` en `R/hf_pkpd.R`). Zie sub-tabblad **1 (*Stapsgewijze Fysiologische Afleiding & Kalibratie*)** voor de inhoudelijke klinisch-farmacologische uitleg en kalibratie-achtergrond. Naast het zelfstandige R-referentieblok in §6 bevat de map **`nonmem/`** in de repository een kant-en-klare NONMEM 7.x control stream (`nonmem/dapa_ped_hf_sim.ctl`) met een 1-kind voorbeelddataset (`nonmem/example_child_hf.csv`) voor directe reproductie in NONMEM.

---

## 1. Overzicht en onderbouwing van de wijzigingen t.o.v. Penland et al. 2023

Het uitgangspunt is het twee-compartimentele orale NONMEM popPK-model van **Penland et al. 2023**
(*J Clin Pharmacol* 2023;63:551–559, [doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196)),
ontwikkeld op volwassen T2DM-, DAPA-CKD- en adolescente T2DM-data (`11–17 jaar`). Omdat deze
applicatie bedoeld is voor **vrije extrapolatie naar elke leeftijd en elke dosis** (waaronder
pasgeborenen, zuigelingen, ernstige obesitas en chronische nierschade), zijn zes empirische
beperkingen van het bronmodel vervangen door mechanistische fysiologische schaling:

| Onderdeel | Empirisch bron-/eerste model | Huidig fysiologisch model (`dapapksim`) | Fysiologische & klinische onderbouwing |
|---|---|---|---|
| **1. Lichaamsgrootte & obesitas** | Totale lichaamsmassa `(BWT/84)^0.75` en `BWT/84` | **Vetvrije massa (`FFM`) & Normal Fat Mass (`NFM`)** voor `CL`, `Q` (`Ffat = 0.25`) en `V2` (`Ffat = 0.20`); totaal `BWT` voor lipofiel `V3` | Vetweefsel draagt niet evenredig bij aan lever-/niermassa. In obese adolescenten (Tirucherai 2016, gem. `99,7 kg`) overschatte `BWT/84` de klaring en het centrale volume; `FFM`/`NFM` dekt tevens het geslachtsverschil zodat `THETA(7)` vervalt. |
| **2. Eiwitbinding (`fu`)** | Impliciet constant (`fu` onafhankelijk van leeftijd) | **Albumine-ontogenie (McNamara & Alcorn 2002)** rond `fu_adult = 0.086` (`91,4%` gebonden bij `44 g/L`); rapportage van **`AUCu,ss`** | Pasgeborenen hebben lager serumalbumine (`~34 g/L`), waardoor `fu` ~27% hoger ligt (`0.109`). Voor een laag-extractie farmacon (`ER < 0.05`) schalen totale `CL/F` en `V3/F` lineair met `fu/fu_adult`, terwijl farmacologisch actieve ongebonden blootstelling (`AUCu,ss`) afhangt van intrinsieke klaring. |
| **3. Nierfunctie (`eGFR` & `nGFR`)** | Vaste bedside Schwartz (`0.413`) op `0–18 j`; onjuiste aanname dat `NHANES` creatinine bevat; ad-hoc fallback `>= 30` | **`CKD-EPI 2021`** (`≥ 18 j`), **`Pierce 2021 CKiD U25`** (`1–<18 j` bij `CREAT`), en **`Rhodin 2009` PMA-rijping** (`< 1 j` en datasets zonder `CREAT`), begrensd op `[5, 180]` | Cross-sectionele `0–18 j` copula-correlatie mist maternale creatinineoverdracht bij pasgeborenen en gaf onfysiologisch hoge neonatale Schwartz-waarden (`> 110`). `CKiD U25` corrigeert leeftijds-/geslachtsbias tussen `1` en `18 j`. |
| **4. Renale klaring & renaal UGT1A9** | Volledige `eGFR`-afhankelijke klaring (`~74%` van `CL`) gemodelleerd als passieve filtratie zonder enzymrijping | **Drie-paden decompositie:** hepatische UGT1A9 (`5.44 L/h`) + **renale tubulaire UGT1A9** (`15.78 L/h`, geschaald met `MAT` én `RF = eGFR/nGFR`) + onveranderde filtratie (`0.011 × eGFR × BSA/2`) | Slechts `< 2%` van dapagliflozine wordt onveranderd renaal uitgescheiden; het steile `eGFR`-effect komt door proximale tubulaire UGT1A9-glucuronidering tot D3OG (Kasichayanula 2014, Jo 2021). Door `RF = eGFR/nGFR` rijpt renaal UGT1A9 mee met `MAT` zonder dubbeltelling van lage neonatale GFR, maar daalt klaring wél proportioneel bij CKD (`RF < 1`). |
| **5. UGT1A9-ontogenie** | Twee opties, waarbij `badee2019` per abuis door `1.305` werd gedeeld (`MAT(18 j) = 0.76`) | **Drie gevalideerde opties (`MAT(volwassen) = 1.0`):** `Farhan 2024` (standaard), gecorrigeerde `Badee 2019` (Jo s2) en `Simcyp s1` (Jo s1) | Dekt het volledige literatuurspectrum van snelle *in-vivo* rijping tot conservatieve langzame *in-vitro* ontogenie zonder discontinuïteit rond `18–20 jaar`. |
| **6. GI-absorptie, `Ka`-BSV & ODE-solver** | `(BWT/84)^-0.25` (snellere `KA` in baby's!), `ω²_Ka = 8.52 FIX`, en pure-R RK4 met `ifelse`- en lagtijd-bugs | **Neonatale maagledigingsrijping `f_GI(AGE)`**, gekalibreerde **`ω²_Ka = 0.36`** (standaard), mean-gecentreerde `IOV`, en **exacte analytische `ADVAN4`-solver** | Maaglediging is bij pasgeborenen trager, niet sneller. `ω²_Ka = 8.52` was een sparse-dalspiegelartefact in Penland (`Cmax`-bias `-51%`, CV `126%`); `ω²_Ka = 0.36` reproduceert klinische `Cmax` en `35%` totale CV. De analytische solver is exact tot machineprecisie en ~50× sneller. |

---

## 2. Compartimentele structuur (`ADVAN4 TRANS4`)

Het structurele model bestaat uit drie compartimenten met hoeveelheden \((A_1, A_2, A_3)\) in **microgram (`µg`)**, tijd \(t\) in **uur (`h`)**, volumes in **liter (`L`)** en klaringen in **`L/h`**:

1. **Absorptiedepot (\(A_1\), darm):** orale dosis \(D\) (in `µg` \(= 1000 \times \text{dosis}_{\text{mg}}\)) komt het depot binnen na absorptielagtijd \(\text{ALAG}_1 = 0.5\text{ h}\) en wordt met eerste-orde snelheidsconstante \(K_a\) (`h⁻¹`) opgenomen.
2. **Centraal compartiment (\(A_2\), volume \(V_2\)):** ontvangt \(K_a A_1\), elimineert lineair met totale schijnbare klaring \(CL\) (\(k_e = CL / V_2\)), en wisselt uit met het perifere compartiment via \(Q\) (\(k_{12} = Q / V_2\), \(k_{21} = Q / V_3\)).
3. **Perifeer compartiment (\(A_3\), volume \(V_3\)):** weefseldistributiecompartiment.

De plasmaconcentratie is \(C(t) = A_2(t) / V_2\) in **`µg/L = ng/mL`**.

$$\frac{dA_1}{dt} = -K_a A_1, \qquad \frac{dA_2}{dt} = K_a A_1 - (k_e + k_{12})A_2 + k_{21}A_3, \qquad \frac{dA_3}{dt} = k_{12}A_2 - k_{21}A_3$$

---

## 3. Volledige parameterspecificatie & covariaatvergelijkingen

### 3.1 Vaste populatieparameters (`nonmem_theta()`)

| Symbool in code | Waarde | Eenheid | Beschrijving & herkomst |
|---|---|---|---|
| `ka` | `2.833` | `h⁻¹` | Volwassen absorptiesnelheidsconstante `THETA(1)` (Penland 2023) |
| `cl_hep` | `5.44` | `L/h` | Hepatische UGT1A9-klaring (`CL/F`) bij ref-grootte, `MAT = 1`, `fu = 0.086` (afgeleid uit Penland Tabel S4 intercept bij `eGFR → 0`) |
| `cl_ren_ugt` | `15.78` | `L/h` | Renale tubulaire UGT1A9-klaring (`CL/F`) bij ref-grootte, `MAT = 1`, `RF = 1`, `fu = 0.086` (zie afleiding §3.5) |
| `cl_filt_s` | `0.011` | `L/h` per `mL/min/1.73m²` | Onveranderde glomerulaire/renale excretieslope bij `BSA_ref = 2.0 m²` (`0.99 L/h` bij `eGFR = 90`, ~4,6% van totale `CL/F`; Kasichayanula 2014) |
| `v2` | `73.47` | `L` | Centraal verdelingsvolume `Vc/F` (`THETA(3)`) bij `NFM_V2,ref = 61.68 kg` |
| `v3` | `168.5` | `L` | Perifeer verdelingsvolume `Vp/F` (`THETA(4)`) bij `BWT_ref = 84.0 kg`, `fu = 0.086` |
| `q` | `9.137` | `L/h` | Intercompartimentale klaring `Q/F` (`THETA(5)`) bij `NFM_CL,ref = 63.08 kg` |
| `alag` | `0.5` | `h` | Absorptielagtijd `ALAG1` (gefixeerd) |
| `bwt_ref` / `bsa_ref` | `84.0` / `2.0` | `kg` / `m²` | Referentiegewicht en -BSA van de Penland-bronpopulatie (`84 kg`, `1.71 m`) |
| `ffm_ref` | `56.1` | `kg` | Geslachtsgewogen referentie-FFM van de Penland-populatie (`~65%` man, `84 kg`, `BMI 28.7`) |
| `nfm_cl_ref` / `nfm_v2_ref` | `63.08` / `61.68` | `kg` | Referentie-NFM: `56.1 + 0.25×(84−56.1) = 63.08` en `56.1 + 0.20×(84−56.1) = 61.68` |
| `fu_adult` / `alb_adult` | `0.086` / `44.0` | `–` / `g/L` | Volwassen ongebonden fractie (`91.4%` gebonden) en referentie-serumalbumine |
| `gi_f0` / `gi_t50` | `0.45` / `0.10` | `–` / `jaar` | Neonatale fractie maagledigingssnelheid bij geboorte en rijpingshalfwaardetijd (`~5 weken`) |
| `tm50` / `hill` / `fbirth` | `0.09` / `1.9` / `0.05` | `jaar` / `–` / `–` | Farhan et al. 2024 *in-vivo* UGT1A9-maturatieparameters |

### 3.2 Covariaten uit de virtuele populatie (`copulapop`)

Per virtueel individu worden gebruikt:
- `BWT = WEIGHT` (`kg`), `AGE` (`jaar`, begrensd op \(\ge 0\)), `SEXISMALE` (`1` = man/jongen, `0` = vrouw/meisje).
- `HEIGHT` in `m` (als een invoertabel `HEIGHT > 3` bevat, wordt automatisch gedeeld door `100` voor `cm → m`).
- `BMI = BWT / max(HEIGHT^2, 0.09)` (\(\text{kg/m}^2\)).
- `BSA = sqrt((HEIGHT * 100 * BWT) / 3600)` (\(\text{m}^2\), Mosteller), tenzij `BSA` al expliciet in de tabel staat.
- *Fallback voor synthetische validatietabellen zonder `HEIGHT`:* als `BSA` bekend is, wordt `HEIGHT` geïnverteerd uit Mosteller (`min(2.15, max(0.35, (BSA^2 * 3600 / max(BWT, 0.5)) / 100))`); anders wordt een gladde groeicurve van `0.50 m` naar `1.76 m` (man) / `1.63 m` (vrouw) gebruikt (`estimate_height_m()`).

### 3.3 Lichaamssamenstelling (`FFM` en `NFM`)

Volgens **Janmahasatian et al. 2005** (\(\ge 18\text{ j}\)) en **Al-Sallami et al. 2015** (\(< 18\text{ j}\)):

$$\text{FFM}_{\text{adult}} = \begin{cases} \dfrac{9270 \cdot \text{BWT}}{6680 + 216 \cdot \text{BMI}} & \text{als SEXISMALE} = 1 \\[10pt] \dfrac{9270 \cdot \text{BWT}}{8780 + 244 \cdot \text{BMI}} & \text{als SEXISMALE} = 0 \end{cases}$$

$$M_{\text{FFM}}(\text{AGE}) = \begin{cases} 1.0 & \text{als AGE} \ge 18 \\[6pt] 0.88 + \dfrac{0.12}{1 + (\text{AGE}/13.4)^{-12.7}} & \text{als AGE} < 18 \text{ en man (}0.88\text{ bij AGE}=0\text{)} \\[10pt] 1.11 - \dfrac{0.11}{1 + (\text{AGE}/7.1)^{-1.1}} & \text{als AGE} < 18 \text{ en vrouw (}1.11\text{ bij AGE}=0\text{)} \end{cases}$$

$$\text{FFM} = \min\!\left(\text{FFM}_{\text{adult}} \cdot M_{\text{FFM}}(\text{AGE}), \; 0.97 \cdot \text{BWT}\right)$$

$$\text{NFM}_{CL} = \text{FFM} + 0.25(\text{BWT} - \text{FFM}), \qquad F_{\text{size},CL} = \left(\frac{\text{NFM}_{CL}}{63.08}\right)^{0.75}$$

$$\text{NFM}_{V_2} = \text{FFM} + 0.20(\text{BWT} - \text{FFM}), \qquad F_{\text{size},V_2} = \frac{\text{NFM}_{V_2}}{61.68}$$

### 3.4 Ontogenie van plasma-eiwitbinding (`fu_mcnamara()`)

Volgens **McNamara & Alcorn 2002** (met \(\text{days} = \max(1.0, \; 365.25 \times \text{AGE})\) wanneer `ALBUMIN` niet gemeten is):

$$[\text{Alb}](\text{AGE}) = \begin{cases} 44.0\text{ g/L} & \text{als AGE} \ge 18 \\ \min\!\left(44.0, \; 1.1287 \ln(\text{days}) + 33.746\right)\text{ g/L} & \text{als AGE} < 18 \end{cases}$$

$$f_u(\text{AGE}) = \frac{1}{1 + \dfrac{[\text{Alb}](\text{AGE})}{44.0} \cdot \dfrac{1 - 0.086}{0.086}}, \qquad f_{u,\text{ratio}} = \frac{f_u(\text{AGE})}{0.086}$$

### 3.5 Nierfunctie-ontogenie (`ngfr_rhodin()`) en individuele `eGFR` (`egfr_from_copula()`)

1. **Verwachte gezonde GFR naar postmenstruele leeftijd (`nGFR` in `mL/min/1.73m²`; Rhodin et al. 2009 + veroudering boven 35 j):**
   - Geschatte zwangerschapsduur bij geboorte (`GA` in weken) voor kleine pasgeborenen:
     $$\text{GA} = \begin{cases} \min\!\left(40, \max(24, \; 22 + 7.2 \cdot \text{BWT})\right) & \text{als AGE} < 0.1\text{ j en BWT} < 2.5\text{ kg} \\ 40 & \text{anders} \end{cases}$$
   - Postmenstruele leeftijd \(\text{PMA} = \text{GA} + 52.1775 \cdot \text{AGE}\) (weken).
   - $$\text{nGFR}(\text{AGE}, \text{BWT}) = 115.0 \cdot \frac{\text{PMA}^{3.4}}{47.7^{3.4} + \text{PMA}^{3.4}} \cdot 0.9925^{\max(0, \; \text{AGE} - 35)}$$

2. **Individuele `eGFR` (`mL/min/1.73m²`):**
   - **Indien `CREAT` (µmol/L) aanwezig is in de dataset (`RUMC`):** zij \(\text{SCr} = \max(\text{CREAT} / 88.4, \; 0.05)\) in `mg/dL`:
     - **Volwassenen (\(\text{AGE} \ge 18\text{ j}\), CKD-EPI 2021 race-free):**
       met \(\kappa = 0.9\) (man) / \(0.7\) (vrouw), \(\alpha = -0.302\) (man) / \(-0.241\) (vrouw), \(f_{\text{sex}} = 1.0\) (man) / \(1.012\) (vrouw):
       $$\text{eGFR} = 142 \cdot \min(\text{SCr}/\kappa, 1)^{\alpha} \cdot \max(\text{SCr}/\kappa, 1)^{-1.200} \cdot 0.9938^{\text{AGE}} \cdot f_{\text{sex}}$$
     - **Kinderen en adolescenten (\(1 \le \text{AGE} < 18\text{ j}\), Pierce et al. 2021 CKiD U25):**
       $$k = \begin{cases} 39.0 \cdot 1.008^{\text{AGE}-12} & \text{man, AGE} < 12 \\ 39.0 \cdot 1.045^{\text{AGE}-12} & \text{man, } 12 \le \text{AGE} < 18 \\ 36.1 \cdot 1.008^{\text{AGE}-12} & \text{vrouw, AGE} < 12 \\ 36.1 \cdot 1.023^{\text{AGE}-12} & \text{vrouw, } 12 \le \text{AGE} < 18 \end{cases}, \qquad \text{eGFR} = \frac{k \cdot \text{HEIGHT}_m}{\text{SCr}}$$
     - **Zuigelingen (\(\text{AGE} < 1\text{ j}\)):** \(\text{eGFR} = \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\) met \(\eta_{\text{GFR}} \sim \mathrm{N}(0, 0.215^2)\).
     - Alle waarden worden begrensd op **`[5, 180] mL/min/1.73m²`**.
   - **Indien `CREAT` ontbreekt (`NHANES`, `HSE`):**
     $$\text{eGFR} = \min\!\left(180, \; \max\!\left(5, \; \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\right)\right), \qquad \eta_{\text{GFR}} \sim \mathrm{N}(0, 0.215^2)$$

3. **Relatieve nierfunctie (`RF`) en afleiding van de 3-paden klaring:**
   $$\text{RF} = \min\!\left(1.75, \; \max\!\left(0, \; \frac{\text{eGFR}}{\max(\text{nGFR}(\text{AGE}, \text{BWT}), 1.0)}\right)\right)$$
   *Afleiding anker:* In Penland et al. 2023 bedraagt populatie-`CL/F = 21.56 L/h` bij `BWT = 84 kg`, `BSA = 2.0 m²`, `eGFR = 90` en mediane leeftijd `62 jaar` (`nGFR(62) = 93.84`, dus `RF = 90 / 93.84 = 0.9591`). Met `CL_hep = 5.44 L/h` en onveranderde filtratie `0.011 × 90 = 0.99 L/h` volgt exact `CL_ren_ugt = (21.56 − 5.44 − 0.99) / 0.9591 = 15.78 L/h`.

### 3.6 Drie UGT1A9-maturatiefuncties (`MAT`)

1. **`farhan2024` (standaard, Farhan et al. 2024 *in-vivo*):**
   $$\text{MAT}(\text{AGE}) = 0.05 + \frac{0.95 \cdot \text{AGE}^{1.9}}{0.09^{1.9} + \text{AGE}^{1.9}}$$
2. **`badee2019` (Badee et al. 2019 / Jo et al. 2021 scenario 2):**
   $$\text{MAT}(\text{AGE}) = \begin{cases} 0 & \text{AGE} \le 0 \\ \dfrac{1.305 \cdot \text{AGE}^{1.1}}{0.086^{1.1} + \text{AGE}^{1.1}} & 0 < \text{AGE} \le 1 \\ \max\!\left(1.0, \; 1.243 - 0.833 \exp(0.2(\text{AGE} - 25))\right) & 1 < \text{AGE} \le 20 \\ 1.0 & \text{AGE} > 20 \end{cases}$$
3. **`simcyp_s1` (Simcyp standaard / Jo et al. 2021 scenario 1):**
   $$\text{MAT}(\text{AGE}) = \begin{cases} 0.086 + \dfrac{(1.041 - 0.086) \cdot \text{AGE}^{1.2}}{1.16^{1.2} + \text{AGE}^{1.2}} & 0 \le \text{AGE} \le 15 \\ 1.0 & \text{AGE} > 15 \end{cases}$$

### 3.7 Typische individuele PK-parameters (`individual_params()`)

$$CL_{\text{tv}} = f_{u,\text{ratio}} \cdot \left[5.44 \cdot \text{MAT} \cdot F_{\text{size},CL} + 15.78 \cdot \text{MAT} \cdot F_{\text{size},CL} \cdot \text{RF} + 0.011 \cdot \text{eGFR} \cdot \left(\frac{\text{BSA}}{2.0}\right)\right]$$

$$V_{2,\text{tv}} = 73.47 \cdot F_{\text{size},V_2}, \qquad V_{3,\text{tv}} = 168.5 \cdot \left(\frac{\text{BWT}}{84.0}\right) \cdot f_{u,\text{ratio}}, \qquad Q_{\text{tv}} = 9.137 \cdot F_{\text{size},CL}$$

$$K_{a,\text{tv}} = 2.833 \cdot \left(0.45 + 0.55 \cdot \frac{\text{AGE}}{0.10 + \text{AGE}}\right), \qquad \text{ALAG}_1 = 0.5\text{ h}$$

---

## 4. Stochastische variabiliteit (`IIV` en `IOV`)

- **Inter-individuele variabiliteit (`IIV`, log-normaal; waarbij \(\mathrm{N}(\mu, \omega^2)\) een normale verdeling met gemiddelde \(\mu\) en variantie \(\omega^2\) aanduidt):**
  $$CL_i = CL_{\text{tv},i} \exp(\eta_{2,i}), \quad \eta_{2,i} \sim \mathrm{N}(0, 0.1136)$$
  $$V_{3,i} = V_{3,\text{tv},i} \exp(\eta_{3,i}), \quad \eta_{3,i} \sim \mathrm{N}(0, 0.8229)$$
  $$K_{a,i} = \min\!\left(50, \; \max\!\left(0.02, \; K_{a,\text{tv},i} \exp(\eta_{1,i})\right)\right), \quad \eta_{1,i} \sim \mathrm{N}(0, \omega^2_{Ka})$$
  met \(\omega^2_{Ka} = 0.36\) (`calibrated`, standaard), \(8.52\) (`penland_raw`), of \(0\) (`none`). Geen `IIV` op \(V_2\) en \(Q\).
- **Inter-occasion variabiliteit (`IOV` op relatieve biologische beschikbaarheid `F1`):**
  Voor opgegeven `iov_cv` (`%CV`) geldt \(\sigma_{\text{IOV}} = \sqrt{\ln\!\left(1 + (\text{IOV}_{\text{CV}}/100)^2\right)}\). Per toediening \(d\) wordt een onafhankelijke mean-gecentreerde multiplier getrokken zodat \(\mathrm{E}[F_{1,i,d}] = 1.0\):
  $$F_{1,i,d} = \exp(\kappa_{i,d}), \qquad \kappa_{i,d} \sim \mathrm{N}\!\left(-\tfrac{1}{2}\sigma_{\text{IOV}}^2, \; \sigma_{\text{IOV}}^2\right), \qquad D_{i,d} = 1000 \cdot \text{Dosis}_{\text{mg},i} \cdot F_{1,i,d}\text{ (µg)}$$
- **Residuele fout (`$SIGMA / EPS`)** wordt bewust **niet** gesimuleerd zodat de 90%-band het zuivere populatie-predictie-interval (`P5–P95`) over individuen weergeeft.

---

## 5. Exacte analytische `ADVAN4`-oplossing, steady state & metrieken

Voor elk individu zijn \(k_e = CL_i / V_{2,i}\), \(k_{12} = Q_i / V_{2,i}\), \(k_{21} = Q_i / V_{3,i}\) en \(S = k_e + k_{12} + k_{21}\):

$$\alpha, \beta = \frac{S \pm \sqrt{S^2 - 4 k_e k_{21}}}{2} \quad (\alpha > \beta > 0)$$

*(Numerieke singulariteitsguard: indien \(|K_{a,i} - \alpha| < 10^{-5}\) of \(|K_{a,i} - \beta| < 10^{-5}\), wordt \(K_{a,i}\) met \(+2\times 10^{-5}\) verschoven.)*

Definieer \(c_a = \frac{1}{(\alpha - K_a)(\beta - K_a)}\), \(c_1 = \frac{1}{(K_a - \alpha)(\beta - \alpha)}\), \(c_2 = \frac{1}{(K_a - \beta)(\alpha - \beta)}\). De exacte toestandspropagator \(\mathbf{A}(t) = \Phi(t)\mathbf{A}(0)\) over elk bolusvrij tijdsinterval \(t \ge 0\) heeft niet-nul elementen (\(E_a = e^{-K_a t}\), \(E_\alpha = e^{-\alpha t}\), \(E_\beta = e^{-\beta t}\)):

$$\Phi_{11}(t) = E_a$$
$$\Phi_{21}(t) = K_a\left[(k_{21} - K_a)c_a E_a + (k_{21} - \alpha)c_1 E_\alpha + (k_{21} - \beta)c_2 E_\beta\right]$$
$$\Phi_{31}(t) = K_a k_{12}\left(c_a E_a + c_1 E_\alpha + c_2 E_\beta\right)$$
$$\Phi_{22}(t) = \frac{(\alpha - k_{21})E_\alpha + (k_{21} - \beta)E_\beta}{\alpha - \beta}, \qquad \Phi_{23}(t) = \frac{k_{21}(E_\beta - E_\alpha)}{\alpha - \beta}$$
$$\Phi_{32}(t) = \frac{k_{12}(E_\beta - E_\alpha)}{\alpha - \beta}, \qquad \Phi_{33}(t) = \frac{(k_{21} - \beta)E_\alpha + (\alpha - k_{21})E_\beta}{\alpha - \beta}$$

- **Toestandsverloop binnen doseringsinterval \([0, \tau]\) met lagtijd \(t_{\text{lag}} = 0.5\text{ h}\) (\(\tau > t_{\text{lag}}\), indicatorfunctie \(\mathbf{1}(t > t_{\text{lag}})\)):**
  $$\mathbf{A}(t) = \Phi(t)\mathbf{A}(0) + \mathbf{1}(t > t_{\text{lag}}) D_d \begin{pmatrix} \Phi_{11}(t - t_{\text{lag}}) \\ \Phi_{21}(t - t_{\text{lag}}) \\ \Phi_{31}(t - t_{\text{lag}}) \end{pmatrix}$$
- **Directe analytische steady state (bij `steady_state = TRUE` en `IOV = 0%`):**
  Uit \(\mathbf{A}_{ss}(0) = \Phi(\tau)\mathbf{A}_{ss}(0) + D\boldsymbol{\phi}_1(\tau - t_{\text{lag}})\) volgt direct:
  $$A_{1,ss}(0) = \frac{D\Phi_{11}(\tau - t_{\text{lag}})}{1 - \Phi_{11}(\tau)}, \qquad \begin{pmatrix} 1 - \Phi_{22}(\tau) & -\Phi_{23}(\tau) \\ -\Phi_{32}(\tau) & 1 - \Phi_{33}(\tau) \end{pmatrix} \begin{pmatrix} A_{2,ss}(0) \\ A_{3,ss}(0) \end{pmatrix} = \begin{pmatrix} D\Phi_{21}(\tau - t_{\text{lag}}) + \Phi_{21}(\tau)A_{1,ss}(0) \\ D\Phi_{31}(\tau - t_{\text{lag}}) + \Phi_{31}(\tau)A_{1,ss}(0) \end{pmatrix}$$
- **Expositiemetrieken per individu:**
  - **Exacte massabalans-`AUCtau`:** \(\text{AUC}_\tau = \dfrac{(A_1(0) + A_2(0) + A_3(0) + D_{\text{last}}) - (A_1(\tau) + A_2(\tau) + A_3(\tau))}{CL}\)
  - **Ongebonden `AUCu,ss`:** \(\text{AUC}_u = f_u(\text{AGE}) \cdot \text{AUC}_\tau\)
  - **`Cavg` / `Cavgu`:** \(\text{AUC}_\tau / \tau\) en \(\text{AUC}_u / \tau\)
  - **`Cmax` en `Tmax`:** bepaald op het `0.1 h` uitvoerrooster en lokaal verfijnd op `0.01 h` resolutie in \([t_{\text{grid,max}} - 0.1, \; t_{\text{grid,max}} + 0.1]\)
  - **`Cmin`:** dalspiegel \(C(\tau) = A_2(\tau) / V_2\)
  - **Terminale halfwaardetijd:** \(t_{1/2,\beta} = \ln(2) / \beta\)

---

## 6. Zelfstandige R- én NONMEM-referentiecode voor directe reproductie (copy-paste)

Hieronder staan zowel de zelfstandige **R-implementatie** als de volledige **NONMEM 7.x control stream (`dapa_ped_hf_sim.ctl`)** en de **1-kind voorbeelddataset (`example_child_hf.csv`)** direct onder elkaar als kopieerbare tekstblokken (elk blok heeft rechtsboven een **Kopieer code**-knop).

### 6.1 Zelfstandig R-referentieblok

Onderstaande R-functie bevat alle formules uit §2–§7 zonder externe pakketafhankelijkheden en reproduceert exact dezelfde individuele PK/PD-parameters als `dapapksim` én NONMEM:

```r
dapapksim_reproduce_one <- function(age, bwt, height_m, sex_male = 1, egfr = NULL, creat_umol = NULL,
                                    dose_mg = 10, tau = 24, mat_option = "farhan2024",
                                    edema_pct = 0, rf_factor = 1.0, cl_hep_factor = 1.0,
                                    alb_gL = NULL, f1_rel = 1.0) {
  bwt_dry <- bwt / (1 + max(0, edema_pct) / 100)
  v_edema <- max(0, bwt - bwt_dry)
  hm      <- if (height_m > 3) height_m / 100 else height_m
  bmi     <- bwt_dry / max(hm^2, 0.09)
  bsa     <- sqrt((hm * 100 * bwt_dry) / 3600)

  ffm_ad <- if (sex_male == 1) 9270 * bwt_dry / (6680 + 216 * bmi) else 9270 * bwt_dry / (8780 + 244 * bmi)
  m_ffm  <- if (age >= 18) 1 else if (sex_male == 1) {
    if (age <= 0) 0.88 else 0.88 + 0.12 / (1 + (age / 13.4)^(-12.7))
  } else {
    if (age <= 0) 1.11 else 1.11 - 0.11 / (1 + (age / 7.1)^(-1.1))
  }
  ffm       <- min(ffm_ad * m_ffm, 0.97 * bwt_dry)
  f_size_cl <- ((ffm + 0.25 * (bwt_dry - ffm)) / 63.08)^0.75
  f_size_v2 <- (ffm + 0.20 * (bwt_dry - ffm)) / 61.68

  alb <- if (!is.null(alb_gL) && !is.na(alb_gL) && alb_gL > 0) alb_gL else {
    if (age >= 18) 44 else min(44, 1.1287 * log(max(1, age * 365.25)) + 33.746)
  }
  fu       <- 1 / (1 + (alb / 44) * ((1 - 0.086) / 0.086))
  fu_ratio <- fu / 0.086

  w_gfr <- max(bwt_dry, 0.4)
  ga    <- if (age < 0.1 && w_gfr < 2.5) min(40, max(24, 22 + 7.2 * w_gfr)) else 40
  pma   <- ga + 52.1775 * max(age, 0)
  ngfr  <- 115.0 * (pma^3.4 / (47.7^3.4 + pma^3.4)) * (0.9925^max(0, age - 35))
  if (is.null(egfr)) {
    if (!is.null(creat_umol) && !is.na(creat_umol) && creat_umol > 0 && age >= 1) {
      scr <- max(creat_umol / 88.4, 0.05)
      egfr <- if (age >= 18) {
        kap <- if (sex_male == 1) 0.9 else 0.7; alp <- if (sex_male == 1) -0.302 else -0.241
        142 * min(scr / kap, 1)^alp * max(scr / kap, 1)^(-1.2) * 0.9938^age * (if (sex_male == 1) 1 else 1.012)
      } else {
        k <- if (sex_male == 1) (if (age < 12) 39.0 * 1.008^(age - 12) else 39.0 * 1.045^(age - 12))
             else               (if (age < 12) 36.1 * 1.008^(age - 12) else 36.1 * 1.023^(age - 12))
        k * hm / scr
      }
      egfr <- min(max(egfr, 5), 180)
    } else egfr <- min(max(ngfr, 5), 180)
  }
  if (abs(rf_factor - 1.0) > 1e-12) egfr <- min(max(egfr * rf_factor, 5), 180)
  rf <- min(max(egfr / max(ngfr, 1.0), 0), 1.75)

  mat <- switch(mat_option,
    farhan2024 = 0.05 + 0.95 * age^1.9 / (0.09^1.9 + age^1.9),
    badee2019  = if (age <= 0) 0 else if (age <= 1) 1.305 * age^1.1 / (0.086^1.1 + age^1.1)
                 else if (age <= 20) max(1.0, 1.243 - 0.833 * exp(0.2 * (age - 25))) else 1.0,
    simcyp_s1  = if (age <= 15) 0.086 + (1.041 - 0.086) * age^1.2 / (1.16^1.2 + age^1.2) else 1.0
  )

  CL <- fu_ratio * (5.44 * mat * f_size_cl * cl_hep_factor +
                    15.78 * mat * f_size_cl * rf +
                    0.011 * egfr * (bsa / 2.0))
  V2 <- 73.47 * f_size_v2 + v_edema
  V3 <- 168.5 * (bwt_dry / 84.0) * fu_ratio
  Q  <- 9.137 * f_size_cl
  KA <- 2.833 * (0.45 + 0.55 * age / (0.10 + age))
  aucu_24h <- fu * ((dose_mg * 1000 * f1_rel) / CL) * (24 / tau)
  list(WDRY = bwt_dry, FFM = ffm, FU = fu, NGFR = ngfr, EGFR = egfr, RF = rf, MAT = mat,
       CL = CL, V2 = V2, V3 = V3, Q = Q, KA = KA,
       AUC_ss = (dose_mg * 1000 * f1_rel) / CL, AUCu_24h = aucu_24h)
}

# Voorbeeld — reproduceert exact hetzelfde hartfalen-voorbeeldkind (ID = 1) als in NONMEM hieronder:
dapapksim_reproduce_one(age = 0.75, bwt = 8.40, height_m = 0.71, sex_male = 1,
                        edema_pct = 5, rf_factor = 0.80, cl_hep_factor = 0.90,
                        alb_gL = 32, dose_mg = 1.20, tau = 24, mat_option = "farhan2024")
```

---

### 6.2 NONMEM 7.x control stream (`dapa_ped_hf_sim.ctl`)

Sla onderstaand blok op als `dapa_ped_hf_sim.ctl` in dezelfde map als `example_child_hf.csv` (§6.3) en voer uit met `nmfe75 dapa_ped_hf_sim.ctl dapa_ped_hf_sim.lst` (of via PsN: `execute dapa_ped_hf_sim.ctl`):

```fortran
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
```

---

### 6.3 Voorbeelddataset met 1 kind voor NONMEM (`example_child_hf.csv`)

Sla onderstaande CSV op als `example_child_hf.csv` naast `dapa_ped_hf_sim.ctl` (zuigeling van `9 maanden`, actueel gewicht `8,40 kg` met `5%` overvulling → drooggewicht `8,00 kg`, `RFFAC = 0,80`, `HEPFAC = 0,90`, `ALBIN = 32 g/L`, `1,20 mg` 1×/24u gedurende 7 dagen). Verwachte NONMEM-uitvoer op `ID = 1`: `CL = 4.564798 L/h`, `V2 = 8.458278 L`, `V3 = 21.37610 L`, `Q = 1.726980 L/h`, `KA = 2.649688 h⁻¹`, `FU = 0.114556`, `AUCU24 = 30.11452 ng·h/mL`, `CTOT(t=1h, Dag 1) = 83.623496 ng/mL`, `CTOT(t=145h, Dag 7) = 85.259729 ng/mL`:

```csv
ID,TIME,TAD,DAY,EVID,AMT,ADDL,II,CMT,DV,MDV,AGE,BWT,HT,SEX,CREAT,EGFRIN,ALBIN,EDEMA,RFFAC,HEPFAC,F1REL,MATOPT,DOSEMG,TAU
1,0,0,1,1,1.2,6,24,1,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,0,0,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,0.25,0.25,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,0.5,0.5,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,0.6,0.6,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,0.75,0.75,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,1,1,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,1.5,1.5,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,2,2,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,3,3,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,4,4,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,6,6,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,8,8,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,12,12,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,16,16,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,20,20,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,24,24,1,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,144,0,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,144.25,0.25,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,144.5,0.5,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,144.6,0.6,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,144.75,0.75,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,145,1,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,145.5,1.5,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,146,2,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,147,3,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,148,4,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,150,6,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,152,8,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,156,12,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,160,16,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,164,20,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
1,168,24,7,0,0,0,0,2,0,1,0.75,8.4,0.71,1,-99,-99,32,5,0.8,0.9,1,1,1.2,24
```

---

## 7. Uitbreiding voor pediatrisch hartfalen: pathofysiologie & duale PK–PD (`R/hf_pkpd.R`)

Voor dosisvoorspelling bij **pediatrisch hartfalen** (van pasgeborene tot adolescent) kan het standaard fysiologische model optioneel worden uitgebreid met hartfalen-specifieke pathofysiologie (standaard uitgeschakeld op neutrale referentiewaarden `1.0` / `0%`) en semi-mechanistische PD-uitkomsten:

1. **Extracellulaire overvulling vs. euvolemisch drooggewicht (`edema_pct`):**
   Bij percentage overvulling \(\text{FO}\) (`%` boven drooggewicht) geldt:
   $$W_{\text{droog}} = \frac{\text{BWT}}{1 + \text{FO}/100}, \qquad V_{\text{oedeem}} = \text{BWT} - W_{\text{droog}}\text{ (L)}$$
   Omdat oedeemvocht geen hepatocyten of niertubuli bevat, worden `FFM`, `NFM_CL`, `NFM_V2` en lipofiel \(V_3\) berekend uit \(W_{\text{droog}}\), terwijl het centrale/interstitiële verdelingsvolume wordt uitgebreid met \(+V_{\text{oedeem}}\).
2. **Cardiorenaal syndroom (`rf_factor`), hepatische stuwing (`cl_hep_factor`), hypoalbuminemie (`alb_override`) en sonde-opname (`f1_rel`):**
   - \(\text{eGFR}_{\text{eff}} = \min\!\left(180, \max(5, \; \text{eGFR} \cdot f_{\text{RF}})\right)\), waarmee zowel proximaal tubulaire UGT1A9-klaring (`15.78 · MAT · Fsize_CL · RF`) als passieve filtratie daalt bij verminderde renale perfusie (ook bij zuigelingen `< 1 j` waar `copulapop` geen creatinine gebruikt).
   - \(CL_{\text{hep}} = f_{u,\text{ratio}} \cdot 5.44 \cdot \text{MAT} \cdot F_{\text{size},CL} \cdot f_{\text{hep}}\) voor verminderde hepatische glucuronidering bij verhoogde centraal-veneuze druk (rechtszijdig falen / Fontan).
   - Handmatig ingesteld serumalbumine verhoogt bij hypoalbuminemie of *protein-losing enteropathy* de vrije fractie \(f_u\), waardoor totale \(C(t)\) en \(\text{AUC}_\tau\) dalen terwijl de farmacologisch actieve vrije concentratie \(C_u(t) = f_u C(t)\) en \(\text{AUC}_{u,24\text{h}}\) ongewijzigd blijven.
   - Gemiddelde relatieve biologische beschikbaarheid \(F_{1,\text{rel}}\) schaalt de toegediende orale/enterale dosis voor sonde- of suspensieverlies.
3. **Semi-mechanistische duale PK–PD relaties:**
   - **Proximale tubulaire luminale SGLT2-receptorbezetting:** met humane SGLT2 \(\text{IC}_{50} = 1.1\text{ nM} = 0.45\text{ ng/mL}\) en proximale S1/S2 waterresorptie-concentratiefactor \(\phi_{\text{PT}} = 2.25\) (Yakovleva et al. 2014) bedraagt de effectieve vrije plasmaconcentratie voor 50% luminale bezetting \(C_{u,50} = 0.20\text{ ng/mL}\):
     $$O_{\text{SGLT2}}(t) = \frac{100 \cdot C_u(t)}{0.20 + C_u(t)}\text{ (\%)}, \qquad \overline{O}_{\text{SGLT2}} = \frac{1}{\tau}\int_0^\tau O_{\text{SGLT2}}(t)\,dt$$
   - **GFR-gekoppelde effectieve renale osmotische/natriuretische PD-index:** corrigeert de SGLT2-blokkade voor de glomerulair gefiltreerde solute-aanvoer naar de proximale tubulus:
     $$\text{PD}_{\text{renaal,eff}} = \overline{O}_{\text{SGLT2}} \cdot \min\!\left(1.5, \; \frac{\text{eGFR}_{\text{eff}}}{100}\right)$$
   - **Systemische/cardiorenaal-pleiotrope \(E_{\max}\)-responsen (Kroonen et al. 2020):** op basis van gepubliceerde klinische blootstelling-responsmodellen omgerekend naar vrije 24-uurs steady-state blootstelling (\(\text{EAUC}_{u,50,\text{TGF}} = 5.2\text{ ng·h/mL}\), \(\text{EAUC}_{u,50,\text{HCT}} = 20.7\text{ ng·h/mL}\)).

---

## 8. Bronnen

1. Penland RC et al. *J Clin Pharmacol* 2023;63:551–559 — [doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196).
2. Al-Sallami HS et al. *Clin Pharmacokinet* 2015;54:1197–1208 — [doi:10.1007/s40262-015-0285-z](https://doi.org/10.1007/s40262-015-0285-z).
3. Janmahasatian S et al. *Clin Pharmacokinet* 2005;44:1051–1065 — [doi:10.2165/00003088-200544100-00004](https://doi.org/10.2165/00003088-200544100-00004).
4. Anderson BJ, Holford NHG. *Annu Rev Pharmacol Toxicol* 2008;48:303–332 (Normal Fat Mass).
5. McNamara PJ, Alcorn J. *AAPS PharmSci* 2002;4:E4 — [doi:10.1208/ps040104](https://doi.org/10.1208/ps040104).
6. Rhodin MM et al. *Pediatr Nephrol* 2009;24:67–76 — [doi:10.1007/s00467-008-0997-5](https://doi.org/10.1007/s00467-008-0997-5).
7. Pierce CB et al. *Kidney Int* 2021;99:948–956 — [doi:10.1016/j.kint.2020.10.047](https://doi.org/10.1016/j.kint.2020.10.047) (CKiD U25).
8. Inker LA et al. *N Engl J Med* 2021;385:1737–1749 (CKD-EPI 2021).
9. Farhan H et al. *J Clin Pharmacol* 2024;64:1224 — [doi:10.1002/jcph.2484](https://doi.org/10.1002/jcph.2484).
10. Jo H et al. *CPT Pharmacometrics Syst Pharmacol* 2021;10:108–118 — [doi:10.1002/psp4.12577](https://doi.org/10.1002/psp4.12577).
11. Tirucherai GS et al. *Diabetes Obes Metab* 2016;18:678–684 — [doi:10.1111/dom.12638](https://doi.org/10.1111/dom.12638).
12. Kasichayanula S et al. *Clin Pharmacokinet* 2014;53:17–27 — [doi:10.1007/s40262-013-0104-3](https://doi.org/10.1007/s40262-013-0104-3).
13. Zwep LHM et al. *Clin Pharmacol Ther* 2024;115:795–804 — [doi:10.1002/cpt.3099](https://doi.org/10.1002/cpt.3099) (`copulapop`).
14. Kroonen MYAM et al. *Br J Clin Pharmacol* 2020;86:1905–1916 — [doi:10.1111/bcp.14318](https://doi.org/10.1111/bcp.14318) (blootstelling-respons van dapagliflozine op cardiorenaal-markers).
15. Yakovleva T et al. *Front Pharmacol* 2014;5:218 — [doi:10.3389/fphar.2014.00218](https://doi.org/10.3389/fphar.2014.00218) (semi-mechanistische luminale SGLT2-modellering).

## 9. Beperkingen

- Extrapolatie onder ~1 jaar berust op mechanistische fysiologische ontogenie (UGT1A9, GFR,
  albumine, GI-absorptie) en PBPK-verificatie; gepubliceerde klinische PK-data van dapagliflozine
  bij neonaten en jonge zuigelingen zijn nog beperkt.
- Virtuele `copulapop`-populaties representeren algemene of ziekenhuispopulaties; de hartfalen-opties
  (`rf_factor`, `cl_hep_factor`, `edema_pct`, `alb_override`) bieden mechanistische scenario-analyse,
  maar vervangen geen prospectieve klinische PK/PD-validatie in pediatrisch hartfalen.
- Onderzoeks- en educatief simulatiemodel — geen gevalideerd klinisch doseringsadvies.
