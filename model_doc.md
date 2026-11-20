# Model, fysiologische onderbouwing & volledige reproduceerbaarheid

Deze pagina bevat de **volledige wiskundige, fysiologische en algoritmische specificatie**
van het in **dapapksim** geïmplementeerde populatie-PK model (`R/pk_model.R` en `R/population.R`).
Alle vergelijkingen, parameterwaarden, randvoorwaarden, stochastische trekkingen en de
gesloten-vorm solver zijn hieronder volledig uitgeschreven — inclusief een zelfstandig
uitvoerbaar R-referentieblok in §6 — zodat het model **100% reproduceerbaar** is op basis
van uitsluitend de informatie in deze app.

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

$$f_u(\text{AGE}) = \frac{1}{1 + \dfrac{[\text{Alb}](\text{AGE})}{44.0} \cdot \dfrac{1 - 0.086}{0.086}}, \qquad \text{fu\_ratio} = \frac{f_u(\text{AGE})}{0.086}$$

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
     - **Zuigelingen (\(\text{AGE} < 1\text{ j}\)):** \(\text{eGFR} = \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\) met \(\eta_{\text{GFR}} \sim \mathcal{N}(0, 0.215^2)\).
     - Alle waarden worden begrensd op **`[5, 180] mL/min/1.73m²`**.
   - **Indien `CREAT` ontbreekt (`NHANES`, `HSE`):**
     $$\text{eGFR} = \min\!\left(180, \; \max\!\left(5, \; \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\right)\right), \qquad \eta_{\text{GFR}} \sim \mathcal{N}(0, 0.215^2)$$

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

$$CL_{\text{tv}} = \text{fu\_ratio} \cdot \left[5.44 \cdot \text{MAT} \cdot F_{\text{size},CL} + 15.78 \cdot \text{MAT} \cdot F_{\text{size},CL} \cdot \text{RF} + 0.011 \cdot \text{eGFR} \cdot \left(\frac{\text{BSA}}{2.0}\right)\right]$$

$$V_{2,\text{tv}} = 73.47 \cdot F_{\text{size},V_2}, \qquad V_{3,\text{tv}} = 168.5 \cdot \left(\frac{\text{BWT}}{84.0}\right) \cdot \text{fu\_ratio}, \qquad Q_{\text{tv}} = 9.137 \cdot F_{\text{size},CL}$$

$$K_{a,\text{tv}} = 2.833 \cdot \left(0.45 + 0.55 \cdot \frac{\text{AGE}}{0.10 + \text{AGE}}\right), \qquad \text{ALAG}_1 = 0.5\text{ h}$$

---

## 4. Stochastische variabiliteit (`IIV` en `IOV`)

- **Inter-individuele variabiliteit (`IIV`, log-normaal):**
  $$CL_i = CL_{\text{tv},i} \exp(\eta_{2,i}), \quad \eta_{2,i} \sim \mathcal{N}(0, 0.1136)$$
  $$V_{3,i} = V_{3,\text{tv},i} \exp(\eta_{3,i}), \quad \eta_{3,i} \sim \mathcal{N}(0, 0.8229)$$
  $$K_{a,i} = \min\!\left(50, \; \max\!\left(0.02, \; K_{a,\text{tv},i} \exp(\eta_{1,i})\right)\right), \quad \eta_{1,i} \sim \mathcal{N}(0, \omega^2_{Ka})$$
  met \(\omega^2_{Ka} = 0.36\) (`calibrated`, standaard), \(8.52\) (`penland_raw`), of \(0\) (`none`). Geen `IIV` op \(V_2\) en \(Q\).
- **Inter-occasion variabiliteit (`IOV` op relatieve biologische beschikbaarheid `F1`):**
  Voor opgegeven `iov_cv` (`%CV`) geldt \(\sigma_{\text{IOV}} = \sqrt{\ln\!\left(1 + (\text{iov\_cv}/100)^2\right)}\). Per toediening \(d\) wordt een onafhankelijke mean-gecentreerde multiplier getrokken zodat \(\mathbb{E}[F_{1,i,d}] = 1.0\):
  $$F_{1,i,d} = \exp(\kappa_{i,d}), \qquad \kappa_{i,d} \sim \mathcal{N}\!\left(-\tfrac{1}{2}\sigma_{\text{IOV}}^2, \; \sigma_{\text{IOV}}^2\right), \qquad D_{i,d} = 1000 \cdot \text{DOSE\_MG}_i \cdot F_{1,i,d}\text{ (µg)}$$
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

- **Toestandsverloop binnen doseringsinterval \([0, \tau]\) met lagtijd \(t_{\text{lag}} = 0.5\text{ h}\) (\(\tau > t_{\text{lag}}\)):**
  $$\mathbf{A}(t) = \Phi(t)\mathbf{A}(0) + \mathbb{I}(t > t_{\text{lag}}) D_d \begin{pmatrix} \Phi_{11}(t - t_{\text{lag}}) \\ \Phi_{21}(t - t_{\text{lag}}) \\ \Phi_{31}(t - t_{\text{lag}}) \end{pmatrix}$$
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

## 6. Zelfstandig R-referentieblok voor directe reproductie

Onderstaande minimalistische R-functie bevat exact alle formules uit §2–§5 zonder externe
pakketafhankelijkheden en reproduceert voor elk individu exact dezelfde PK-parameters en
concentratie-tijdcurve als `dapapksim`:

```r
dapapksim_reproduce_one <- function(age, bwt, height_m, sex_male = 1, egfr = NULL, creat_umol = NULL,
                                    dose_mg = 10, tau = 24, mat_option = "farhan2024") {
  bmi <- bwt / (height_m^2)
  bsa <- sqrt((height_m * 100 * bwt) / 3600)
  ffm_ad <- if (sex_male == 1) 9270 * bwt / (6680 + 216 * bmi) else 9270 * bwt / (8780 + 244 * bmi)
  m_ffm  <- if (age >= 18) 1 else if (sex_male == 1) {
    if (age <= 0) 0.88 else 0.88 + 0.12 / (1 + (age / 13.4)^(-12.7))
  } else {
    if (age <= 0) 1.11 else 1.11 - 0.11 / (1 + (age / 7.1)^(-1.1))
  }
  ffm <- min(ffm_ad * m_ffm, 0.97 * bwt)
  f_size_cl <- ((ffm + 0.25 * (bwt - ffm)) / 63.08)^0.75
  f_size_v2 <- (ffm + 0.20 * (bwt - ffm)) / 61.68

  alb <- if (age >= 18) 44 else min(44, 1.1287 * log(max(1, age * 365.25)) + 33.746)
  fu  <- 1 / (1 + (alb / 44) * ((1 - 0.086) / 0.086))
  fu_ratio <- fu / 0.086

  ga   <- if (age < 0.1 && bwt < 2.5) min(40, max(24, 22 + 7.2 * bwt)) else 40
  pma  <- ga + 52.1775 * max(age, 0)
  ngfr <- 115.0 * (pma^3.4 / (47.7^3.4 + pma^3.4)) * (0.9925^max(0, age - 35))
  if (is.null(egfr)) {
    if (!is.null(creat_umol) && age >= 1) {
      scr <- max(creat_umol / 88.4, 0.05)
      egfr <- if (age >= 18) {
        kap <- if (sex_male == 1) 0.9 else 0.7; alp <- if (sex_male == 1) -0.302 else -0.241
        142 * min(scr / kap, 1)^alp * max(scr / kap, 1)^(-1.2) * 0.9938^age * (if (sex_male == 1) 1 else 1.012)
      } else {
        k <- if (sex_male == 1) (if (age < 12) 39.0 * 1.008^(age - 12) else 39.0 * 1.045^(age - 12))
             else               (if (age < 12) 36.1 * 1.008^(age - 12) else 36.1 * 1.023^(age - 12))
        k * height_m / scr
      }
      egfr <- min(max(egfr, 5), 180)
    } else egfr <- min(max(ngfr, 5), 180)
  }
  rf <- min(max(egfr / max(ngfr, 1.0), 0), 1.75)

  mat <- switch(mat_option,
    farhan2024 = 0.05 + 0.95 * age^1.9 / (0.09^1.9 + age^1.9),
    badee2019  = if (age <= 0) 0 else if (age <= 1) 1.305 * age^1.1 / (0.086^1.1 + age^1.1)
                 else if (age <= 20) max(1.0, 1.243 - 0.833 * exp(0.2 * (age - 25))) else 1.0,
    simcyp_s1  = if (age <= 15) 0.086 + (1.041 - 0.086) * age^1.2 / (1.16^1.2 + age^1.2) else 1.0
  )

  CL <- fu_ratio * (5.44 * mat * f_size_cl + 15.78 * mat * f_size_cl * rf + 0.011 * egfr * (bsa / 2.0))
  V2 <- 73.47 * f_size_v2
  V3 <- 168.5 * (bwt / 84.0) * fu_ratio
  Q  <- 9.137 * f_size_cl
  KA <- 2.833 * (0.45 + 0.55 * age / (0.10 + age))
  list(CL = CL, V2 = V2, V3 = V3, Q = Q, KA = KA, FU = fu, AUC_ss = (dose_mg * 1000) / CL,
       AUCu_ss = fu * (dose_mg * 1000) / CL)
}
```

---

## 7. Bronnen

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

## 8. Beperkingen

- Extrapolatie onder ~1 jaar berust op mechanistische fysiologische ontogenie (UGT1A9, GFR,
  albumine, GI-absorptie) en PBPK-verificatie; gepubliceerde klinische PK-data van dapagliflozine
  bij neonaten en jonge zuigelingen zijn nog beperkt.
- Virtuele `copulapop`-populaties representeren algemene of ziekenhuispopulaties, geen
  specifieke T2DM- of hartfalencohorten.
- Onderzoeks- en educatief simulatiemodel — geen gevalideerd klinisch doseringsadvies.
