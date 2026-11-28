# Complete Mathematical, Algorithmic & NONMEM-Reproducible Model Specification

This sub-page contains the **exact mathematical equations, parameter values, boundary conditions, stochastic draws, and closed-form `ADVAN4` solver** implemented in **dapapksim** (`R/pk_model.R`, `R/population.R`, and `R/hf_pkpd.R`). See sub-tab **1 (*Step-by-Step Physiological Derivation & Calibration*)** for the clinical-pharmacological rationale and calibration background. Alongside the standalone R reference block in §6, the **`nonmem/`** directory in the repository provides a ready-to-run NONMEM 7.x control stream (`nonmem/dapa_ped_hf_sim.ctl`) and a 1-child example dataset (`nonmem/example_child_hf.csv`) for direct reproduction in NONMEM.

---

## 1. Summary of Mechanistic Refinements vs. Penland et al. 2023

The starting point is the two-compartment oral NONMEM popPK model of **Penland et al. 2023**
(*J Clin Pharmacol* 2023;63:551–559, [doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196)),
developed on adult T2DM, DAPA-CKD, and adolescent T2DM data (`11–17 years`). Because this
application is designed for **free extrapolation across any age and dose** (including
neonates, infants, severe obesity, chronic kidney disease, and paediatric heart failure), six
empirical limitations of the source model were replaced by mechanistic physiological scaling:

| Component | Empirical source / initial model | Current physiological model (`dapapksim`) | Physiological & clinical rationale |
|---|---|---|---|
| **1. Body size & obesity** | Total body mass `(BWT/84)^0.75` and `BWT/84` | **Fat-Free Mass (`FFM`) & Normal Fat Mass (`NFM`)** for `CL`, `Q` (`Ffat = 0.25`) and `V2` (`Ffat = 0.20`); total `BWT` for lipophilic `V3` | Adipose tissue does not contribute proportionally to hepatic/renal mass. In obese adolescents (Tirucherai 2016, mean `99.7 kg`), `BWT/84` overpredicted clearance and central volume; `FFM`/`NFM` also accounts for sex differences so `THETA(7)` is omitted. |
| **2. Protein binding (`fu`)** | Implicitly constant (`fu` independent of age) | **Albumin ontogeny (McNamara & Alcorn 2002)** anchored at `fu_adult = 0.086` (`91.4%` bound at `44 g/L`); explicit **`AUCu,ss`** reporting | Neonates have lower serum albumin (`~34 g/L`), raising `fu` by ~27% (`0.109`). For a low-extraction drug (`ER < 0.05`), total `CL/F` and `V3/F` scale linearly with `fu/fu_adult`, while pharmacologically active unbound exposure (`AUCu,ss`) is governed by intrinsic clearance. |
| **3. Renal function (`eGFR` & `nGFR`)** | Fixed bedside Schwartz (`0.413`) across `0–18 y`; missing neonatal maturation | **`CKD-EPI 2021`** (`≥ 18 y`), **`Pierce 2021 CKiD U25`** (`1–<18 y` when `CREAT` present), and **`Rhodin 2009` PMA maturation** (`< 1 y` and datasets without `CREAT`), bounded to `[5, 180]` | Cross-sectional `0–18 y` copula creatinine correlations miss maternal creatinine transfer in newborns and yielded unphysiologically high neonatal Schwartz values (`> 110`). `CKiD U25` removes age/sex bias across `1–18 y`. |
| **4. Renal clearance & renal UGT1A9** | Entire `eGFR`-dependent clearance (`~74%` of `CL`) treated as passive filtration without enzyme maturation | **Three-pathway decomposition:** hepatic UGT1A9 (`5.44 L/h`) + **renal tubular UGT1A9** (`15.78 L/h`, scaled by both `MAT` and `RF = eGFR/nGFR`) + unchanged filtration (`0.011 × eGFR × BSA/2`) | Less than `2%` of dapagliflozin is excreted unchanged in urine; the steep `eGFR` effect reflects proximal tubular UGT1A9 glucuronidation to D3OG (Kasichayanula 2014, Jo 2021). Using `RF = eGFR/nGFR` lets renal UGT1A9 mature via `MAT` without double-counting low neonatal GFR while declining proportionally in CKD (`RF < 1`). |
| **5. UGT1A9 ontogeny** | Two options with adult normalisation discontinuity | **Three validated options (`MAT(adult) = 1.0`):** `Farhan 2024` (default), `Badee 2019` (Jo s2), and `Simcyp s1` (Jo s1) | Covers the full literature spectrum from rapid *in-vivo* maturation to conservative slow *in-vitro* ontogeny without discontinuity at `18–20 years`. |
| **6. GI absorption, `Ka`-BSV & solver** | `(BWT/84)^-0.25` (faster `KA` in babies!), `ω²_Ka = 8.52 FIX`, and fixed-step numerical RK4 | **Neonatal gastric-emptying maturation `f_GI(AGE)`**, calibrated **`ω²_Ka = 0.36`** (default), mean-centred `IOV`, and **exact analytical `ADVAN4` solver** | Gastric emptying is slower in newborns, not faster. `ω²_Ka = 8.52` was a sparse-trough artefact in Penland (`Cmax` bias `-51%`, CV `126%`); `ω²_Ka = 0.36` reproduces clinical `Cmax` and `35%` total CV. The analytical solver is exact to machine precision and ~50× faster. |

---

## 2. Compartmental Structure (`ADVAN4 TRANS4`)

The structural model comprises three compartments with amounts \((A_1, A_2, A_3)\) in **micrograms (`µg`)**, time \(t\) in **hours (`h`)**, volumes in **litres (`L`)**, and clearances in **`L/h`**:

1. **Absorption depot (\(A_1\), gut):** oral dose \(D\) (in `µg` \(= 1000 \times \text{dose}_{\text{mg}}\)) enters the depot after absorption lag time \(\text{ALAG}_1 = 0.5\text{ h}\) and is absorbed via first-order rate constant \(K_a\) (`h⁻¹`).
2. **Central compartment (\(A_2\), volume \(V_2\)):** receives \(K_a A_1\), eliminates linearly with apparent total clearance \(CL\) (\(k_e = CL / V_2\)), and exchanges with the peripheral compartment via \(Q\) (\(k_{12} = Q / V_2\), \(k_{21} = Q / V_3\)).
3. **Peripheral compartment (\(A_3\), volume \(V_3\)):** tissue distribution compartment.

Plasma concentration is \(C(t) = A_2(t) / V_2\) in **`µg/L = ng/mL`**.

$$\frac{dA_1}{dt} = -K_a A_1, \qquad \frac{dA_2}{dt} = K_a A_1 - (k_e + k_{12})A_2 + k_{21}A_3, \qquad \frac{dA_3}{dt} = k_{12}A_2 - k_{21}A_3$$

---

## 3. Complete Parameter Specification & Covariate Equations

### 3.1 Fixed population parameters (`nonmem_theta()`)

| Code symbol | Value | Unit | Description & provenance |
|---|---|---|---|
| `ka` | `2.833` | `h⁻¹` | Adult absorption rate constant `THETA(1)` (Penland 2023) |
| `cl_hep` | `5.44` | `L/h` | Hepatic UGT1A9 clearance (`CL/F`) at reference size, `MAT = 1`, `fu = 0.086` (derived from Penland Table S4 intercept as `eGFR → 0`) |
| `cl_ren_ugt` | `15.78` | `L/h` | Renal proximal tubular UGT1A9 clearance (`CL/F`) at reference size, `MAT = 1`, `RF = 1`, `fu = 0.086` (derivation in §3.5) |
| `cl_filt_s` | `0.011` | `L/h` per `mL/min/1.73m²` | Unchanged glomerular/renal excretion slope at `BSA_ref = 2.0 m²` (`0.99 L/h` at `eGFR = 90`, ~4.6% of total `CL/F`; Kasichayanula 2014) |
| `v2` | `73.47` | `L` | Central volume of distribution `Vc/F` (`THETA(3)`) at `NFM_V2,ref = 61.68 kg` |
| `v3` | `168.5` | `L` | Peripheral volume of distribution `Vp/F` (`THETA(4)`) at `BWT_ref = 84.0 kg`, `fu = 0.086` |
| `q` | `9.137` | `L/h` | Intercompartmental clearance `Q/F` (`THETA(5)`) at `NFM_CL,ref = 63.08 kg` |
| `alag` | `0.5` | `h` | Absorption lag time `ALAG1` (fixed) |
| `bwt_ref` / `bsa_ref` | `84.0` / `2.0` | `kg` / `m²` | Reference weight and BSA of the Penland source population (`84 kg`, `1.71 m`) |
| `ffm_ref` | `56.1` | `kg` | Sex-weighted reference FFM of the Penland population (`~65%` male, `84 kg`, `BMI 28.7`) |
| `nfm_cl_ref` / `nfm_v2_ref` | `63.08` / `61.68` | `kg` | Reference NFM: `56.1 + 0.25×(84−56.1) = 63.08` and `56.1 + 0.20×(84−56.1) = 61.68` |
| `fu_adult` / `alb_adult` | `0.086` / `44.0` | `–` / `g/L` | Adult unbound fraction (`91.4%` bound) and reference serum albumin |
| `gi_f0` / `gi_t50` | `0.45` / `0.10` | `–` / `years` | Neonatal fraction of adult gastric emptying at birth and maturation half-time (`~5 weeks`) |
| `tm50` / `hill` / `fbirth` | `0.09` / `1.9` / `0.05` | `years` / `–` / `–` | Farhan et al. 2024 *in-vivo* UGT1A9 maturation parameters |

### 3.2 Covariates from the virtual population (`copulapop`)

For each virtual individual:
- `BWT = WEIGHT` (`kg`), `AGE` (`years`, bounded \(\ge 0\)), `SEXISMALE` (`1` = male/boy, `0` = female/girl).
- `HEIGHT` in `m` (if an input table supplies `HEIGHT > 3`, values are automatically divided by `100` to convert `cm → m`).
- `BMI = BWT / max(HEIGHT^2, 0.09)` (\(\text{kg/m}^2\)).
- `BSA = sqrt((HEIGHT * 100 * BWT) / 3600)` (\(\text{m}^2\), Mosteller), unless `BSA` is explicitly supplied.
- *Fallback for synthetic validation tables lacking `HEIGHT`:* if `BSA` is known, `HEIGHT` is inverted from Mosteller (`min(2.15, max(0.35, (BSA^2 * 3600 / max(BWT, 0.5)) / 100))`); otherwise a smooth growth curve from `0.50 m` to `1.76 m` (male) / `1.63 m` (female) is used (`estimate_height_m()`).

### 3.3 Body composition (`FFM` and `NFM`)

Following **Janmahasatian et al. 2005** (\(\ge 18\text{ y}\)) and **Al-Sallami et al. 2015** (\(< 18\text{ y}\)):

$$\text{FFM}_{\text{adult}} = \begin{cases} \dfrac{9270 \cdot \text{BWT}}{6680 + 216 \cdot \text{BMI}} & \text{if SEXISMALE} = 1 \\[10pt] \dfrac{9270 \cdot \text{BWT}}{8780 + 244 \cdot \text{BMI}} & \text{if SEXISMALE} = 0 \end{cases}$$

$$M_{\text{FFM}}(\text{AGE}) = \begin{cases} 1.0 & \text{if AGE} \ge 18 \\[6pt] 0.88 + \dfrac{0.12}{1 + (\text{AGE}/13.4)^{-12.7}} & \text{if AGE} < 18 \text{ and male (}0.88\text{ at AGE}=0\text{)} \\[10pt] 1.11 - \dfrac{0.11}{1 + (\text{AGE}/7.1)^{-1.1}} & \text{if AGE} < 18 \text{ and female (}1.11\text{ at AGE}=0\text{)} \end{cases}$$

$$\text{FFM} = \min\!\left(\text{FFM}_{\text{adult}} \cdot M_{\text{FFM}}(\text{AGE}), \; 0.97 \cdot \text{BWT}\right)$$

$$\text{NFM}_{CL} = \text{FFM} + 0.25(\text{BWT} - \text{FFM}), \qquad F_{\text{size},CL} = \left(\frac{\text{NFM}_{CL}}{63.08}\right)^{0.75}$$

$$\text{NFM}_{V_2} = \text{FFM} + 0.20(\text{BWT} - \text{FFM}), \qquad F_{\text{size},V_2} = \frac{\text{NFM}_{V_2}}{61.68}$$

### 3.4 Ontogeny of plasma protein binding (`fu_mcnamara()`)

Following **McNamara & Alcorn 2002** (with \(\text{days} = \max(1.0, \; 365.25 \times \text{AGE})\) when `ALBUMIN` is unmeasured):

$$[\text{Alb}](\text{AGE}) = \begin{cases} 44.0\text{ g/L} & \text{if AGE} \ge 18 \\ \min\!\left(44.0, \; 1.1287 \ln(\text{days}) + 33.746\right)\text{ g/L} & \text{if AGE} < 18 \end{cases}$$

$$f_u(\text{AGE}) = \frac{1}{1 + \dfrac{[\text{Alb}](\text{AGE})}{44.0} \cdot \dfrac{1 - 0.086}{0.086}}, \qquad f_{u,\text{ratio}} = \frac{f_u(\text{AGE})}{0.086}$$

### 3.5 Renal function ontogeny (`ngfr_rhodin()`) and individual `eGFR` (`egfr_from_copula()`)

1. **Healthy age-expected GFR by postmenstrual age (`nGFR` in `mL/min/1.73m²`; Rhodin et al. 2009 + senescence above 35 y):**
   - Estimated gestational age at birth (`GA` in weeks) for small neonates:
     $$\text{GA} = \begin{cases} \min\!\left(40, \max(24, \; 22 + 7.2 \cdot \text{BWT})\right) & \text{if AGE} < 0.1\text{ y and BWT} < 2.5\text{ kg} \\ 40 & \text{otherwise} \end{cases}$$
   - Postmenstrual age \(\text{PMA} = \text{GA} + 52.1775 \cdot \text{AGE}\) (weeks).
   - $$\text{nGFR}(\text{AGE}, \text{BWT}) = 115.0 \cdot \frac{\text{PMA}^{3.4}}{47.7^{3.4} + \text{PMA}^{3.4}} \cdot 0.9925^{\max(0, \; \text{AGE} - 35)}$$

2. **Individual `eGFR` (`mL/min/1.73m²`):**
   - **When `CREAT` (µmol/L) is present in the dataset (`RUMC`):** let \(\text{SCr} = \max(\text{CREAT} / 88.4, \; 0.05)\) in `mg/dL`:
     - **Adults (\(\text{AGE} \ge 18\text{ y}\), CKD-EPI 2021 race-free):**
       with \(\kappa = 0.9\) (male) / \(0.7\) (female), \(\alpha = -0.302\) (male) / \(-0.241\) (female), \(f_{\text{sex}} = 1.0\) (male) / \(1.012\) (female):
       $$\text{eGFR} = 142 \cdot \min(\text{SCr}/\kappa, 1)^{\alpha} \cdot \max(\text{SCr}/\kappa, 1)^{-1.200} \cdot 0.9938^{\text{AGE}} \cdot f_{\text{sex}}$$
     - **Children and adolescents (\(1 \le \text{AGE} < 18\text{ y}\), Pierce et al. 2021 CKiD U25):**
       $$k = \begin{cases} 39.0 \cdot 1.008^{\text{AGE}-12} & \text{male, AGE} < 12 \\ 39.0 \cdot 1.045^{\text{AGE}-12} & \text{male, } 12 \le \text{AGE} < 18 \\ 36.1 \cdot 1.008^{\text{AGE}-12} & \text{female, AGE} < 12 \\ 36.1 \cdot 1.023^{\text{AGE}-12} & \text{female, } 12 \le \text{AGE} < 18 \end{cases}, \qquad \text{eGFR} = \frac{k \cdot \text{HEIGHT}_m}{\text{SCr}}$$
     - **Infants (\(\text{AGE} < 1\text{ y}\)):** \(\text{eGFR} = \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\) with \(\eta_{\text{GFR}} \sim \mathrm{N}(0, 0.215^2)\).
     - All values are bounded to **`[5, 180] mL/min/1.73m²`**.
   - **When `CREAT` is absent (`NHANES`, `HSE`):**
     $$\text{eGFR} = \min\!\left(180, \; \max\!\left(5, \; \text{nGFR}(\text{AGE}, \text{BWT}) \cdot \exp(\eta_{\text{GFR}})\right)\right), \qquad \eta_{\text{GFR}} \sim \mathrm{N}(0, 0.215^2)$$

3. **Relative renal function (`RF`) and 3-pathway clearance anchor:**
   $$\text{RF} = \min\!\left(1.75, \; \max\!\left(0, \; \frac{\text{eGFR}}{\max(\text{nGFR}(\text{AGE}, \text{BWT}), 1.0)}\right)\right)$$
   *Anchor derivation:* In Penland et al. 2023, population `CL/F = 21.56 L/h` at `BWT = 84 kg`, `BSA = 2.0 m²`, `eGFR = 90`, and median age `62 years` (`nGFR(62) = 93.84`, so `RF = 90 / 93.84 = 0.9591`). With `CL_hep = 5.44 L/h` and unchanged filtration `0.011 × 90 = 0.99 L/h`, exact mass balance gives `CL_ren_ugt = (21.56 − 5.44 − 0.99) / 0.9591 = 15.78 L/h`.

### 3.6 Three UGT1A9 maturation functions (`MAT`)

1. **`farhan2024` (default, Farhan et al. 2024 *in-vivo*):**
   $$\text{MAT}(\text{AGE}) = 0.05 + \frac{0.95 \cdot \text{AGE}^{1.9}}{0.09^{1.9} + \text{AGE}^{1.9}}$$
2. **`badee2019` (Badee et al. 2019 / Jo et al. 2021 scenario 2):**
   $$\text{MAT}(\text{AGE}) = \begin{cases} 0 & \text{AGE} \le 0 \\ \dfrac{1.305 \cdot \text{AGE}^{1.1}}{0.086^{1.1} + \text{AGE}^{1.1}} & 0 < \text{AGE} \le 1 \\ \max\!\left(1.0, \; 1.243 - 0.833 \exp(0.2(\text{AGE} - 25))\right) & 1 < \text{AGE} \le 20 \\ 1.0 & \text{AGE} > 20 \end{cases}$$
3. **`simcyp_s1` (Simcyp default / Jo et al. 2021 scenario 1):**
   $$\text{MAT}(\text{AGE}) = \begin{cases} 0.086 + \dfrac{(1.041 - 0.086) \cdot \text{AGE}^{1.2}}{1.16^{1.2} + \text{AGE}^{1.2}} & 0 \le \text{AGE} \le 15 \\ 1.0 & \text{AGE} > 15 \end{cases}$$

### 3.7 Typical individual PK parameters (`individual_params()`)

$$CL_{\text{tv}} = f_{u,\text{ratio}} \cdot \left[5.44 \cdot \text{MAT} \cdot F_{\text{size},CL} + 15.78 \cdot \text{MAT} \cdot F_{\text{size},CL} \cdot \text{RF} + 0.011 \cdot \text{eGFR} \cdot \left(\frac{\text{BSA}}{2.0}\right)\right]$$

$$V_{2,\text{tv}} = 73.47 \cdot F_{\text{size},V_2}, \qquad V_{3,\text{tv}} = 168.5 \cdot \left(\frac{\text{BWT}}{84.0}\right) \cdot f_{u,\text{ratio}}, \qquad Q_{\text{tv}} = 9.137 \cdot F_{\text{size},CL}$$

$$K_{a,\text{tv}} = 2.833 \cdot \left(0.45 + 0.55 \cdot \frac{\text{AGE}}{0.10 + \text{AGE}}\right), \qquad \text{ALAG}_1 = 0.5\text{ h}$$

---

## 4. Stochastic Variability (`IIV` and `IOV`)

- **Inter-individual variability (`IIV`, log-normal; where \(\mathrm{N}(\mu, \omega^2)\) denotes a normal distribution with mean \(\mu\) and variance \(\omega^2\)):**
  $$CL_i = CL_{\text{tv},i} \exp(\eta_{2,i}), \quad \eta_{2,i} \sim \mathrm{N}(0, 0.1136)$$
  $$V_{3,i} = V_{3,\text{tv},i} \exp(\eta_{3,i}), \quad \eta_{3,i} \sim \mathrm{N}(0, 0.8229)$$
  $$K_{a,i} = \min\!\left(50, \; \max\!\left(0.02, \; K_{a,\text{tv},i} \exp(\eta_{1,i})\right)\right), \quad \eta_{1,i} \sim \mathrm{N}(0, \omega^2_{Ka})$$
  with \(\omega^2_{Ka} = 0.36\) (`calibrated`, default), \(8.52\) (`penland_raw`), or \(0\) (`none`). No `IIV` on \(V_2\) or \(Q\).
- **Inter-occasion variability (`IOV` on relative bioavailability `F1`):**
  For specified `iov_cv` (`%CV`), \(\sigma_{\text{IOV}} = \sqrt{\ln\!\left(1 + (\text{IOV}_{\text{CV}}/100)^2\right)}\). At each dosing occasion \(d\), an independent mean-centred multiplier is drawn such that \(\mathrm{E}[F_{1,i,d}] = 1.0\):
  $$F_{1,i,d} = \exp(\kappa_{i,d}), \qquad \kappa_{i,d} \sim \mathrm{N}\!\left(-\tfrac{1}{2}\sigma_{\text{IOV}}^2, \; \sigma_{\text{IOV}}^2\right), \qquad D_{i,d} = 1000 \cdot \text{Dose}_{\text{mg},i} \cdot F_{1,i,d}\text{ (µg)}$$
- **Residual error (`$SIGMA / EPS`)** is intentionally **omitted** so that the 90% band represents the true inter-individual prediction interval (`P5–P95`).

---

## 5. Exact Analytical `ADVAN4` Solution, Steady State & Exposure Metrics

For each individual, let \(k_e = CL_i / V_{2,i}\), \(k_{12} = Q_i / V_{2,i}\), \(k_{21} = Q_i / V_{3,i}\), and \(S = k_e + k_{12} + k_{21}\):

$$\alpha, \beta = \frac{S \pm \sqrt{S^2 - 4 k_e k_{21}}}{2} \quad (\alpha > \beta > 0)$$

*(Numerical singularity guard: if \(|K_{a,i} - \alpha| < 10^{-5}\) or \(|K_{a,i} - \beta| < 10^{-5}\), \(K_{a,i}\) is shifted by \(+2\times 10^{-5}\)).*

Define \(c_a = \frac{1}{(\alpha - K_a)(\beta - K_a)}\), \(c_1 = \frac{1}{(K_a - \alpha)(\beta - \alpha)}\), \(c_2 = \frac{1}{(K_a - \beta)(\alpha - \beta)}\). The exact state-transition propagator \(\mathbf{A}(t) = \Phi(t)\mathbf{A}(0)\) over any bolus-free interval \(t \ge 0\) has non-zero entries (\(E_a = e^{-K_a t}\), \(E_\alpha = e^{-\alpha t}\), \(E_\beta = e^{-\beta t}\)):

$$\Phi_{11}(t) = E_a$$
$$\Phi_{21}(t) = K_a\left[(k_{21} - K_a)c_a E_a + (k_{21} - \alpha)c_1 E_\alpha + (k_{21} - \beta)c_2 E_\beta\right]$$
$$\Phi_{31}(t) = K_a k_{12}\left(c_a E_a + c_1 E_\alpha + c_2 E_\beta\right)$$
$$\Phi_{22}(t) = \frac{(\alpha - k_{21})E_\alpha + (k_{21} - \beta)E_\beta}{\alpha - \beta}, \qquad \Phi_{23}(t) = \frac{k_{21}(E_\beta - E_\alpha)}{\alpha - \beta}$$
$$\Phi_{32}(t) = \frac{k_{12}(E_\beta - E_\alpha)}{\alpha - \beta}, \qquad \Phi_{33}(t) = \frac{(k_{21} - \beta)E_\alpha + (\alpha - k_{21})E_\beta}{\alpha - \beta}$$

- **State trajectory within dosing interval \([0, \tau]\) with lag time \(t_{\text{lag}} = 0.5\text{ h}\) (\(\tau > t_{\text{lag}}\), indicator \(\mathbf{1}(t > t_{\text{lag}})\)):**
  $$\mathbf{A}(t) = \Phi(t)\mathbf{A}(0) + \mathbf{1}(t > t_{\text{lag}}) D_d \begin{pmatrix} \Phi_{11}(t - t_{\text{lag}}) \\ \Phi_{21}(t - t_{\text{lag}}) \\ \Phi_{31}(t - t_{\text{lag}}) \end{pmatrix}$$
- **Direct analytical steady state (when `steady_state = TRUE` and `IOV = 0%`):**
  Solving \(\mathbf{A}_{ss}(0) = \Phi(\tau)\mathbf{A}_{ss}(0) + D\boldsymbol{\phi}_1(\tau - t_{\text{lag}})\) in closed form gives:
  $$A_{1,ss}(0) = \frac{D\Phi_{11}(\tau - t_{\text{lag}})}{1 - \Phi_{11}(\tau)}, \qquad \begin{pmatrix} 1 - \Phi_{22}(\tau) & -\Phi_{23}(\tau) \\ -\Phi_{32}(\tau) & 1 - \Phi_{33}(\tau) \end{pmatrix} \begin{pmatrix} A_{2,ss}(0) \\ A_{3,ss}(0) \end{pmatrix} = \begin{pmatrix} D\Phi_{21}(\tau - t_{\text{lag}}) + \Phi_{21}(\tau)A_{1,ss}(0) \\ D\Phi_{31}(\tau - t_{\text{lag}}) + \Phi_{31}(\tau)A_{1,ss}(0) \end{pmatrix}$$
- **Individual exposure metrics:**
  - **Exact mass-balance `AUCtau`:** \(\text{AUC}_\tau = \dfrac{(A_1(0) + A_2(0) + A_3(0) + D_{\text{last}}) - (A_1(\tau) + A_2(\tau) + A_3(\tau))}{CL}\)
  - **Unbound `AUCu,ss`:** \(\text{AUC}_u = f_u(\text{AGE}) \cdot \text{AUC}_\tau\)
  - **`Cavg` / `Cavgu`:** \(\text{AUC}_\tau / \tau\) and \(\text{AUC}_u / \tau\)
  - **`Cmax` and `Tmax`:** evaluated on the `0.1 h` output grid and locally refined at `0.01 h` resolution over \([t_{\text{grid,max}} - 0.1, \; t_{\text{grid,max}} + 0.1]\)
  - **`Cmin`:** trough concentration \(C(\tau) = A_2(\tau) / V_2\)
  - **Terminal half-life:** \(t_{1/2,\beta} = \ln(2) / \beta\)

---

## 6. Standalone R & NONMEM Reference Code for Direct Reproduction (Copy-Paste)

Below are the standalone **R implementation**, the complete **NONMEM 7.x control stream (`dapa_ped_hf_sim.ctl`)**, and the **1-child example dataset (`example_child_hf.csv`)** as directly copyable code blocks (each block includes a **Copy code** button in its top-right corner).

### 6.1 Standalone R reference block

The R function below implements all equations from §2–§7 without external package dependencies and reproduces the exact same individual PK/PD parameters as `dapapksim` and NONMEM:

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

# Example — reproduces the exact same paediatric HF example child (ID = 1) as in NONMEM below:
dapapksim_reproduce_one(age = 0.75, bwt = 8.40, height_m = 0.71, sex_male = 1,
                        edema_pct = 5, rf_factor = 0.80, cl_hep_factor = 0.90,
                        alb_gL = 32, dose_mg = 1.20, tau = 24, mat_option = "farhan2024")
```

---

### 6.2 NONMEM 7.x control stream (`dapa_ped_hf_sim.ctl`)

Save the block below as `dapa_ped_hf_sim.ctl` in the same folder as `example_child_hf.csv` (§6.3) and run with `nmfe75 dapa_ped_hf_sim.ctl dapa_ped_hf_sim.lst` (or via PsN: `execute dapa_ped_hf_sim.ctl`):

```fortran
$SIZES PD=50
$PROBLEM Dapagliflozin 6-pillar physiological popPK/PD model + paediatric heart failure (dapapksim)
; =============================================================================
; Complete NONMEM implementation of the dapapksim model (R/pk_model.R & R/hf_pkpd.R)
; Basis: Penland et al. J Clin Pharmacol 2023;63:551-559 (doi:10.1002/jcph.2196),
; extended with:
;   1. FFM/NFM allometry on euvolaemic dry weight (Al-Sallami 2015; Janmahasatian 2005)
;   2. Plasma albumin ontogeny & unbound fraction fu (McNamara & Alcorn 2002; fu_adult = 0.086)
;   3. Renal function: Rhodin 2009 PMA-nGFR (<1 y / no CREAT), Pierce 2021 CKiD U25 (1-<18 y),
;      and CKD-EPI 2021 (>=18 y), bounded to [5, 180] mL/min/1.73m2
;   4. Three-pathway clearance: hepatic UGT1A9 + renal tubular UGT1A9 + passive filtration
;   5. UGT1A9 ontogeny switch via MATOPT: 1 = Farhan 2024 (default), 2 = Badee 2019, 3 = Simcyp s1
;   6. Developmental GI absorption on KA + heart failure physiology (RFFAC, HEPFAC, EDEMA, ALBIN, F1REL)
;   7. Semi-mechanistic dual PD in $ERROR: Cu(t), luminal SGLT2 occupancy (%) and TGF/HCT Emax (%)
;
; Units:
;   TIME, TAD, II : hours (h)
;   AMT, DOSEMG   : milligrams (mg) -> via S2 = V2/1000, F = A(2)/S2 is directly in ng/mL (= ug/L)
;   AGE           : years
;   BWT           : actual body weight (kg; when EDEMA > 0, WDRY = BWT / (1 + EDEMA/100))
;   HT            : height in metres (m; values > 3 are automatically converted from cm -> m)
;   SEX           : 1 = male/boy, 0 = female/girl
;   CREAT         : serum creatinine in umol/L (set to -99 if unavailable)
;   EGFRIN        : direct eGFR in mL/min/1.73m2 (set to -99 to compute from Rhodin/CREAT)
;   ALBIN         : serum albumin in g/L (set to -99 for physiological McNamara age equation)
;   EDEMA         : percentage fluid overload above dry weight (%; 0 = euvolaemic)
;   RFFAC         : cardiorenal GFR preservation (fraction; 1.0 = normal, e.g. 0.80 in stable HF)
;   HEPFAC        : hepatic UGT1A9 preservation in venous congestion (fraction; 1.0 = normal)
;   F1REL         : relative enteral bioavailability (fraction; 1.0 = 100%)
;   MATOPT        : UGT1A9 ontogeny (1 = Farhan 2024, 2 = Badee 2019, 3 = Simcyp s1)
;   TAU           : dosing interval in hours (e.g. 24) used for AUCU24 calculation
; =============================================================================

$INPUT ID TIME TAD DAY EVID AMT ADDL II CMT DV MDV AGE BWT HT SEX CREAT EGFRIN ALBIN EDEMA RFFAC HEPFAC F1REL MATOPT DOSEMG TAU

$DATA example_child_hf.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

$PK
  ; --- 0. Euvolaemic dry weight, height, BMI, and Mosteller BSA ---
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

  ; --- 1. Body composition: FFM (Al-Sallami 2015 / Janmahasatian 2005) & NFM ---
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

  ; --- 2. Plasma albumin ontogeny (McNamara & Alcorn 2002) & protein binding fu ---
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

  ; --- 3. Healthy age-expected GFR (Rhodin 2009) & individual eGFR ---
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
      ; Adult CKD-EPI 2021 (race-free)
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
        ; Paediatric Pierce 2021 CKiD U25 (1 - <18 years)
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
        ; Infants (<1 year) or missing creatinine: Rhodin PMA maturation
        EGFR0 = NGFR
        IF (EGFR0.LT.5.0) EGFR0 = 5.0
        IF (EGFR0.GT.180.0) EGFR0 = 180.0
      ENDIF
    ENDIF
  ENDIF

  ; Cardiorenal GFR preservation (RFFAC)
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

  ; --- 4 & 5. UGT1A9 ontogeny (MAT) & three-pathway clearance decomposition ---
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
      ; Default (MATOPT = 1): Farhan et al. 2024 in-vivo UGT1A9
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

  ; --- 6. Volumes (incl. extracellular oedema in V2), Q, and GI absorption maturation ---
  TVV2 = THETA(5) * FSZV2 + VEDEM
  TVV3 = THETA(6) * (WDRY / 84.0) * FURAT
  TVQ  = THETA(7) * FSZCL
  FGI  = THETA(9) + (1.0 - THETA(9)) * (AGE / (THETA(10) + AGE))
  TVKA = THETA(1) * FGI

  ; --- Individual PK parameters (ADVAN4 TRANS4) ---
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

  ; Analytical steady-state 24-hour unbound exposure & systemic Emax markers
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
  5.44  FIX  ; 2  CL_hep_UGT1A9 (L/h) at ref-NFM, MAT=1, fu=0.086
  15.78 FIX  ; 3  CL_ren_UGT1A9 (L/h) at ref-NFM, MAT=1, RF=1, fu=0.086
  0.011 FIX  ; 4  CL_filt_slope (L/h per mL/min/1.73m2 at BSA=2.0 m2)
  73.47 FIX  ; 5  V2_ref (L at NFM_V2 = 61.68 kg) - Penland 2023 THETA(3)
  168.5 FIX  ; 6  V3_ref (L at BWT = 84 kg, fu = 0.086) - Penland 2023 THETA(4)
  9.137 FIX  ; 7  Q_ref (L/h at NFM_CL = 63.08 kg) - Penland 2023 THETA(5)
  0.5   FIX  ; 8  ALAG1 absorption lag time (h)
  0.45  FIX  ; 9  GI_F0 relative neonatal gastric emptying rate at birth
  0.10  FIX  ; 10 GI_T50 maturation half-time of gastric emptying (years)
  0.09  FIX  ; 11 TM50 UGT1A9 maturation (years) - Farhan et al. 2024
  1.9   FIX  ; 12 HILL UGT1A9 maturation - Farhan et al. 2024
  0.05  FIX  ; 13 FBIRTH UGT1A9 fraction at birth - Farhan et al. 2024

; Default set to 0 FIX for exact 1-to-1 deterministic reproduction of the example child.
; For population Monte Carlo simulation with inter-individual variability (IIV):
; replace the three lines below with:
;   0.36   FIX  ; ETA(1) BSV on KA (clinically calibrated; or 8.52 FIX from Penland 2023)
;   0.1136 FIX  ; ETA(2) BSV on CL/F
;   0.8229 FIX  ; ETA(3) BSV on V3/F
$OMEGA
  0 FIX  ; ETA(1) BSV KA (deterministic mode: ETA = 0)
  0 FIX  ; ETA(2) BSV CL
  0 FIX  ; ETA(3) BSV V3

$SIGMA
  0 FIX  ; No residual error (pure PK/PD prediction)

$SIMULATION (20261125) ONLYSIM SUBPROBLEMS=1

$TABLE ID TIME TAD DAY EVID AMT WDRY FFM ALBU FU NGFR EGFR RF MAT
       CL CLHEP CLRUGT CLFILT V2 V3 Q KA CTOT CU SGLT2 AUCU24 ETGF EHCT
       NOAPPEND NOPRINT ONEHEADER FORMAT=s1PE14.6 FILE=sdtab_dapa_ped_hf.tab
```

---

### 6.3 Example 1-child dataset for NONMEM (`example_child_hf.csv`)

Save the CSV below as `example_child_hf.csv` next to `dapa_ped_hf_sim.ctl` (`9-month-old` infant, actual weight `8.40 kg` with `5%` fluid overload → dry weight `8.00 kg`, `RFFAC = 0.80`, `HEPFAC = 0.90`, `ALBIN = 32 g/L`, `1.20 mg` q24h for 7 days). Expected NONMEM output on `ID = 1`: `CL = 4.564798 L/h`, `V2 = 8.458278 L`, `V3 = 21.37610 L`, `Q = 1.726980 L/h`, `KA = 2.649688 h⁻¹`, `FU = 0.114556`, `AUCU24 = 30.11452 ng·h/mL`, `CTOT(t=1h, Day 1) = 83.623496 ng/mL`, `CTOT(t=145h, Day 7) = 85.259729 ng/mL`:

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

## 7. Paediatric Heart Failure Extension: Pathophysiology & Dual PK–PD (`R/hf_pkpd.R`)

For dose prediction in **paediatric heart failure** (from neonates to adolescents), the baseline physiological model can optionally be extended with heart failure-specific pathophysiology (defaulting to neutral baseline values `1.0` / `0%`) and semi-mechanistic PD readouts:

1. **Extracellular fluid overload vs. euvolaemic dry weight (`edema_pct`):**
   Given fluid overload percentage \(\text{FO}\) (`%` above dry weight):
   $$W_{\text{dry}} = \frac{\text{BWT}}{1 + \text{FO}/100}, \qquad V_{\text{edema}} = \text{BWT} - W_{\text{dry}}\text{ (L)}$$
   Because oedema fluid contains no hepatocytes or renal tubules, `FFM`, `NFM_CL`, `NFM_V2`, and lipophilic \(V_3\) are calculated from \(W_{\text{dry}}\), whereas the central/interstitial distribution volume expands by \(+V_{\text{edema}}\).
2. **Cardiorenal syndrome (`rf_factor`), hepatic congestion (`cl_hep_factor`), hypoalbuminaemia (`alb_override`), and enteral bioavailability (`f1_rel`):**
   - \(\text{eGFR}_{\text{eff}} = \min\!\left(180, \max(5, \; \text{eGFR} \cdot f_{\text{RF}})\right)\), reducing both proximal tubular UGT1A9 clearance (`15.78 · MAT · Fsize_CL · RF`) and passive filtration in renal hypoperfusion (including in infants `< 1 y` where `copulapop` does not use creatinine).
   - \(CL_{\text{hep}} = f_{u,\text{ratio}} \cdot 5.44 \cdot \text{MAT} \cdot F_{\text{size},CL} \cdot f_{\text{hep}}\) accounts for impaired hepatic glucuronidation under elevated central venous pressure (right-sided failure / Fontan).
   - Manually specified serum albumin raises \(f_u\) in hypoalbuminaemia or *protein-losing enteropathy*, lowering total \(C(t)\) and \(\text{AUC}_\tau\) while leaving pharmacologically active unbound \(C_u(t) = f_u C(t)\) and \(\text{AUC}_{u,24\text{h}}\) unchanged.
   - Mean relative bioavailability \(F_{1,\text{rel}}\) scales the administered oral/enteral dose for feeding-tube or suspension loss.
3. **Semi-mechanistic dual PK–PD relationships:**
   - **Proximal tubular luminal SGLT2 receptor occupancy:** from human SGLT2 \(\text{IC}_{50} = 1.1\text{ nM} = 0.45\text{ ng/mL}\) and proximal S1/S2 water-reabsorption concentration factor \(\phi_{\text{PT}} = 2.25\) (Yakovleva et al. 2014), the effective unbound plasma concentration for 50% luminal occupancy is \(C_{u,50} = 0.20\text{ ng/mL}\):
     $$O_{\text{SGLT2}}(t) = \frac{100 \cdot C_u(t)}{0.20 + C_u(t)}\text{ (\%)}, \qquad \overline{O}_{\text{SGLT2}} = \frac{1}{\tau}\int_0^\tau O_{\text{SGLT2}}(t)\,dt$$
   - **GFR-coupled effective renal osmotic/natriuretic PD index:** scales SGLT2 blockade by filtered solute delivery to the proximal tubule:
     $$\text{PD}_{\text{renal,eff}} = \overline{O}_{\text{SGLT2}} \cdot \min\!\left(1.5, \; \frac{\text{eGFR}_{\text{eff}}}{100}\right)$$
   - **Systemic/cardiorenal pleiotropic \(E_{\max}\) responses (Kroonen et al. 2020):** derived from published clinical exposure–response models converted to unbound 24-hour steady-state exposure (\(\text{EAUC}_{u,50,\text{TGF}} = 5.2\text{ ng·h/mL}\), \(\text{EAUC}_{u,50,\text{HCT}} = 20.7\text{ ng·h/mL}\)).

---

## 8. References

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
14. Kroonen MYAM et al. *Br J Clin Pharmacol* 2020;86:1905–1916 — [doi:10.1111/bcp.14318](https://doi.org/10.1111/bcp.14318) (dapagliflozin exposure–response on cardiorenal markers).
15. Yakovleva T et al. *Front Pharmacol* 2014;5:218 — [doi:10.3389/fphar.2014.00218](https://doi.org/10.3389/fphar.2014.00218) (semi-mechanistic luminal SGLT2 modelling).

## 9. Limitations

- Extrapolation below ~1 year relies on mechanistic physiological ontogeny (UGT1A9, GFR,
  albumin, GI absorption) and PBPK verification; published clinical PK data of dapagliflozin
  in neonates and young infants remain limited.
- Virtual `copulapop` cohorts represent general or hospitalised populations; the heart failure
  controls (`rf_factor`, `cl_hep_factor`, `edema_pct`, `alb_override`) provide mechanistic
  scenario analysis rather than replacing prospective clinical PK/PD validation in paediatric HF.
- Research and educational simulation model — not validated clinical dosing advice.
