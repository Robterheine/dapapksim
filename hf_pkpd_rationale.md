## Wetenschappelijke rationale: PK, PD en dosisextrapolatie bij pediatrisch hartfalen

Bij het voorspellen van een dapagliflozinedosering voor **hartfalen bij jonge kinderen** (pasgeborenen, zuigelingen en peuters) vanuit volwassen studies (DAPA-HF, DELIVER: `10 mg` 1×/dag) vereist het pediatrische extrapolatiekader (**ICH E11A**) een kritische toets op zowel **farmacokinetische (PK)** als **farmacodynamische (PD)** aannames. Deze module maakt beide assen expliciet zichtbaar en kwantificeerbaar.

---

### 1. Farmacokinetische (PK) rationale specifiek voor pediatrisch hartfalen

1. **Sturen op ongebonden blootstelling (`AUCu,24h,ss`), nóóit op totale plasmaconcentratie:**
   - Dapagliflozine is voor `91,4%` aan albumine gebonden (`fu_adult = 0.086`) en heeft een zeer lage extractieratio (`ER < 0.05`). Totale schijnbare klaring \(CL/F\) schaalt daarom lineair met `fu_ratio = fu / 0.086`, terwijl **ongebonden klaring** \(CL_u/F = (CL/F)/f_u\) uitsluitend wordt bepaald door intrinsieke metabole en glomerulaire capaciteit.
   - Pasgeborenen en jonge zuigelingen hebben fysiologisch lager serumalbumine (McNamara & Alcorn 2002), en kinderen met ernstig hartfalen hebben vaak **hypoalbuminemie** door veneuze stuwing, cardiale cachexie of *protein-losing enteropathy* (PLE, bv. bij Fontan-circulatie). Hierdoor daalt de *totale* plasmaconcentratie \(C(t)\), terwijl de farmacologisch actieve *vrije* concentratie \(C_u(t) = f_u C(t)\) gelijk blijft. Doseren of TDM-bijsturen op totale `AUCtau` zou bij hypoalbuminemie leiden tot **15–40% overdosering**.
2. **Drie-paden klaring en UGT1A9-ontogenie-onzekerheid onder 1 jaar:**
   - Ongeveer `73%` van de volwassen dapagliflozineklaring verloopt via **renale proximale tubulaire UGT1A9-glucuronidering** (`15.78 L/h`), `~25%` via **hepatische UGT1A9** (`5.44 L/h`) en slechts `< 5%` via onveranderde glomerulaire filtratie (`0.011 × eGFR × BSA/2`).
   - Boven `1–2 jaar` convergeren de drie literatuur-ontogeniefuncties (`Farhan 2024`, `Badee 2019` en `Simcyp s1`): door allometrische schaling op vetvrije massa (\((\text{NFM}_{CL}/63.08)^{0.75}\)) is de klaring *per kg lichaamsgewicht* bij peuters en kleuters hoger dan bij volwassenen, zodat **`~0.10 mg/kg/dag`** overeenkomt met een volwassen `5 mg`-blootstelling (`AUCu,24h ≈ 20 ng·h/mL`) en **`~0.18–0.20 mg/kg/dag`** met een volwassen `10 mg`-blootstelling (`AUCu,24h ≈ 38–48 ng·h/mL`).
   - **Onder 6–12 maanden** lopen snelle *in-vivo* rijping (`Farhan 2024`, afgeleid uit hepatisch propofolmetabolisme) en conservatieve *in-vitro* ontogenie (`Simcyp s1` uit Jo et al. 2021) echter een **factor 3 tot 4 uiteen** (bv. op `2 maanden`: `0.20 mg/kg` vs. `0.05 mg/kg` voor dezelfde `AUCu,24h`). Bovendien is niet klinisch bewezen dat renaal tubulair UGT1A9 exact even snel rijpt als hepatisch UGT1A9. Dit dwingt bij jonge zuigelingen tot een conservatieve **"start low, go slow"** strategie op basis van het `Simcyp s1`-scenario.
3. **Hartfalen-pathofysiologie: cardiorenaal syndroom, leverstuwing, oedeem en sonde-toediening:**
   - **Cardiorenale hypoperfusie (`rf_factor`):** Standaard kent de virtuele populatie aan zuigelingen `< 1 jaar` een gezonde leeftijds-GFR toe (`RF = eGFR / nGFR ≈ 1.0`), en bij oudere kinderen met spieratrofie/cachexie kan een laag serumcreatinine de werkelijke GFR via `CKiD U25` zelfs overschatten. Bij *low cardiac output* en lisdiuretica/ACEi-gebruik daalt de werkelijke GFR (`rf_factor < 1.0`), waardoor zowel tubulaire UGT1A9-klaring als passieve filtratie evenredig afnemen.
   - **Hepatische stuwing (`cl_hep_factor`):** Verhoogde centraal-veneuze druk (rechtszijdig/biventriculair falen, Fontan) vermindert hepatische glucuronidering (`cl_hep_factor < 1.0`; bij matig-ernstige leverfunctiestoornis stijgt dapagliflozine-`AUC` klinisch tot `+40–67%`).
   - **Extracellulair oedeem (`edema_pct`):** Bij overvulling (\(W_{\text{actueel}} = W_{\text{droog}}(1 + \text{FO}/100)\)) vergroot extracellulair vocht wél het centrale verdelingsvolume \(V_2\) (`+1 L/kg` oedeemvocht), maar voegt **geen metabool actief lever- of nierparenchym** toe. Zowel `FFM`/`NFM` als `mg/kg`-doseringen horen daarom op het **euvolemische drooggewicht** te worden gebaseerd.
   - **Pediatrische toedieningsvorm (`f1_rel` & `IOV`):** Verpulverde/gedispergeerde tabletten of suspensie via een neusmaagsonde en splanchnisch darmwandoedeem kunnen zowel een lagere gemiddelde biologische beschikbaarheid (`f1_rel < 1.0`) als hogere dag-tot-dag variabiliteit (`IOV %CV`) geven.

---

### 2. Kritische farmacodynamische (PD) rationale: duale PK–PD vs. veiligheidsgrenzen

Waarom volstaat blinde PK *exposure-matching* niet zonder PD-interpretatie? Omdat dapagliflozine bij hartfalen via **twee verschillende farmacodynamische assen** werkt die bij een lage of onrijpe GFR **tegengesteld** reageren:

| Mechanistische as | Drijvende concentratie / blootstelling | Modelvergelijking in `dapapksim` | Gedrag bij lage GFR (neonaat of cardiorenaal HF) |
|---|---|---|---|
| **1. Luminale SGLT2-receptorbezetting (proximale tubulus S1/S2)** | Gefiltreerde vrije fractie in tubuluslumen: \(C_{\text{lumen}}(t) \approx \phi_{\text{PT}} \cdot f_u C(t)\) met \(\phi_{\text{PT}} = 2.25\) en \(\text{IC}_{50} = 0.45\text{ ng/mL}\) (`1.1 nM`) \(\Rightarrow C_{u,50} = 0.20\text{ ng/mL}\) | $$O_{\text{SGLT2}}(t) = \frac{C_u(t)}{0.20 + C_u(t)}$$ *(Yakovleva et al. 2014; FDA NDA 202293)* | **Blijft hoog of stijgt (`~80–92%`):** bij gematchte `AUCu,24h` is SGLT2 al grotendeels verzadigd (*ceiling effect*). |
| **2. Effectieve renale osmotische / natriuretische PD-drive** | Product van SGLT2-blokkade én glomerulair gefiltreerde natrium-/glucosevracht (`eGFR` t.o.v. rijpe referentie `100 mL/min/1.73m²`) | $$\text{PD}_{\text{renaal,eff}} = \overline{O}_{\text{SGLT2}} \cdot \min\!\left(1.5, \frac{\text{eGFR}}{100}\right)$$ | **Daalt sterk (`~25–45%`):** weinig filtraat bereikt S1/S2. **Valkuil:** hoger doseren geeft *géén* extra diurese (SGLT2 is al bezet!), maar wel overblootstelling. |
| **3. Systemische / cardiorenaal-pleiotrope respons (TGF & HCT/erytropoëse)** | Systemische ongebonden 24-uurs blootstelling \(\text{AUC}_{u,24\text{h,ss}}\) (o.a. NHE1/NHE3, keton-energetica, erytropoëse, afferente tonus) | $$E = \frac{\text{AUC}_{u,24\text{h,ss}}}{\text{EAUC}_{u,50} + \text{AUC}_{u,24\text{h,ss}}}$$ met \(\text{EAUC}_{u,50,\text{TGF}} = 5.2\) en \(\text{EAUC}_{u,50,\text{HCT}} = 20.7\text{ ng·h/mL}\) *(Kroonen et al. 2020)* | **Stijgt bij dalende klaring:** verklaart waarom cardiovasculaire/renale bescherming bij volwassenen behouden blijft tot `eGFR 20–25`, ook zonder sterke glucosurie. |

#### Kritische kanttekeningen voor de kliniek bij jonge kinderen:
- **Ontbreken van volwassen dosis-responsdata in HF:** In DAPA-HF en DELIVER is alleen `10 mg` 1×/dag onderzocht (`AUCu,24h,ss` mediaan `~38.5 ng·h/mL` bij gezonde nierfunctie tot `~48.2 ng·h/mL` bij de gemiddelde DAPA-HF `eGFR` van `66 mL/min/1.73m²`; `P5–P95` bij milde-matige nierfunctiestoornis `~24–75 ng·h/mL`). Een volwassen `5 mg` dosis (`AUCu,24h,ss ≈ 19–22 ng·h/mL`) bereikt al `~79%` 24-uurs SGLT2-remming en `~80%` TGF-respons en vormt daarom een rationeel **veilig startdoel**.
- **Verschoven veiligheids-PD bij zuigelingen (`< 1–2 jaar`):**
  1. **Normoglykemische glucosurie → hypoglykemie & ketose:** Zuigelingen hebben kleine hepatische glycogeenvoorraden en een hoge glucosebehoefte per kg. Hanteer strikte *sick-day rules* (direct pauzeren bij vasten, braken, gastro-enteritis of verminderde drinkintake).
  2. **Interactie met lisdiuretica → volumedepletie & AKI:** De TGF-gemedieerde daling in intraglomerulaire druk telt op bij hooggedoseerde furosemide + ACEi/ARNI en vochtbeperking.
  3. **Calorisch verlies (*failure to thrive*):** Volg gewichts- en lengtegroei nauwgezet bij zuigelingen met cardiale cachexie of grote links-rechts shunts.

> **Disclaimer:** Dit model ondersteunt wetenschappelijke hypothesevorming, studie-ontwerp en farmacometrische scenario-analyse (*Model-Informed Drug Development*), en vervangt geen geregistreerd klinisch doseringsadvies.
