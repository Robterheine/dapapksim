# About this app

This Shiny app (**dapapksim**) was created by **Rob ter Heine**
(e-mail: [r.terheine@radboudumc.nl](mailto:r.terheine@radboudumc.nl)), of the
**Radboud Applied Pharmacometrics group**:
[radboudumc.nl/en/research/research-groups/applied-pharmacometrics](https://www.radboudumc.nl/en/research/research-groups/applied-pharmacometrics)

The latest version is available at:
[github.com/Robterheine/dapapksim](https://github.com/Robterheine/dapapksim)

## What does this app do?

Monte Carlo population simulation and physiological extrapolation of dapagliflozin
pharmacokinetics and pharmacodynamics (based on Penland et al. 2023,
[doi:10.1002/jcph.2196](https://doi.org/10.1002/jcph.2196), extended with `FFM`/`NFM`
body-composition allometry, plasma albumin-binding ontogeny, Rhodin/CKiD-U25 renal
function, hepatic + renal proximal tubular UGT1A9 maturation, and paediatric heart failure
pathophysiology) in virtual populations from neonates to older adults generated with
[copulapop](https://github.com/Robterheine/copulapop).
See the **"Model & assumptions"** tab for the complete equations and references, and
**`docs/VALIDATION_REPORT.md`** in the repository for the quantitative validation against
published adult and paediatric PK studies.

Research and educational simulation tool — not clinical dosing advice.
