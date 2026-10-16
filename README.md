# Code for *Bayesian variable and hazard structure selection in the General Hazard model*

## Requirements

R (>= 4.2) with the packages `mvtnorm`, `survival`, `glmnet`, `numDeriv` and `MASS`. In addition:

- reference methods for NKI70: `mombf`, `BVSNLP`;
- simulation study on a cluster: `foreach`, `doMPI`, `Rmpi`, `doRNG`;
- figures: `ggplot2`, `scales`, `dplyr`, `tidyr`, `ggh4x`.

## Contents

| Folder | Contents |
|---|---|
| `core/` | Extended ADS sampler and model-space prior (`ADS_sampler.R`), marginal likelihoods under the LCM and Product priors, likelihood routines and helper functions. Sourced by all other scripts. |
| `real_application/` | FLC (`flchain.R`, `flchain_summary.R`), NKI70 (`nki_proposed.R`, `nki_reference.R`, `nki_summary.R`) and the ipilimumab trial (`ipilimumab.R`); data in `data/`. |
| `simulation/` | Simulation studies (`sim_*.R`), summaries (`summarise_sim.R`) and figures (`make_figures.R`). |

Each script states its arguments in its first line. Run scripts from their own folder; results are
written to `results/` within that folder.

## Results in the paper

| Result | Script |
|---|---|
| Section 6, Figure 2; Supplementary Section 10 | `simulation/sim_main.R`, `sim_conv.R`, `sim_bic.R`, `sim_freq.R`, then `summarise_sim.R` and `make_figures.R` |
| Section 7.1, Table 2 | `real_application/flchain.R`, `flchain_summary.R` |
| Section 7.2, Table 3; Supplementary Section 12 | `real_application/nki_proposed.R`, `nki_reference.R`, `nki_summary.R` |
| Supplementary Section 1 | `real_application/ipilimumab.R` |
| Supplementary Section 13 | `simulation/sim_polarised.R` |
| Supplementary Section 14 | `simulation/sim_pm.R` |

The simulation scripts run under MPI (`mpirun -np <N> Rscript sim_main.R ...`).
