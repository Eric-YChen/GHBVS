# GHBVS: Bayesian Variable and Hazard Structure Selection in the General Hazard Model

## Overview

This repository contains the R code and data used to reproduce the results in:

> Chen, Y., Griffin, J. and Rubio, F.J. (2026). Bayesian variable and hazard
> structure selection in the General Hazard model. *Submitted.*
> [arXiv:2602.03756](https://arxiv.org/abs/2602.03756)

Standard Bayesian variable selection approaches for survival data typically assume a fixed hazard structure — most commonly Proportional Hazards (PH) or Accelerated Failure Time (AFT). This paper proposes a **joint Bayesian framework for simultaneous selection of relevant covariates and the hazard structure** within the **General Hazard (GH)** model, which contains PH, AFT, and Accelerated Hazards (AH) as special cases. The method is implemented via MCMC using the `BVS.MCMC()` function.

## Requirements

The following R packages are required:

```r
install.packages(c("survival", "mvtnorm", "coda", "Rcpp"))
```

> Verify exact dependencies against the `library()` calls at the top of each
> script, as the list above may be incomplete.

## Repository structure

```
GHBVS/
├── BVS.MCMC.R          # Core MCMC function for joint variable and
│                       #   hazard structure selection
├── simulation.R        # Example simulation code
├── mpi.R               # RNG stream controls for reproducible simulation results
├── flchain.R           # Real data application: FLC (flchain) dataset
├── nki70.R             # Real data application: NKI70 breast cancer dataset
└── [data files]        # Deidentified datasets for the real data applications
```

## Usage

The main function is `BVS.MCMC()`. Two implementations are available:

| Setting | Argument | Status | Notes |
|---|---|---|---|
| Pure R | `use_cpp = FALSE` | ✓ Default | Used to produce all results in the paper |
| C++ (optimisation, likelihood, Hessian) | `use_cpp = TRUE` | ⚠ In development — **do not use** | — |

For the simulation study, run `mpi.R` first to set RNG stream controls before
running `simulation.R`, in order to reproduce the exact results reported in the
paper. Hyperparameters for priors and the sampler can be adjusted via the
function arguments — see comments within the scripts for guidance.

## Data applications

| Script | Dataset | Description |
|---|---|---|
| `flchain.R` | FLC (`flchain`) | Free light chain study |
| `nki70.R` | NKI70 | Breast cancer gene expression (70-gene signature) |

Deidentified data for both applications are included in the repository.

## Related resources

- [HazReg](https://github.com/FJRubio67/HazReg) — R package for parametric
  GH regression (maximum likelihood); the GH framework underlying this work
- [BVSSurv](https://github.com/FJRubio67/BVSSurv) — short course on Bayesian
  variable selection for survival data, covering related methodology

## Citation

If you use this code, please cite:

```bibtex
@article{chen2026ghbvs,
  author  = {Chen, Y. and Griffin, J. and Rubio, F.J.},
  title   = {Bayesian variable and hazard structure selection in the
             General Hazard model},
  journal = {Submitted},
  year    = {2026},
  note    = {arXiv:2602.03756}
}
```

