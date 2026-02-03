# GHBVS
Code for: 

> Chen, Y., Griffin, J. and Rubio, F.J. (2026). Bayesian variable and hazard structure selection in the General Hazard model. Submitted.

- `simulation.R`: Example code to carry out simulation. To reproduce the same result, run `mpi.R` for RNG stream controls.
- `flchain.R`: Reproducing the FLC application.
- `nki70.R`: Reproducing the NKI70 application.

Default setting: Pure R implementation (`BVS.MCMC(..., use_cpp = F)`), used to produce the results in the paper.
For C++ implementation of optimisation, likelihood, hessian, etc. (currently under working, DO NOT use it), set `BVS.MCMC(..., use_cpp = T)`.

For normal use, adjust the hyperparameters of the priors and sampler according to the comments as needed. 
