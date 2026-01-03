# GHBVS
Code for: 

> Chen, Y., Griffin, J. and Rubio, F.J. (2025). Bayesian variable and hazard structure selection in the General Hazard model. Submitted.

- `simulation.R`: Reproducing the simulation results in the paper.
- `flchain.R`: Reproducing the FLC application.
- `nki70.R`: Reproducing the NKI70 application. 

Default setting: C++ implementation of likelihood, gradient, and Hessian with optimisation in R (`BVS.MCMC(..., use_cpp = F)`); for additional C++ implementation of optimisation (currently under working), set `BVS.MCMC(..., use_cpp = T)`.

Adjust the hyperparameters of the priors and sampler according to the comments if needed. 
