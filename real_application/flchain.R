# Usage (from real_application/): Rscript flchain.R [orig|trans] [LCM|Product|LCM-hyper]   (no arguments: all six runs)
BVS_CORE <- normalizePath("../core")
suppressMessages({library(mvtnorm); library(survival); library(glmnet)})
source(file.path(BVS_CORE, "ADS_sampler.R"))
source("prepare_data.R")

args <- commandArgs(trailingOnly = TRUE)
designs <- if (length(args) >= 1) args[1] else c("orig", "trans")
priors  <- if (length(args) >= 2) args[2] else c("LCM", "Product", "LCM-hyper")
quick   <- nzchar(Sys.getenv("BVS_QUICK"))

mcmcpar <- if (quick) list(n.iter = 2000, n.burnin = 500, n.thin = 5, n.chain = 1, printevery = 1000) else
  list(n.iter = 50000, n.burnin = 10000, n.thin = 10, n.chain = 1, printevery = 10000)
res.dir <- if (quick) "results_quick" else "results"
dir.create(res.dir, showWarnings = FALSE)

for (design in designs) for (prior in priors){
  set.seed(1234)
  dat <- prep.flchain(design)
  modelpar <- list(inits = list(nu = rnorm(1), theta0 = rnorm(1),
                                theta = rnorm(ncol(dat$xt), 0, 0.5), eta = rnorm(ncol(dat$x), 0, 0.5),
                                g.theta = 1, g.eta = 1,
                                g = 1,
                                a_g = NULL,
                                gam = init.gam(dat, history.keep = FALSE),
                                nu.alpha = 0.01, nu.beta = 0.01),
                   data = dat)
  cat("\n==== FLC design:", design, "| prior:", prior, "====\n")
  t0 <- Sys.time()
  chains <- BVS.MCMC(modelpar, mcmcpar, verbose = TRUE, gprior.mat = prior, gamprior = "default")
  runtime <- difftime(Sys.time(), t0, units = "mins")
  print(runtime)
  print(head(chains$post.prob, 5))
  print(chains$haz.post)
  if (length(chains$opt.log)) { cat("optimisation problems:\n"); print(chains$opt.log) }
  save(chains, dat, modelpar, mcmcpar, design, prior, runtime,
       file = file.path(res.dir, paste0("flchain_", design, "_", prior, ".RData")))
}
