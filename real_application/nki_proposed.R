# Usage (from real_application/): Rscript nki_proposed.R <LCM|Product>
suppressMessages({library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
gprior.mat <- if (length(args) >= 1) args[1] else "LCM"
stopifnot(gprior.mat %in% c("LCM", "Product"))
quick <- nzchar(Sys.getenv("BVS_QUICK"))
res.dir <- if (quick) "results_quick" else "results"
dir.create(res.dir, showWarnings = FALSE)

n.chains <- if (quick) 2 else 10
n.cores  <- min(n.chains, detectCores())
rng.seed <- 1234
mcmcpar  <- if (quick) list(n.iter = 400, n.burnin = 200, n.thin = 2, n.chain = 1, printevery = 200) else
  list(n.iter = 110000, n.burnin = 10000, n.thin = 100, n.chain = 1, printevery = 10000)
log.dir  <- file.path(res.dir, paste0("logs_nki_", gprior.mat))
dir.create(log.dir, showWarnings = FALSE)
wd <- getwd()

run.chain <- function(chain.id){
  setwd(wd)
  BVS_CORE <<- normalizePath("../core")
  suppressMessages({library(mvtnorm); library(survival); library(glmnet)})
  source(file.path(BVS_CORE, "ADS_sampler.R"))
  source("prepare_data.R")
  dat <- prep.nki()

  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng.seed)
  seed.i <- .Random.seed
  for (i in seq_len(chain.id)) seed.i <- parallel::nextRNGStream(seed.i)
  assign(".Random.seed", seed.i, envir = .GlobalEnv)

  con <- file(file.path(log.dir, sprintf("chain%02d.log", chain.id)), open = "wt")
  sink(con); sink(con, type = "message")
  on.exit({sink(type = "message"); sink(); close(con)}, add = TRUE)

  gam.init <- init.gam(dat, history.keep = FALSE)
  modelpar <- list(inits = list(nu = rnorm(1), theta0 = rnorm(1),
                                theta = rnorm(ncol(dat$xt), 0, 0.5), eta = rnorm(ncol(dat$x), 0, 0.5),
                                g.theta = 1, g.eta = 1,
                                g = 1,
                                a_g = NULL, gam = gam.init,
                                nu.alpha = 0.01, nu.beta = 0.01),
                   data = dat)
  t0 <- Sys.time()
  chains <- tryCatch(BVS.MCMC(modelpar, mcmcpar, verbose = TRUE,
                              gprior.mat = gprior.mat, gamprior = "default"),
                     error = function(e) list(error = conditionMessage(e)))
  chains$chain.id <- chain.id
  chains$gam.init <- gam.init
  chains$runtime  <- difftime(Sys.time(), t0, units = "mins")
  chains
}

cat("Running", n.chains, "chains on", n.cores, "cores; logs in", log.dir, "\n")
t.start <- Sys.time()
cl <- makeCluster(n.cores)
chains.list <- tryCatch({
  clusterExport(cl, c("mcmcpar", "gprior.mat", "wd", "log.dir", "rng.seed"))
  parLapplyLB(cl, seq_len(n.chains), run.chain, chunk.size = 1)
}, finally = stopCluster(cl))
names(chains.list) <- paste0("chain", seq_len(n.chains))
cat("Wall time:", format(difftime(Sys.time(), t.start, units = "mins")), "\n")
failed <- which(sapply(chains.list, function(x) !is.null(x$error)))
if (length(failed)) for (i in failed) cat("FAILED", names(chains.list)[i], ":", chains.list[[i]]$error, "\n")

source("prepare_data.R"); dat <- prep.nki()
save(chains.list, mcmcpar, gprior.mat, dat, file = file.path(res.dir, paste0("nki_", gprior.mat, ".RData")))
cat("Saved", length(chains.list), "chains to", file.path(res.dir, paste0("nki_", gprior.mat, ".RData")), "\n")
