library(mvtnorm)
library(survival)
library(glmnet)
library(HazReg)
library(parallel)

set.seed(1234)
nki70 <- read.csv("./data/nki70.csv", header = TRUE, stringsAsFactors = T)
X <- nki70[,-c(1:2)]
X$Diam <- as.numeric(X$Diam) - 1
X$N <- as.numeric(X$N)
X$ER <- as.numeric(X$ER) - 1
X$Grade <- as.numeric(X$Grade)
X <- apply(X,2,scale)
dat <- list(y = nki70$time, status = nki70$event, x=X, xt=X)

########################################################
##### Proposed Method
########################################################
source("ADS_sampler.R")

n.chains <- 10                            # independent chains
n.cores  <- min(n.chains, detectCores())  # one worker per chain
rng.seed <- 1234
out.file <- "nki_new_chains.RData"
log.dir  <- "logs"
wd       <- getwd()

# n.chain stays 1: the chain loop inside BVS.MCMC is unfinished, so each worker
# runs a single chain and the parallelism is over workers instead.
mcmcpar <- list(n.iter = 110000, n.burnin = 10000, n.thin = 100,
                n.chain = 1, printevery = 10000)
gamparam <- list(layer1.betabin.par = c(1,1),
                 layer2.mulnom.par = c(0.25,0.25,0.25,0.25), # (AH,PH,GH,AFT)
                 layer3.beta.par = 1/3)
gprior.mat    <- "LCM"    # one of ("LCM", "LCM-hyper", "Product")
gamprior.type <- "default"

# One chain, run in its own process. BVS.MCMC keeps state in global variables
# (model.dict, move.list, current_move, gamprior, gam.param), so chains must not
# share a session -- separate processes give each chain its own copy.
run.chain <- function(chain.id){
  setwd(wd)
  library(mvtnorm); library(survival); library(glmnet); library(HazReg)
  source("ADS_sampler.R")
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng.seed)
  seed.i <- .Random.seed
  for (i in seq_len(chain.id)) seed.i <- parallel::nextRNGStream(seed.i)
  assign(".Random.seed", seed.i, envir = .GlobalEnv)

  con <- file(file.path(log.dir, sprintf("chain%02d.log", chain.id)), open = "wt")
  sink(con); sink(con, type = "message")
  on.exit({sink(type = "message"); sink(); close(con)}, add = TRUE)

  gam.init <- init.gam(dat, history.keep = F)

  modelpar <- list(inits = list(nu = rnorm(1), theta0 = rnorm(1),
                                theta = rnorm(ncol(dat$xt),0,0.5),
                                eta = rnorm(ncol(dat$x),0,0.5),
                                g.theta = 1, g.eta = 1, g.aft = 1, # for Product prior
                                g = 1, # for LCM prior
                                a_g = NULL, # for hyper-g LCM prior, NULL for robust prior, 2 < a_g <= 4 for hyper-g/n prior
                                gam = gam.init,
                                nu.alpha = 0.01, nu.beta = 0.01),
                   gam.param = gamparam,
                   data = dat)

  t0 <- Sys.time()
  chains <- tryCatch(BVS.MCMC(modelpar, mcmcpar, verbose = T,
                              gprior.mat = gprior.mat,
                              gamprior = gamprior.type),
                     error = function(e) list(error = conditionMessage(e)))
  chains$chain.id <- chain.id
  chains$gam.init <- gam.init
  chains$runtime  <- difftime(Sys.time(), t0, units = "mins")
  chains
}

dir.create(log.dir, showWarnings = FALSE)
cat("Running", n.chains, "chains on", n.cores, "cores; per-chain logs in", log.dir, "\n")
t.start <- Sys.time()

cl <- makeCluster(n.cores)
chains.list <- tryCatch({
  clusterExport(cl, c("dat", "mcmcpar", "gamparam", "gprior.mat",
                      "gamprior.type", "wd", "log.dir", "rng.seed"))
  parLapplyLB(cl, seq_len(n.chains), run.chain, chunk.size = 1)
}, finally = stopCluster(cl))
names(chains.list) <- paste0("chain", seq_len(n.chains))

cat("Wall time:", format(difftime(Sys.time(), t.start, units = "mins")), "\n")
failed <- which(sapply(chains.list, function(x) !is.null(x$error)))
if (length(failed)){
  cat("FAILED chains:", paste(names(chains.list)[failed], collapse = ", "), "\n")
  for (i in failed) cat("  ", names(chains.list)[i], ":", chains.list[[i]]$error, "\n")
}

save(chains.list, mcmcpar, gamparam, gprior.mat, gamprior.type, dat,
     file = out.file)
cat("Saved", length(chains.list), "chains to", out.file, "\n")

##################################################################
# result summary
##################################################################
vars <- colnames(dat$x)
lse <- function(x){ m <- max(x); m + log(sum(exp(x - m))) }
role <- c("1" = "^time", "2" = "^haz", "3" = "", "4" = "")
gam2label <- function(s){
  g <- strsplit(s, "")[[1]]; i <- which(g != "0")
  paste0("{", paste0(vars[i], role[g[i]], collapse = ", "), "}^", gam2haz(s))
}
draws <- lapply(chains.list, function(x)
  apply(x$chain[[1]], 1, paste0, collapse = ""))
visited <- unique(unlist(draws))

# log-posterior of each visited model from the model dictionaries
lp <- do.call(rbind, lapply(chains.list, function(x)
  data.frame(gam = x$model.prob$gam, lp = as.numeric(x$model.prob$postprob))))
lp <- lp[is.finite(lp$lp), ]
cat("max cross-chain |log-post diff| for the same model:",
    signif(max(tapply(lp$lp, lp$gam, function(v) diff(range(v)))), 3), "\n")
lpv <- as.vector(tapply(lp$lp, lp$gam, function(v) v[1])[visited])
stopifnot(!anyNA(lpv))

res <- data.frame(model = visited, lp = lpv)
res$pp  <- exp(res$lp - lse(res$lp))                          # pooled, normalised
res$mf  <- as.vector(table(factor(unlist(draws), levels = visited))) / length(unlist(draws))
res$haz <- sapply(res$model, gam2haz)
res$size <- sapply(res$model, function(s) sum(strsplit(s, "")[[1]] != "0"))
res$label <- sapply(res$model, gam2label)
res <- res[order(-res$pp), ]
res$rank <- seq_len(nrow(res))

per.chain <- sapply(draws, function(d){
  u <- unique(d); l <- res$lp[match(u, res$model)]
  setNames(exp(l - lse(l)), u)[res$model[1:10]]
})
per.chain[is.na(per.chain)] <- 0
cat("\n-- top 15 models: pooled normalised post. prob (pp), MCMC frequency (mf)\n")
top <- res[1:15, c("rank", "label", "size", "pp", "mf")]
top$pp <- round(top$pp, 3); top$mf <- round(top$mf, 3)
print(top, row.names = FALSE)