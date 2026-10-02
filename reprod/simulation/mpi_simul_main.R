## ============================================================================
##
## Main simulation study
##
## Usage (via srun on Slurm, or mpirun):
##   srun Rscript mpi_simul_main.R <p> [chunk_size]
##
##   srun Rscript mpi_simul_main.R 10             # p = 10,  default chunk
##   srun Rscript mpi_simul_main.R 50   300       # p = 50,  chunk_size 300
##   srun Rscript mpi_simul_main.R 100  600       # p = 100, chunk_size 600
##
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
p    <- suppressWarnings(as.integer(args[1]))
if (is.na(p) || !(p %in% c(10L, 50L, 100L)))
  stop("Usage: Rscript mpi_simul_main.R <p> [chunk_size]   where p in {10, 50, 100}")

chunk_arg <- if (length(args) >= 2L && nzchar(args[2])) {
  v <- suppressWarnings(as.integer(args[2]))
  if (is.na(v) || v < 1L)
    stop("chunk_size must be a positive integer (got: '", args[2], "')")
  v
} else NA_integer_

source("./HazReg/routines.R")
library(mvtnorm)
library(survival)
library(glmnet)
library(tictoc)
library(Rmpi)
library(doMPI)
library(doRNG)
library(Rcpp)
library(RcppArmadillo)

source("./utils/helper_functions.R")

## ---- Simulation design -----------------------------------------------------
n.vec    <- c(250L, 500L, 1000L)
n.simset <- 250L
censrat  <- c("0.25default", "0.5default")
priorstr <- "LCM"  # or "Product"
filename <- paste0("mpi_simulation_raw_p", p, ".RData")

WALL_LIMIT_SECS <- 47L * 3600L

sigma <- 0.7 ^ abs(matrix(rep(seq_len(p), p), nrow = p, byrow = TRUE) -
                   seq_len(p))
hstr     <- c("PH", "AFT", "GH")
baseline <- c("LN", "PGW")

str.grid <- expand.grid(n.vec, p, hstr, baseline, censrat)
str.grid <- str.grid[!(str.grid[, 3] == "GH" & str.grid[, 4] == "PGW"), ]
row.names(str.grid) <- NULL
names(str.grid) <- c("n", "p", "hstr", "baseline")

true.models <- apply(str.grid, 1, function(x) {
  pv  <- as.integer(x[["p"]]); hv <- x[["hstr"]]
  idx <- as.integer(quantile(seq_len(pv), c(0.2, 0.4, 0.6, 0.8)))
  tmod <- rep(0L, pv)
  if      (hv == "PH")  tmod[idx] <- 2L
  else if (hv == "AFT") tmod[idx] <- 4L
  else if (hv == "GH")  tmod[idx] <- 3L
  else if (hv == "AH")  tmod[idx] <- 1L
  paste0(tmod, collapse = "")
})
str.grid <- cbind(str.grid, true.models)
names(str.grid) <- c("n", "p", "hstr", "baseline", "cens", "true.mod")

rm(list = ls()[!ls() %in% c("n.simset", "str.grid", "sigma", "filename",
                             "priorstr", "WALL_LIMIT_SECS", "p", "chunk_arg")])

mcmcpar <- list(n.iter = 20000, n.burnin = 10000, n.thin = 2,
                n.chain = 1, printevery = 10000)

sim_indices <- expand.grid(i = seq_len(nrow(str.grid)), j = seq_len(n.simset))

source("./utils/checkpoint.R")
CK_OBJECTS <- c("str.grid", "results_parallel", "n.simset")
save_ck <- function(extra = character(0)) {
  atomic_save(c(CK_OBJECTS, extra), file = filename, envir = globalenv())
}

## ---- Checkpoint ------------------------------------------------------------
ck <- load_checkpoint_results(filename, "results_parallel")
report_checkpoint(ck, filename)
results_parallel <- ck$results
rm(ck)

if (length(results_parallel)) {
  results_parallel <- Filter(
    function(r) is.list(r) && !is.null(r[["fullchain"]]),
    results_parallel)
  done_keys <- vapply(results_parallel,
                      function(r) paste0(r$i, "_", r$j), character(1L))
  if (is_rank0()) save_ck()
  cat(sprintf("Checkpoint loaded: %d / %d tasks successfully completed.\n",
              length(done_keys), nrow(sim_indices)))
} else {
  done_keys <- character(0L)
}

all_keys  <- paste0(sim_indices$i, "_", sim_indices$j)
todo_rows <- which(!(all_keys %in% done_keys))
cat(sprintf("Tasks remaining: %d / %d\n", length(todo_rows), nrow(sim_indices)))

chunk_size <- if (!is.na(chunk_arg)) chunk_arg else
              if (p == 10L) 1000L else if (p == 50L) 300L else 600L
cat(sprintf("chunk_size = %d (%s)\n", chunk_size,
            if (is.na(chunk_arg)) "built-in default" else "from command line"))

## ---- MPI cluster -----------------------------------------------------------
cat("Starting MPI cluster\n")
cl        <- startMPIcluster()
n_workers <- cl$workerCount
cat("Using", n_workers, "MPI workers\n")
registerDoMPI(cl)

t_start <- proc.time()["elapsed"]
tic(msg = "###########################\nTotal time elapsed")

## ---- Main chunk loop -------------------------------------------------------
chunk_starts <- if (length(todo_rows) > 0L) {
  seq(1L, length(todo_rows), by = chunk_size)
} else integer(0)

chunk_secs <- numeric(0)

for (cs in chunk_starts) {

  elapsed <- as.numeric(proc.time()["elapsed"] - t_start)
  if (elapsed > WALL_LIMIT_SECS) {
    cat(sprintf("Wall-time guard reached (%.1f h). Checkpoint saved; resubmit.\n",
                elapsed / 3600))
    break
  }
  if (length(chunk_secs)) {
    projected <- elapsed + 1.15 * max(chunk_secs)
    if (projected > WALL_LIMIT_SECS) {
      cat(sprintf(paste0("Stopping at %.1f h: the next chunk would be expected ",
                         "to finish at ~%.1f h, past the %.1f h guard.\n",
                         "Checkpoint saved; resubmit.\n"),
                  elapsed / 3600, projected / 3600, WALL_LIMIT_SECS / 3600))
      break
    }
  }
  t_chunk <- proc.time()["elapsed"]

  ce        <- min(cs + chunk_size - 1L, length(todo_rows))
  chunk_sim <- sim_indices[todo_rows[cs:ce], , drop = FALSE]

  chunk_res <- foreach(
    idx = seq_len(nrow(chunk_sim)),
    .packages = c("mvtnorm", "survival", "glmnet",
                  "tictoc", "RcppArmadillo", "Rcpp"),
    .options.RNG = 1234L + cs
  ) %dorng% {

    i <- chunk_sim$i[idx]
    j <- chunk_sim$j[idx]

    source("./ADS_sampler.R")
    source("./HazReg/routines.R")

    ## -- Data generation (seed = j ensures reproducibility) -----------------
    dat_full <- simdat(
      n = str.grid[i, 1], p = str.grid[i, 2], sigma = sigma,
      hstr = str.grid[i, 3], baseline = str.grid[i, 4],
      n.simset = 1, cens = str.grid[i, 5], seeds = j)
    dat      <- dat_full$data[[1]]
    cens_rat <- dat_full$cens.rat[1]

    n   <- nrow(dat$x)
    pxt <- ncol(dat$xt)
    px  <- ncol(dat$x)

    ## -- Initialisation ------------------------------------------------------
    t_init       <- Sys.time()
    initgam      <- init.gam(dat, history.keep = FALSE)
    init.runtime <- as.numeric(difftime(Sys.time(), t_init, units = "secs"))

    modelpar <- list(
      inits = list(nu = rnorm(1), theta0 = rnorm(1),
                   theta = rnorm(pxt, 0, 0.5), eta = rnorm(px, 0, 0.5),
                   g.theta = 1, g.eta = 1, g.aft = 1, g = 1,
                   gam = initgam, a_g = NULL,
                   nu.alpha = 0.01, nu.beta = 0.01),
      data = dat
    )

    ## -- Sampler -------------------------------------------------------------
    tryCatch({
      t_samp          <- Sys.time()
      chains          <- BVS.MCMC(modelpar, mcmcpar, verbose = FALSE,
                                  gprior.mat = priorstr, gamprior = "default")
      sampler_runtime <- as.numeric(difftime(Sys.time(), t_samp, units = "secs"))

      acc_idx   <- seq(mcmcpar$n.burnin + mcmcpar$n.thin,
                       mcmcpar$n.iter, mcmcpar$n.thin)
      gam_strs  <- apply(chains$chain[[1]], 1,
                         function(x) paste0(x, collapse = ""))
      dict      <- chains$model.prob
      pp_map    <- setNames(as.numeric(dict$postprob),
                            as.character(dict$gam))
      fullchain <- data.frame(
        gam       = gam_strs,
        acc.gam   = as.numeric(chains$acc[acc_idx, ]),
        post.prob = pp_map[gam_strs],
        stringsAsFactors = FALSE
      )

      list(i = i, j = j, fullchain = fullchain, cens_rat = cens_rat,
           runtime = c(init = init.runtime, sampler = sampler_runtime))
    }, error = function(e) {
      msg <- paste0(e$message, "\nn=", str.grid[i, 1], ", p=", str.grid[i, 2],
                    " ", str.grid[i, 3], "_", str.grid[i, 4])
      cat(msg, "\n")
      list(i = i, j = j, error = msg, fullchain = NULL)
    })
  }

  results_parallel <- c(results_parallel, chunk_res)
  save_ck()
  chunk_secs <- c(chunk_secs, as.numeric(proc.time()["elapsed"] - t_chunk))
  cat(sprintf("Chunk rows %d-%d done in %.2f h. Total saved: %d / %d tasks.\n",
              cs, ce, tail(chunk_secs, 1) / 3600,
              length(results_parallel), nrow(sim_indices)))
}

## ---- Finalise --------------------------------------------------------------
closeCluster(cl)
time.elaps <- toc()
time.elaps <- time.elaps$toc - time.elaps$tic

n_done <- length(Filter(function(r) is.list(r) && !is.null(r[["fullchain"]]),
                        results_parallel))
save_ck("time.elaps")
cat(sprintf("Finished: %d / %d tasks stored in %s.\n",
            n_done, nrow(sim_indices), filename))
if (n_done < nrow(sim_indices))
  cat(sprintf("%d task(s) still outstanding -- resubmit the jobscript.\n",
              nrow(sim_indices) - n_done))
mpi.quit()
