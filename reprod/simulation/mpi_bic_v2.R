## ============================================================================
## Usage (via srun on Slurm, or mpirun):
##   srun Rscript mpi_bic_v2.R <p> [chunk_size]
##
##   srun Rscript mpi_bic_v2.R 10                 # p = 10, default chunk
##   srun Rscript mpi_bic_v2.R 50   400           # p = 50, chunk_size 400
## ============================================================================

print("current wd:")
print(getwd())
source("./HazReg/routines.R")

args <- commandArgs(trailingOnly = TRUE)
p    <- suppressWarnings(as.integer(args[1]))
if (is.na(p) || !(p %in% c(10L, 50L)))
  stop("Usage: Rscript mpi_bic_v2.R <p> [chunk_size]   where p in {10, 50}")

chunk_arg <- if (length(args) >= 2L && nzchar(args[2])) {
  v <- suppressWarnings(as.integer(args[2]))
  if (is.na(v) || v < 1L)
    stop("chunk_size must be a positive integer (got: '", args[2], "')")
  v
} else NA_integer_

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
n.vec    <- c(250L, 1000L)
n.simset <- 250L
censrat  <- c("0.25default", "0.5default")
priorstr <- "LCM"
filename <- paste0("mpi_bic_p", p, ".RData")

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

## ---- Method configuration --------------------------------------------------
##  RUN_LCM     : include the proposed LCM method?
##  BIC_METHODS : BIC MCMC variants (same keys as original mpi_bic.R)
##  freq_gh_bic : always run; is a point-estimate comparator, not MCMC
RUN_LCM <- FALSE

BIC_METHODS <- list(
  bic_flat_n = list(prior = "flat", n_type = "n",  label = "BIC-flat (n)"),
  bic_bb_n   = list(prior = "bb",   n_type = "n",  label = "BIC-BB (n)")
)

rm(list = ls()[!ls() %in% c("n.simset", "str.grid", "sigma", "filename",
                             "priorstr", "WALL_LIMIT_SECS", "p", "chunk_arg",
                             "RUN_LCM", "BIC_METHODS")])

mcmcpar <- list(n.iter = 20000, n.burnin = 10000, n.thin = 2,
                n.chain = 1, printevery = 10000)

## Number of saved iterations (used to size the degenerate freq chain)
N_SAVE <- floor((mcmcpar$n.iter - mcmcpar$n.burnin) / mcmcpar$n.thin)

sim_indices <- expand.grid(i = seq_len(nrow(str.grid)), j = seq_len(n.simset))

source("./utils/checkpoint.R")
CK_OBJECTS <- c("str.grid", "results_parallel", "n.simset",
                "RUN_LCM", "BIC_METHODS")
save_ck <- function(extra = character(0)) {
  atomic_save(c(CK_OBJECTS, extra), file = filename, envir = globalenv())
}

## ---- Checkpoint ------------------------------------------------------------
ck <- load_checkpoint_results(filename, "results_parallel")
report_checkpoint(ck, filename)
results_parallel <- ck$results
rm(ck)

if (length(results_parallel)) {
  n_loaded <- length(results_parallel)
  is_current <- vapply(results_parallel,
                       function(r) is.list(r) && "freq_info" %in% names(r),
                       logical(1L))
  n_stale <- sum(!is_current)
  if (n_stale > 0L && is_rank0()) {
    backup <- sub("\\.RData$", "_prev.RData", filename)
    if (!file.exists(backup)) {
      file.copy(filename, backup)
      cat(sprintf("Found %d result(s) predating the rebuilt frequentist\n", n_stale))
      cat(sprintf("  comparator; original backed up to %s.\n", backup))
    } else {
      cat(sprintf("Found %d stale result(s) (backup %s already exists).\n",
                  n_stale, backup))
    }
  }

  keep <- vapply(seq_along(results_parallel), function(k) {
    r <- results_parallel[[k]]
    if (!is_current[k]) return(FALSE)
    if (is.null(r[["i"]]) || is.null(r[["j"]])) return(FALSE)
    has_lcm <- !is.null(r[["fullchain"]])
    has_bic <- !is.null(r[["fullchain_bic"]]) &&
               any(!vapply(r[["fullchain_bic"]], is.null, logical(1L)))
    has_lcm || has_bic
  }, logical(1L))
  results_parallel <- results_parallel[keep]

  done_keys <- vapply(results_parallel,
                      function(r) paste0(r$i, "_", r$j), character(1L))

  if (is_rank0()) save_ck()
  cat(sprintf("Checkpoint loaded: %d of %d stored result(s) usable (%d / %d tasks).\n",
              length(done_keys), n_loaded, length(done_keys), nrow(sim_indices)))
} else {
  done_keys        <- character(0L)
}

all_keys  <- paste0(sim_indices$i, "_", sim_indices$j)
todo_rows <- which(!(all_keys %in% done_keys))
cat(sprintf("Tasks remaining: %d / %d\n", length(todo_rows), nrow(sim_indices)))

chunk_size <- if (!is.na(chunk_arg)) chunk_arg else
              if (p == 10L) 500L else 400L
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
    cat(sprintf("Wall-time guard reached (%.1f h). Saving checkpoint.\n",
                elapsed / 3600))
    break
  }

  if (length(chunk_secs)) {
    projected <- elapsed + 1.15 * max(chunk_secs)
    if (projected > WALL_LIMIT_SECS) {
      cat(sprintf(paste0("Stopping at %.1f h: the next chunk would be expected ",
                         "to finish at ~%.1f h, past the %.1f h guard.
",
                         "Checkpoint saved; resubmit.
"),
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
    source("./utils/EB_gprior_fixed.R")
    source("./utils/bic_functions.R")
    source("./utils/freq_stepwise_bic.R")

    ## -- Data generation -----------------------------------------------------
    dat_full <- simdat(
      n = str.grid[i, 1], p = str.grid[i, 2], sigma = sigma,
      hstr = str.grid[i, 3], baseline = str.grid[i, 4],
      n.simset = 1, cens = str.grid[i, 5], seeds = j)
    dat      <- dat_full$data[[1]]
    cens_rat <- dat_full$cens.rat[1]

    n   <- nrow(dat$x)
    pxt <- ncol(dat$xt)
    px  <- ncol(dat$x)

    ## -- Shared initialisation -----------------------------------------------
    t_init <- Sys.time()
    initgam <- init.gam(dat, history.keep = FALSE)
    init.runtime <- as.numeric(difftime(Sys.time(), t_init, units = "secs"))

    modelpar <- list(
      inits = list(nu = rnorm(1), theta0 = rnorm(1),
                   theta = rnorm(pxt, 0, 0.5), eta = rnorm(px, 0, 0.5),
                   g.theta = 1, g.eta = 1, g.aft = 1, g = 1,
                   gam = initgam, a_g = NULL,
                   nu.alpha = 0.01, nu.beta = 0.01),
      data = dat
    )

    reset_mle_cache()

    ## ================================================================== ##
    ##  LCM (proposed) -- only if RUN_LCM = TRUE                          ##
    ## ================================================================== ##
    lcm_fullchain <- NULL
    lcm_runtime   <- NA_real_
    lcm_error     <- NA_character_
    if (isTRUE(RUN_LCM)) {
      lcm_out <- tryCatch({
        t0 <- Sys.time()
        chains_lcm <- BVS.MCMC(modelpar, mcmcpar, verbose = FALSE,
                               gprior.mat = priorstr, gamprior = "default")
        t1 <- Sys.time()
        list(fullchain = build_fullchain(chains_lcm, mcmcpar),
             runtime   = as.numeric(difftime(t1, t0, units = "secs")),
             err       = NA_character_)
      }, error = function(e) {
        cat(paste0("LCM error: ", e$message,
                   "\n  n=", str.grid[i, 1], ", p=", str.grid[i, 2], " ",
                   str.grid[i, 3], "_", str.grid[i, 4], "\n"))
        list(fullchain = NULL, runtime = NA_real_,
             err = paste0(e$message, " (LCM)"))
      })
      lcm_fullchain <- lcm_out$fullchain
      lcm_runtime   <- lcm_out$runtime
      lcm_error     <- lcm_out$err
    }

    ## ================================================================== ##
    ##  BIC MCMC variants                                                  ##
    ## ================================================================== ##
    fullchain_bic <- vector("list", length(BIC_METHODS))
    names(fullchain_bic) <- names(BIC_METHODS)
    bic_runtimes <- setNames(rep(NA_real_,      length(BIC_METHODS)),
                             names(BIC_METHODS))
    bic_errors   <- setNames(rep(NA_character_, length(BIC_METHODS)),
                             names(BIC_METHODS))

    for (m_name in names(BIC_METHODS)) {
      meth <- BIC_METHODS[[m_name]]
      target.fn <- make.targetdist.bic(prior_type = meth$prior,
                                       n_type     = meth$n_type)
      gp_type <- if (meth$prior == "bb") "default" else "unif"
      one_bic <- tryCatch({
        t0 <- Sys.time()
        chains_b <- BVS.MCMC.bic(modelpar, mcmcpar,
                                  target.fn     = target.fn,
                                  gamprior_type = gp_type,
                                  verbose       = FALSE)
        t1 <- Sys.time()
        list(fullchain = build_fullchain(chains_b, mcmcpar),
             runtime   = as.numeric(difftime(t1, t0, units = "secs")),
             err       = NA_character_)
      }, error = function(e) {
        cat(paste0("BIC error (", m_name, "): ", e$message,
                   "\n  n=", str.grid[i, 1], ", p=", str.grid[i, 2], " ",
                   str.grid[i, 3], "_", str.grid[i, 4], "\n"))
        list(fullchain = NULL, runtime = NA_real_,
             err = paste0(e$message, " (", m_name, ")"))
      })
      fullchain_bic[[m_name]] <- one_bic$fullchain
      bic_runtimes[m_name]    <- one_bic$runtime
      bic_errors[m_name]      <- one_bic$err
    }

    freq_t0  <- Sys.time()
    freq_out <- tryCatch({
      fsel <- step.bic.gh(dat, dist = "LogNormal", maxit = 5000)
      list(
        fullchain = data.frame(
          gam       = rep(fsel$gam.str, N_SAVE),
          acc.gam   = rep(1.0, N_SAVE),
          post.prob = rep(NA_real_, N_SAVE),
          stringsAsFactors = FALSE),
        info = list(gam = fsel$gam.str, hstr = fsel$hstr, bic = fsel$bic,
                    bic_armA = fsel$bic.armA, bic_armB = fsel$bic.armB,
                    n_fits = fsel$n.fits, steps = fsel$steps),
        err = NA_character_)
    }, error = function(e) {
      cat(paste0("freq_gh_bic error: ", conditionMessage(e),
                 "\n  n=", str.grid[i, 1], ", p=", str.grid[i, 2], " ",
                 str.grid[i, 3], "_", str.grid[i, 4], "\n"))
      list(fullchain = NULL, info = NULL, err = conditionMessage(e))
    })
    freq_fullchain <- freq_out$fullchain
    freq_info      <- freq_out$info
    freq_err       <- freq_out$err
    freq_runtime <- as.numeric(difftime(Sys.time(), freq_t0, units = "secs"))

    fullchain_bic[["freq_gh_bic"]] <- freq_fullchain
    bic_runtimes["freq_gh_bic"]    <- freq_runtime
    bic_errors["freq_gh_bic"]      <- freq_err

    ## -- Return --------------------------------------------------------------
    list(
      i             = i,
      j             = j,
      cens_rat      = cens_rat,
      fullchain     = lcm_fullchain,
      fullchain_bic = fullchain_bic,
      freq_info     = freq_info,   # selected model, BIC per arm, fit count
      runtime = c(init = init.runtime, lcm = lcm_runtime, bic_runtimes),
      errors  = c(lcm = lcm_error, bic_errors)
    )
  }

  results_parallel <- c(results_parallel, chunk_res)
  save_ck()
  chunk_secs <- c(chunk_secs, as.numeric(proc.time()["elapsed"] - t_chunk))
  cat(sprintf("Chunk %d-%d done. Total saved: %d / %d\n",
              cs, ce, length(results_parallel), nrow(sim_indices)))
}

## ---- Finalise --------------------------------------------------------------
closeCluster(cl)
time.elaps <- toc()
time.elaps <- time.elaps$toc - time.elaps$tic

n_done <- length(results_parallel)
save_ck("time.elaps")
cat(sprintf("Finished: %d / %d tasks stored in %s.\n",
            n_done, nrow(sim_indices), filename))
if (n_done < nrow(sim_indices))
  cat(sprintf("%d task(s) still outstanding -- resubmit the jobscript.\n",
              nrow(sim_indices) - n_done))
mpi.quit()
