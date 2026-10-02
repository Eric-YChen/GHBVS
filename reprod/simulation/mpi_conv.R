## ============================================================================
## MCMC convergence study.
##
## Runs 8 chains per (scenario, replicate) from deliberately different starting
## models, stores the COMPLETE per-iteration trace of each chain (including
## burn-in and iterations that would normally be thinned out), and reports
## per-variable Gelman-Rubin potential scale reduction factors (PSRFs).
##
## Usage (via srun on Slurm, or mpirun):
##   srun Rscript mpi_conv.R <p> [chunk_size]
##
##   srun Rscript mpi_conv.R 10                   # p = 10,  default chunk
##   srun Rscript mpi_conv.R 50   150             # p = 50,  chunk_size 150
##   srun Rscript mpi_conv.R 100                  # p = 100, default chunk
##
## ---------------------------------------------------------------------------
## Chain initialisation (8 chains)
##   chains 1-2 : init.gam()  -- the paper's default (Cox-LASSO + BIC stepwise).
##                Identical starting gamma; they differ only in RNG stream,
##                which is the relevant comparison because the MCMC state is
##                gamma alone (the continuous parameters are profiled out by
##                the Laplace approximation inside targetdist).
##   chains 3-6 : random starting model at the LASSO-implied sparsity level,
##                i.e. the same model size as init.gam() but a randomly chosen
##                active set and a randomly chosen hazard structure.
##   chain    7 : under-dispersed -- a single active covariate (close to null).
##   chain    8 : over-dispersed  -- twice the active count of init.gam().
##
## Every random start is drawn to be a VALID hazard structure (AH / PH / AFT /
## GH as recognised by gam2haz); sampling gamma entries independently from
## {1,2,3,4} can yield "Invalid model" states that the ADS proposal is not
## designed to start from.
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
p    <- suppressWarnings(as.integer(args[1]))
if (is.na(p) || !(p %in% c(10L, 50L, 100L)))
  stop("Usage: Rscript mpi_conv.R <p> [chunk_size]   where p in {10, 50, 100}")

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

## ---- Simulation design (same grid as the main simulation) ------------------
n.vec    <- c(250L, 500L, 1000L)
n.simset <- 10L
censrat  <- c("0.25default", "0.5default")
priorstr <- "LCM"
filename  <- paste0("mpi_conv_p", p, ".RData")
trace_dir <- paste0("traces_p", p)

SAVE_TRACES <- TRUE

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

mcmcpar <- list(n.iter = 50000, n.burnin = 25000, n.thin = 2,
                n.chain = 1,          
                printevery = 50000)

mcmcpar_full <- modifyList(mcmcpar,
                           list(n.burn = 0L, n.burnin = 0L, n.thin = 1L))

## ---- Chain / initialisation configuration ---------------------------------
INIT_TYPES <- c("init.gam", "init.gam",
                "rand_sparse", "rand_sparse", "rand_sparse", "rand_sparse",
                "underdispersed", "overdispersed")
N_CHAINS   <- length(INIT_TYPES)

PRACTICAL_CHAINS <- which(INIT_TYPES %in% c("init.gam", "rand_sparse"))

rm(list = ls()[!ls() %in% c("n.simset", "str.grid", "sigma", "filename",
                             "trace_dir", "SAVE_TRACES", "priorstr",
                             "WALL_LIMIT_SECS", "p", "mcmcpar", "mcmcpar_full",
                             "chunk_arg",
                             "INIT_TYPES", "N_CHAINS", "PRACTICAL_CHAINS")])

sim_indices <- expand.grid(i = seq_len(nrow(str.grid)), j = seq_len(n.simset))

source("./utils/checkpoint.R") 
CK_OBJECTS <- c("str.grid", "results_parallel", "n.simset",
                "mcmcpar", "mcmcpar_full", "N_CHAINS", "INIT_TYPES",
                "PRACTICAL_CHAINS", "trace_dir")
save_ck <- function(extra = character(0)) {
  atomic_save(c(CK_OBJECTS, extra), file = filename, envir = globalenv())
}

## ---------------------------------------------------------------------------
##  Random starting model at a prescribed size
rand_valid_gam <- function(px, m) {
  m   <- max(1L, min(as.integer(m), px))
  idx <- sample.int(px, m)
  gam <- rep(0L, px)
  strc <- sample(c("AH", "PH", "AFT", "GH"), 1L)
  if (strc == "AH") {
    gam[idx] <- 1L
  } else if (strc == "PH") {
    gam[idx] <- 2L
  } else if (strc == "AFT") {
    gam[idx] <- 4L
  } else {
    gam[idx] <- sample(1:3, m, replace = TRUE)
    if (all(gam[idx] == 1L) || all(gam[idx] == 2L)) gam[idx[1]] <- 3L
  }
  gam
}

## ---------------------------------------------------------------------------
##  Per-variable Gelman-Rubin R-hat on binary inclusion indicators.
compute_rhat <- function(chain_list) {
  valid <- Filter(Negate(is.null), lapply(chain_list, function(m) {
    if (is.null(m)) return(NULL)
    m <- m[complete.cases(m), , drop = FALSE]
    if (nrow(m) < 20L) return(NULL)
    m
  }))
  if (length(valid) < 2L) return(NULL)

  n     <- min(vapply(valid, nrow, integer(1L)))
  valid <- lapply(valid, function(x) x[seq_len(n), , drop = FALSE])
  m     <- length(valid)
  pv    <- ncol(valid[[1L]])

  vapply(seq_len(pv), function(k) {
    chains_k    <- vapply(valid, function(x) x[, k], numeric(n))
    chain_means <- colMeans(chains_k)
    grand_mean  <- mean(chain_means)
    B <- n / (m - 1) * sum((chain_means - grand_mean)^2)
    W <- mean(apply(chains_k, 2L, var))
    if (W < .Machine$double.eps) return(NA_real_)
    sqrt(((n - 1) / n * W + B / n) / W)
  }, numeric(1L))
}

is_current_schema <- function(r) {
  is.list(r) &&
    !is.null(r[["init_models"]]) &&
    length(r[["init_models"]]) == N_CHAINS &&
    "trace_file" %in% names(r)
}

ck <- load_checkpoint_results(filename, "results_parallel")
report_checkpoint(ck, filename)
results_parallel <- ck$results
rm(ck)

if (length(results_parallel)) {
  n_loaded  <- length(results_parallel)
  cur       <- vapply(results_parallel, is_current_schema, logical(1L))
  n_stale   <- sum(!cur)

  if (n_stale > 0L && is_rank0()) {
    backup <- sub("\\.RData$", "_prev.RData", filename)
    if (!file.exists(backup)) {
      file.copy(filename, backup)
      cat(sprintf("Found %d result(s) from an earlier design; original file\n",
                  n_stale))
      cat(sprintf("  backed up to %s before being replaced.\n", backup))
    } else {
      cat(sprintf("Found %d result(s) from an earlier design (backup %s\n",
                  n_stale, backup))
      cat("  already exists; leaving it untouched).\n")
    }
  }

  results_parallel <- Filter(
    function(r) is_current_schema(r) && isTRUE(r[["n_chains_ok"]] >= 2L),
    results_parallel)
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

if (SAVE_TRACES) dir.create(trace_dir, showWarnings = FALSE, recursive = TRUE)

chunk_size <- if (!is.na(chunk_arg)) chunk_arg else
              if (p == 10L) 20L else if (p == 50L) 150L else 2L
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
    .export = c("rand_valid_gam", "compute_rhat",
                "str.grid", "sigma", "mcmcpar", "mcmcpar_full", "priorstr",
                "chunk_sim", "INIT_TYPES", "N_CHAINS", "PRACTICAL_CHAINS",
                "SAVE_TRACES", "trace_dir"),
    .options.RNG = 1234L + cs
  ) %dorng% {

    i <- chunk_sim$i[idx]
    j <- chunk_sim$j[idx]

    source("./ADS_sampler.R")
    source("./HazReg/routines.R")
    source("./utils/EB_gprior_fixed.R")

    ## -- Data generation -----------------------------------------------------
    dat_full <- simdat(
      n = str.grid[i, 1], p = str.grid[i, 2], sigma = sigma,
      hstr = str.grid[i, 3], baseline = str.grid[i, 4],
      n.simset = 1, cens = str.grid[i, 5], seeds = j)
    dat      <- dat_full$data[[1]]
    cens_rat <- dat_full$cens.rat[1]

    pxt <- ncol(dat$xt)
    px  <- ncol(dat$x)

    ## -- init.gam(): Cox-LASSO screening + BIC forward stepwise --------------
    t_init  <- Sys.time()
    initgam <- tryCatch(init.gam(dat, history.keep = FALSE),
                        error = function(e) {
                          tmp <- rep(0L, px); tmp[1L] <- 2L; tmp
                        })
    init.runtime <- as.numeric(difftime(Sys.time(), t_init, units = "secs"))

    ## LASSO-implied sparsity level, used for chains 3-6 and (doubled) chain 8
    m_init <- max(1L, sum(initgam != 0))

    ## -- Build the 8 starting models -----------------------------------------
    init_gams <- vector("list", N_CHAINS)
    for (ch in seq_len(N_CHAINS)) {
      init_gams[[ch]] <- switch(
        INIT_TYPES[ch],
        "init.gam"       = as.integer(initgam),
        "rand_sparse"    = rand_valid_gam(px, m_init),
        "underdispersed" = rand_valid_gam(px, 1L),
        "overdispersed"  = rand_valid_gam(px, min(2L * m_init, px))
      )
    }
    init_models <- vapply(init_gams,
                          function(g) paste0(g, collapse = ""), character(1L))

    ## -- Run the chains ------------------------------------------------------
    traces      <- vector("list", N_CHAINS)
    bin_chains  <- vector("list", N_CHAINS)   # post-burn-in thinned, binary
    pip_mat     <- matrix(NA_real_, nrow = N_CHAINS, ncol = px)
    acc_rate    <- rep(NA_real_, N_CHAINS)
    chain_err   <- rep(NA_character_, N_CHAINS)
    n_chains_ok <- 0L

    keep_idx <- seq(mcmcpar$n.burnin + mcmcpar$n.thin,
                    mcmcpar$n.iter, mcmcpar$n.thin)

    for (ch in seq_len(N_CHAINS)) {
      modelpar_ch <- list(
        inits = list(nu = rnorm(1), theta0 = rnorm(1),
                     theta = rnorm(pxt, 0, 0.5), eta = rnorm(px, 0, 0.5),
                     g.theta = 1, g.eta = 1, g.aft = 1, g = 1,
                     gam = init_gams[[ch]], a_g = NULL,
                     nu.alpha = 0.01, nu.beta = 0.01),
        data = dat
      )

      out <- tryCatch(
        BVS.MCMC(modelpar_ch, mcmcpar_full, verbose = FALSE,
                 gprior.mat = priorstr, gamprior = "default"),
        error = function(e) list(err = conditionMessage(e))
      )

      if (!is.null(out$err)) {
        chain_err[ch] <- out$err
        next
      }

      n_chains_ok <- n_chains_ok + 1L

      gam_all <- apply(out$chain[[1]], 1,
                       function(x) paste0(x, collapse = ""))
      pp_map  <- setNames(as.numeric(out$model.prob$postprob),
                          as.character(out$model.prob$gam))
      traces[[ch]] <- data.frame(gam       = gam_all,
                                 post.prob = unname(pp_map[gam_all]),
                                 stringsAsFactors = FALSE)
      acc_rate[ch] <- mean(as.numeric(out$acc[, 1]), na.rm = TRUE)
      ## Drop the n.iter x p numeric data.frame before starting the next chain
      rm(out)
      ## Binary inclusion indicators on the retained subset
      gam_keep <- traces[[ch]]$gam[keep_idx]
      bin_mat  <- do.call(rbind,
                          lapply(strsplit(gam_keep, ""),
                                 function(s) as.integer(s) != 0L)) * 1.0
      bin_chains[[ch]] <- bin_mat
      pip_mat[ch, ]    <- colMeans(bin_mat)
    }

    ## -- Write the full traces ----------------------------------------------
    trace_file <- NA_character_
    if (SAVE_TRACES && n_chains_ok > 0L) {
      trace_file <- file.path(trace_dir,
                              sprintf("trace_i%03d_j%03d.rds", i, j))
      saveRDS(list(i = i, j = j,
                   scenario    = str.grid[i, ],
                   mcmcpar     = mcmcpar,       # nominal burn-in / thinning
                   n_iter      = mcmcpar$n.iter,
                   init_types  = INIT_TYPES,
                   init_models = init_models,
                   traces      = traces), 
              file = trace_file)
    }

    ## -- Gelman-Rubin --------------------------------------------------------
    rhat_all  <- compute_rhat(bin_chains)
    rhat_prac <- compute_rhat(bin_chains[PRACTICAL_CHAINS])

    safe_mean <- function(v) if (is.null(v)) NA_real_ else mean(v, na.rm = TRUE)
    safe_max  <- function(v) {
      if (is.null(v) || all(is.na(v))) return(NA_real_)
      max(v, na.rm = TRUE)
    }

    ## -- PIP agreement across initialisation groups --------------------------
    grp <- unique(INIT_TYPES)
    pip_by_init <- t(vapply(grp, function(g) {
      rows <- which(INIT_TYPES == g)
      sub  <- pip_mat[rows, , drop = FALSE]
      if (all(is.na(sub))) rep(NA_real_, px) else colMeans(sub, na.rm = TRUE)
    }, numeric(px)))
    rownames(pip_by_init) <- grp

    max_pip_gap <- if (nrow(pip_by_init) < 2L) NA_real_ else {
      rng <- apply(pip_by_init, 2L, function(v) {
        v <- v[is.finite(v)]
        if (length(v) < 2L) NA_real_ else max(v) - min(v)
      })
      if (all(is.na(rng))) NA_real_ else max(rng, na.rm = TRUE)
    }

    list(i = i, j = j, cens_rat = cens_rat,
         mean_psrf            = safe_mean(rhat_all),
         max_psrf             = safe_max(rhat_all),
         mean_psrf_practical  = safe_mean(rhat_prac),
         max_psrf_practical   = safe_max(rhat_prac),
         psrf_per_var         = rhat_all,
         pip_by_init          = pip_by_init,
         max_pip_gap          = max_pip_gap,
         acc_rate             = acc_rate,
         init_models          = init_models,
         m_init               = m_init,
         n_chains_ok          = n_chains_ok,
         chain_err            = chain_err,
         init_runtime         = init.runtime,
         trace_file           = trace_file)
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
