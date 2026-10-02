## ============================================================================
##
## Frequentist per-structure forward stepwise selection over the simulation grid
## of mpi_simul_main.R
##
## Method (utils/freq_stepwise_multi.R, step.gh.freq):
##   Stage 1  four independent forward stepwise searches -- AH, PH, AFT, GH --
##            each step adding exactly one parameter.
##   Stage 2  the four winners are scored on the same criterion and the global
##            minimiser is returned as the single selected model.
##
## Usage (via srun on Slurm, or mpirun):
##   srun Rscript mpi_freq_stepwise.R <p> [crit] [chunk_size]
##
##   srun Rscript mpi_freq_stepwise.R 10                  # p = 10,  BIC
##   srun Rscript mpi_freq_stepwise.R 50  BIC 300         # p = 50,  chunk 300
##
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
p    <- suppressWarnings(as.integer(args[1]))
if (is.na(p) || !(p %in% c(10L, 50L, 100L)))
  stop("Usage: Rscript mpi_freq_stepwise.R <p> [crit]   where p in {10, 50, 100}")

CRIT <- if (length(args) >= 2L && nzchar(args[2])) toupper(args[2]) else "BIC"
if (!(CRIT %in% c("BIC", "AIC", "LASSO")))
  stop("crit must be one of BIC, AIC, LASSO (got: ", CRIT, ")")

chunk_arg <- if (length(args) >= 3L && nzchar(args[3])) {
  v <- suppressWarnings(as.integer(args[3]))
  if (is.na(v) || v < 1L)
    stop("chunk_size must be a positive integer (got: '", args[3], "')")
  v
} else NA_integer_

library(mvtnorm)
library(survival)
library(glmnet)
library(tictoc)
library(Rmpi)
library(doMPI)
library(doRNG)

source("./HazReg/routines.R")           # GHMLE and the baseline hazards
source("./utils/helper_functions.R")    # simdat, gam2haz, metrics
source("./utils/freq_stepwise_multi.R") # step.gh.freq

## ---- Simulation design -----------------------------------------------------
n.vec    <- c(250L, 500L, 1000L)
n.simset <- 250L
censrat  <- c("0.25default", "0.5default")
filename <- sprintf("mpi_freq_%s_p%d.RData", tolower(CRIT), p)

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
FIT_DIST   <- "LogNormal"
MAX_STEPS  <- min(2L * p, 40L)
MAX_ACTIVE <- min(p, 20L)
GH_STRICT  <- FALSE
LAMBDA     <- "cox"
MAXIT      <- 5000L

rm(list = ls()[!ls() %in% c("n.simset", "str.grid", "sigma", "filename",
                            "WALL_LIMIT_SECS", "p", "CRIT", "FIT_DIST",
                            "MAX_STEPS", "MAX_ACTIVE", "GH_STRICT",
                            "LAMBDA", "MAXIT", "chunk_arg")])

sim_indices <- expand.grid(i = seq_len(nrow(str.grid)), j = seq_len(n.simset))


source("./utils/checkpoint.R")
CK_OBJECTS <- c("str.grid", "results_parallel", "n.simset",
                "CRIT", "FIT_DIST", "MAX_STEPS", "MAX_ACTIVE",
                "GH_STRICT", "LAMBDA")
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
    function(r) is.list(r) && !is.null(r[["sel"]]),
    results_parallel)
  done_keys <- vapply(results_parallel,
                      function(r) paste0(r$i, "_", r$j), character(1L))

  if (is_rank0()) save_ck()
  cat(sprintf("Checkpoint loaded: %d / %d tasks successfully completed.\n",
              length(done_keys), nrow(sim_indices)))
} else {
  done_keys        <- character(0L)
}

all_keys  <- paste0(sim_indices$i, "_", sim_indices$j)
todo_rows <- which(!(all_keys %in% done_keys))
cat(sprintf("Criterion: %s | p = %d | output: %s\n", CRIT, p, filename))
cat(sprintf("Tasks remaining: %d / %d\n", length(todo_rows), nrow(sim_indices)))

chunk_size <- if (!is.na(chunk_arg)) chunk_arg else
              if (p == 10L) 500L else if (p == 50L) 100L else 36L
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
    .packages = c("mvtnorm", "survival", "glmnet", "tictoc"),
    .export   = c("chunk_sim", "str.grid", "sigma", "CRIT", "FIT_DIST",
                  "MAX_STEPS", "MAX_ACTIVE", "GH_STRICT", "LAMBDA", "MAXIT"),
    .options.RNG = 1234L + cs
  ) %dorng% {

    i <- chunk_sim$i[idx]
    j <- chunk_sim$j[idx]

    source("./HazReg/routines.R")
    source("./utils/helper_functions.R")
    source("./utils/freq_stepwise_multi.R")

    ## -- Data generation (seed = j ensures reproducibility, and gives the
    ##    SAME data set as mpi_simul_main.R for the same (i, j)) -------------
    dat_full <- simdat(
      n = str.grid[i, 1], p = str.grid[i, 2], sigma = sigma,
      hstr = str.grid[i, 3], baseline = str.grid[i, 4],
      n.simset = 1, cens = str.grid[i, 5], seeds = j)
    dat      <- dat_full$data[[1]]
    cens_rat <- dat_full$cens.rat[1]

    true.mod <- as.character(str.grid[i, 6])
    true.hz  <- gam2haz(true.mod)

    ## -- Selection ------------------------------------------------------------
    tryCatch({
      t_sel <- Sys.time()
      fsel  <- step.gh.freq(dat,
                            crit        = CRIT,
                            lambda      = LAMBDA,
                            dist        = FIT_DIST,
                            maxit       = MAXIT,
                            max.steps   = MAX_STEPS,
                            max.active  = MAX_ACTIVE,
                            gh.strict   = GH_STRICT,
                            verbose     = FALSE)
      sel_runtime <- as.numeric(difftime(Sys.time(), t_sel, units = "secs"))

      b  <- fsel$best
      mm <- metrics(b$gam, true.mod)

      hit.cap <- any(vapply(fsel$arms, function(a) a$steps >= MAX_STEPS,
                            logical(1L))) ||
                 any(fsel$table$n.active >= MAX_ACTIVE, na.rm = TRUE)

      list(
        i        = i,
        j        = j,
        cens_rat = cens_rat,
        ## Selected model. 
        sel = list(gam.str   = b$gam.str,
                   structure = b$structure,   # which arm won
                   hstr      = b$hstr,        # gam2haz of the selected model
                   score     = b$score,
                   loglik    = b$loglik,
                   k         = b$k,
                   AIC       = b$AIC,
                   BIC       = b$BIC),
        ## Per-arm winners with loglik / AIC / BIC side by side
        arm.table = fsel$table,
        coef      = b$coef,
        metrics   = c(Sens = mm$Sens, Spec = mm$Spec, FDR = mm$FDR,
                      phi = mm$phi,
                      haz_correct = as.numeric(b$hstr == true.hz)),
        lambda    = fsel$lambda,
        n.fits    = fsel$n.fits,
        hit.cap   = hit.cap,
        runtime   = c(selection = sel_runtime),
        error     = NA_character_
      )
    }, error = function(e) {
      msg <- paste0(conditionMessage(e), "\nn=", str.grid[i, 1],
                    ", p=", str.grid[i, 2], " ", str.grid[i, 3], "_",
                    str.grid[i, 4])
      cat(msg, "\n")
      list(i = i, j = j, error = msg, sel = NULL)
    })
  }

  results_parallel <- c(results_parallel, chunk_res)
  save_ck()
  chunk_secs <- c(chunk_secs, as.numeric(proc.time()["elapsed"] - t_chunk))

  ok  <- vapply(chunk_res, function(r) !is.null(r[["sel"]]), logical(1L))
  rts <- vapply(chunk_res[ok], function(r) r$runtime[["selection"]], numeric(1))
  cap <- sum(vapply(chunk_res[ok], function(r) isTRUE(r$hit.cap), logical(1L)))
  cat(sprintf(paste0("Chunk rows %d-%d done (%d ok, %d failed, %d hit a guard). ",
                     "Mean %.1f s/task. Total saved: %d / %d tasks.\n"),
              cs, ce, sum(ok), sum(!ok), cap,
              if (length(rts)) mean(rts) else NA_real_,
              length(results_parallel), nrow(sim_indices)))
}

## ---- Finalise --------------------------------------------------------------
closeCluster(cl)
time.elaps <- toc()
time.elaps <- time.elaps$toc - time.elaps$tic

n_done <- length(Filter(function(r) is.list(r) && !is.null(r[["sel"]]),
                        results_parallel))
save_ck("time.elaps")
cat(sprintf("Finished: %d / %d tasks stored in %s.\n",
            n_done, nrow(sim_indices), filename))
if (n_done < nrow(sim_indices))
  cat(sprintf("%d task(s) still outstanding -- resubmit the jobscript.\n",
              nrow(sim_indices) - n_done))
mpi.quit()
