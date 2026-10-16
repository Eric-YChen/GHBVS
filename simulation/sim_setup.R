# Environment variables: SIM_LOCAL=1 (run without MPI), SIM_QUICK=1 (short test run),
# SIM_WALL_HOURS (wall-time limit in hours, default 47).
SIM_LOCAL <- nzchar(Sys.getenv("SIM_LOCAL"))
SIM_QUICK <- nzchar(Sys.getenv("SIM_QUICK"))

sim_paths <- function() {
  sim <- normalizePath(".", winslash = "/")
  list(sim = sim, core = normalizePath(file.path(sim, "..", "core"), winslash = "/"),
       tasks = file.path(sim, "sim_tasks.R"), results = file.path(sim, "results"))
}

source(sim_paths()$tasks)

make_grid <- function(p, n.vec = c(250L, 500L, 1000L), hstr = c("PH", "AFT", "GH"),
                      baseline = c("LN", "PGW"), cens = c("0.25default", "0.5default")) {
  g <- expand.grid(n = n.vec, p = p, hstr = hstr, baseline = baseline, cens = cens,
                   stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  g <- g[!(g$hstr == "GH" & g$baseline == "PGW"), ]
  g$true.mod <- vapply(seq_len(nrow(g)), function(r) true_model(g$p[r], g$hstr[r]), character(1))
  rownames(g) <- NULL
  g
}

atomic_save <- function(names, file, envir) {
  tmp <- paste0(file, ".tmp-", Sys.info()[["nodename"]], "-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  save(list = names, file = tmp, envir = envir)
  if (!file.rename(tmp, file)) {
    if (!file.copy(tmp, file, overwrite = TRUE)) stop("could not write ", file)
  }
  invisible(file)
}

load_checkpoint <- function(file) {
  if (!file.exists(file)) { cat("No checkpoint found; starting from scratch.\n"); return(list()) }
  ck <- new.env(parent = emptyenv())
  err <- tryCatch({ load(file, envir = ck); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(err)) {
    quar <- paste0(file, ".corrupt-", format(Sys.time(), "%Y%m%d-%H%M%S"))
    file.rename(file, quar)
    cat("*** Checkpoint unreadable (", err, "); moved to ", quar, "; starting from scratch.\n", sep = "")
    return(list())
  }
  if (exists("results_parallel", envir = ck, inherits = FALSE)) get("results_parallel", envir = ck) else list()
}

# ctx: study settings (task, str.grid, n.simset, mcmcpar, paths); file: result file;
# chunk_size: tasks per checkpoint (NA: chunk_mult x number of workers); is_done(r): TRUE if result r is complete;
# sim_indices: (i, j) pairs of scenario and replicate (default: all)
run_study <- function(ctx, file, chunk_size, is_done, sim_indices = NULL, chunk_mult = 4L) {
  if (is.null(sim_indices))
    sim_indices <- expand.grid(i = seq_len(nrow(ctx$str.grid)), j = seq_len(ctx$n.simset))
  wall_secs <- as.numeric(Sys.getenv("SIM_WALL_HOURS", "47")) * 3600

  if (SIM_LOCAL) {
    suppressMessages({ library(foreach); library(doRNG) })
    registerDoSEQ(); n_workers <- 1L
  } else {
    suppressMessages({ library(Rmpi); library(doMPI); library(doRNG) })
    cl <- startMPIcluster()
    registerDoMPI(cl); n_workers <- cl$workerCount
  }
  if (is.na(chunk_size)) chunk_size <- max(1L, as.integer(chunk_mult * n_workers))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  results_parallel <- Filter(function(r) is.list(r) && isTRUE(is_done(r)), load_checkpoint(file))
  design <- ctx
  save_ck <- function() atomic_save(c("results_parallel", "design", "sim_indices"), file, environment())
  done <- vapply(results_parallel, function(r) paste0(r$i, "_", r$j), character(1))
  todo <- which(!(paste0(sim_indices$i, "_", sim_indices$j) %in% done))
  cat(sprintf("%s | %d workers | %d of %d tasks done | chunk size %d | output %s\n",
              ctx$label, n_workers, length(done), nrow(sim_indices), chunk_size, file))

  t_start <- proc.time()[["elapsed"]]; chunk_secs <- numeric(0)
  for (cs in if (length(todo)) seq(1L, length(todo), by = chunk_size) else integer(0)) {
    elapsed <- proc.time()[["elapsed"]] - t_start
    projected <- elapsed + if (length(chunk_secs)) 1.15 * max(chunk_secs) else 0
    if (projected > wall_secs) {
      cat(sprintf("Stopping at %.1f h (next chunk expected to end at %.1f h). Resubmit the job.\n",
                  elapsed / 3600, projected / 3600))
      break
    }
    t_chunk <- proc.time()[["elapsed"]]
    chunk <- sim_indices[todo[cs:min(cs + chunk_size - 1L, length(todo))], , drop = FALSE]
    chunk_res <- foreach(idx = seq_len(nrow(chunk)), .export = c("ctx", "chunk"),
                         .packages = c("mvtnorm", "survival", "glmnet"),
                         .options.RNG = 1234L + cs) %dorng% {
      source(ctx$paths$tasks, local = FALSE)
      sim_task(chunk$i[idx], chunk$j[idx], ctx)
    }
    results_parallel <- c(results_parallel, unname(chunk_res))
    save_ck()
    chunk_secs <- c(chunk_secs, proc.time()[["elapsed"]] - t_chunk)
    n_err <- sum(vapply(chunk_res, function(r) !isTRUE(is_done(r)), logical(1)))
    cat(sprintf("[%s] tasks %d-%d of %d remaining done in %.2f h (%d failed). Stored: %d / %d.\n",
                format(Sys.time(), "%Y-%m-%d %H:%M"), cs, cs + nrow(chunk) - 1L, length(todo),
                tail(chunk_secs, 1) / 3600, n_err,
                sum(vapply(results_parallel, function(r) isTRUE(is_done(r)), logical(1))), nrow(sim_indices)))
  }
  n_ok <- sum(vapply(results_parallel, function(r) isTRUE(is_done(r)), logical(1)))
  save_ck()
  cat(sprintf("Finished this job: %d / %d tasks stored in %s.%s\n", n_ok, nrow(sim_indices), file,
              if (n_ok < nrow(sim_indices)) " Resubmit the job to continue." else ""))
  if (!SIM_LOCAL) { closeCluster(cl); Rmpi::mpi.quit() }
  invisible(results_parallel)
}
