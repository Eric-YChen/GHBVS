mpi_rank <- function() {
  for (v in c("SLURM_PROCID", "OMPI_COMM_WORLD_RANK", "PMI_RANK",
              "MV2_COMM_WORLD_RANK")) {
    r <- Sys.getenv(v)
    if (nzchar(r)) {
      n <- suppressWarnings(as.integer(r))
      if (!is.na(n)) return(n)
    }
  }
  ## Fall back to Rmpi if it is loaded and initialised.
  if ("Rmpi" %in% loadedNamespaces()) {
    n <- suppressWarnings(try(Rmpi::mpi.comm.rank(0L), silent = TRUE))
    if (!inherits(n, "try-error") && is.numeric(n)) return(as.integer(n))
  }
  0L
}

is_rank0 <- function() mpi_rank() == 0L

## Write `names` (a character vector of object names in `envir`) to `file`,
## atomically. Returns the path invisibly.
atomic_save <- function(names, file, envir = parent.frame(), compress = TRUE) {
  tmp <- paste0(file, ".tmp-", Sys.info()[["nodename"]], "-", Sys.getpid(),
                "-", basename(tempfile("")))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  save(list = names, file = tmp, envir = envir, compress = compress)

  if (!file.rename(tmp, file)) {
    warning("atomic_save(): rename failed, falling back to a non-atomic copy.")
    if (!file.copy(tmp, file, overwrite = TRUE))
      stop("atomic_save(): could not move ", tmp, " into place as ", file)
    unlink(tmp)
  }
  invisible(file)
}


load_checkpoint_results <- function(file, what = "results_parallel") {
  if (!file.exists(file))
    return(list(results = list(), status = "missing",
                msg = NA_character_, quarantined = NA_character_))

  ck  <- new.env(parent = emptyenv())
  err <- tryCatch({ load(file, envir = ck); NULL },
                  error = function(e) conditionMessage(e))

  if (!is.null(err)) {
    quar <- paste0(file, ".corrupt-", format(Sys.time(), "%Y%m%d-%H%M%S"))
    if (!file.rename(file, quar)) quar <- NA_character_
    return(list(results = list(), status = "corrupt",
                msg = err, quarantined = quar))
  }

  res <- if (exists(what, envir = ck, inherits = FALSE))
           get(what, envir = ck) else list()
  list(results = res, status = "ok",
       msg = NA_character_, quarantined = NA_character_)
}

report_checkpoint <- function(ck, file) {
  if (identical(ck$status, "corrupt")) {
    cat(strrep("*", 78), "\n", sep = "")
    cat("*** CHECKPOINT UNREADABLE: ", file, "\n", sep = "")
    cat("***   ", ck$msg, "\n", sep = "")
    if (!is.na(ck$quarantined))
      cat("***   moved aside to: ", ck$quarantined, "\n", sep = "")
    cat("*** This usually means a previous job was killed mid-save.\n")
    cat("*** Starting from scratch; completed work in that file is lost.\n")
    cat(strrep("*", 78), "\n", sep = "")
  } else if (identical(ck$status, "missing")) {
    cat("No checkpoint found; starting from scratch.\n")
  }
  invisible(ck)
}
