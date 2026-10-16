# Usage (from simulation/, under mpirun/srun): Rscript sim_bic.R <10|50> [n values, comma-separated; default 250,1000] [chunk_size]
source("sim_setup.R")
args  <- commandArgs(trailingOnly = TRUE)
p     <- as.integer(if (length(args) >= 1) args[1] else 10)
n.vec <- if (length(args) >= 2) as.integer(strsplit(args[2], ",", fixed = TRUE)[[1]]) else c(250L, 1000L)
chunk <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
stopifnot(p %in% c(10L, 50L), length(n.vec) > 0, all(n.vec %in% c(250L, 500L, 1000L)))
n.tag <- if (identical(n.vec, c(250L, 1000L))) "" else paste0("_n", paste(n.vec, collapse = "-"))

methods <- if (p == 10L) list(bic_flat = "unif", bic_gamma = "default") else list(bic_gamma = "default")
ctx <- list(label = sprintf("sim_bic p = %d, n = %s (%s)", p, paste(n.vec, collapse = ", "), paste(names(methods), collapse = ", ")),
            task = "task_bic", str.grid = make_grid(p, n.vec = n.vec),
            n.simset = if (SIM_QUICK) 2L else 250L,
            mcmcpar = if (SIM_QUICK) list(n.iter = 400, n.burnin = 200, n.thin = 2, n.chain = 1, printevery = 1e9)
                      else list(n.iter = 20000, n.burnin = 10000, n.thin = 2, n.chain = 1, printevery = 1e9),
            methods = methods, paths = sim_paths())
file <- file.path(ctx$paths$results, sprintf("sim_bic_p%d%s%s.RData", p, n.tag, if (SIM_QUICK) "_quick" else ""))
run_study(ctx, file, chunk,
          is_done = function(r) is.null(r[["error"]]) && length(r[["fits"]]) == length(methods))
