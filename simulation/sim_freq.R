# Usage (from simulation/, under mpirun/srun): Rscript sim_freq.R <10|50|100> [BIC|AIC|LASSO] [chunk_size]
source("sim_setup.R")
args  <- commandArgs(trailingOnly = TRUE)
p     <- as.integer(if (length(args) >= 1) args[1] else 10)
crit  <- toupper(if (length(args) >= 2 && nzchar(args[2])) args[2] else "BIC")
chunk <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
stopifnot(p %in% c(10L, 50L, 100L), crit %in% c("BIC", "AIC", "LASSO"))

ctx <- list(label = sprintf("sim_freq %s p = %d", crit, p), task = "task_freq", str.grid = make_grid(p),
            n.simset = if (SIM_QUICK) 2L else 250L, crit = crit, paths = sim_paths())
file <- file.path(ctx$paths$results, sprintf("sim_freq_%s_p%d%s.RData", tolower(crit), p,
                                             if (SIM_QUICK) "_quick" else ""))
run_study(ctx, file, chunk, is_done = function(r) is.null(r[["error"]]) && !is.null(r[["sel"]]))
