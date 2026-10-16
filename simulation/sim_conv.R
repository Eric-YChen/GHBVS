# Usage (from simulation/, under mpirun/srun): Rscript sim_conv.R <10|50|100> [chunk_size]
source("sim_setup.R")
args  <- commandArgs(trailingOnly = TRUE)
p     <- as.integer(if (length(args) >= 1) args[1] else 10)
chunk <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_
stopifnot(p %in% c(10L, 50L, 100L))

n.iter <- if (SIM_QUICK) 600L else 50000L
ctx <- list(label = sprintf("sim_conv p = %d", p), task = "task_conv", str.grid = make_grid(p),
            n.simset = if (SIM_QUICK) 1L else 10L, prior = "LCM",
            mcmcpar = list(n.iter = n.iter, n.burnin = 0, n.thin = 1, n.chain = 1, printevery = 1e9),
            gr_settings = list(s50k = list(n.iter = n.iter, burnin = n.iter / 2, thin = 2),
                               s20k = list(n.iter = 0.4 * n.iter, burnin = 0.2 * n.iter, thin = 2)),
            init_types = c("init.gam", "init.gam", "rand_sparse", "rand_sparse", "rand_sparse",
                           "rand_sparse", "underdispersed", "overdispersed"),
            paths = sim_paths())
file <- file.path(ctx$paths$results, sprintf("sim_conv_p%d%s.RData", p, if (SIM_QUICK) "_quick" else ""))
run_study(ctx, file, chunk, is_done = function(r) is.null(r[["error"]]) && !is.null(r[["gr"]]), chunk_mult = 1L)
