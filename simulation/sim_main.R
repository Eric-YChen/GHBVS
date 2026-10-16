# Usage (from simulation/, under mpirun/srun): Rscript sim_main.R <LCM|Product> <p> [default|unif|restricted] [chunk_size]
source("sim_setup.R")
args    <- commandArgs(trailingOnly = TRUE)
prior   <- if (length(args) >= 1) args[1] else "LCM"
p       <- as.integer(if (length(args) >= 2) args[2] else 10)
variant <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "default"
chunk   <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
stopifnot(prior %in% c("LCM", "Product"), p %in% c(10L, 50L, 100L),
          variant %in% c("default", "unif", "restricted"))
if (variant == "restricted" && p != 10L) stop("the restricted-structure study uses p = 10")

str.grid <- make_grid(p, n.vec = if (variant == "restricted") c(250L, 500L) else c(250L, 500L, 1000L))
if (variant == "restricted") {
  keep <- with(str.grid, hstr %in% c("PH", "AFT") &
                 ((baseline == "LN" & cens == "0.25default") | (baseline == "PGW" & cens == "0.5default") |
                  (baseline == "PGW" & cens == "0.25default" & hstr == "AFT")))
  str.grid <- str.grid[keep, ]; rownames(str.grid) <- NULL
}

ctx <- list(label = sprintf("sim_main %s p = %d (%s)", prior, p, variant), task = "task_main",
            str.grid = str.grid, n.simset = if (SIM_QUICK) 2L else 250L,
            mcmcpar = if (SIM_QUICK) list(n.iter = 400, n.burnin = 200, n.thin = 2, n.chain = 1, printevery = 1e9)
                      else list(n.iter = 20000, n.burnin = 10000, n.thin = 2, n.chain = 1, printevery = 1e9),
            prior = prior, variant = variant, paths = sim_paths())
file <- file.path(ctx$paths$results, sprintf("sim_main_%s_p%d_%s%s.RData", prior, p, variant,
                                             if (SIM_QUICK) "_quick" else ""))
run_study(ctx, file, chunk, is_done = function(r) is.null(r[["error"]]) && !is.null(r[["hpm"]]))
