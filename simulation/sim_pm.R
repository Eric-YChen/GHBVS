# Usage (from simulation/): Rscript sim_pm.R [n_cores]   or   Rscript sim_pm.R plot   (figures from results/sim_pm.csv)
suppressMessages({ library(mvtnorm); library(parallel) })
paths <- list(sim = normalizePath(".", winslash = "/"), core = normalizePath("../core", winslash = "/"))
quick <- nzchar(Sys.getenv("SIM_QUICK"))
args <- commandArgs(trailingOnly = TRUE)
plot.only <- identical(args[1], "plot")
n.cores <- if (length(args) >= 1 && !plot.only) as.integer(args[1]) else
  as.integer(Sys.getenv("NSLOTS", max(1L, parallel::detectCores() - 1L)))

struct.grid <- c("PH", "AFT", "AH", "GH")
cens.grid   <- plogis(c(1, 2, 3, 4, 5) + 0.5)
n.data      <- 2000L
n.seed      <- if (quick) 3L else 100L
priors      <- c("LCM", "Product")
scalings    <- c(orig = FALSE, calib = TRUE)
anchor.struct <- "PH"; beta.active <- 0.5
calib.ref.cens <- 0.99
gh.ratio <- 0.5
mu0 <- 1.55; sig0 <- 0.7; gE <- 1; n.quad <- 120L
res.dir <- file.path(paths$sim, "results"); dir.create(res.dir, showWarnings = FALSE)
suffix <- if (quick) "_quick" else ""

if (!plot.only) {
struct2gam <- c(PH = 2, AFT = 4, AH = 1, GH = 3)
beta_dir <- function(struct, b) switch(struct, PH = c(0, b), AFT = c(b, b), AH = c(b, 0),
                                       GH = c(b, gh.ratio * b))
gauss_hermite <- function(K) {
  i <- seq_len(K - 1L); J <- matrix(0, K, K)
  J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2 / sum(e$vectors[1, ]^2))
}
GHQ <- gauss_hermite(n.quad)
S0 <- function(t) pnorm((log(t) - mu0) / sig0, lower.tail = FALSE)
cens_marg <- function(struct, b, tau) {
  bb <- beta_dir(struct, b); s <- pmin(pmax(S0(tau * exp(bb[1] * GHQ$x)), 1e-300), 1)
  sum(GHQ$w * exp(exp((bb[2] - bb[1]) * GHQ$x) * log(s)))
}
solve_beta <- function(struct, target, tau, b.max = 4, n.scan = 80) {
  f <- function(b) cens_marg(struct, b, tau) - target
  bg <- seq(0, b.max, length.out = n.scan); fg <- vapply(bg, f, numeric(1))
  for (k in seq_len(n.scan - 1L)) if (fg[k] * fg[k + 1L] <= 0) return(uniroot(f, bg[k:(k + 1L)], tol = 1e-10)$root)
  stop("calibration target unreachable for ", struct)
}
tau.ref <- qlnorm(1 - calib.ref.cens, mu0, sig0)
target.ref <- cens_marg(anchor.struct, beta.active, tau.ref)
beta.fixed <- vapply(struct.grid, function(s) solve_beta(s, target.ref, tau.ref), numeric(1))
cat("calibrated coefficients b (GH: beta_t = b, beta_h = 0.5 b):\n"); print(round(beta.fixed, 3))
calib <- expand.grid(struct = struct.grid, cens = cens.grid, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
calib$cutoff <- qlnorm(1 - calib$cens, mu0, sig0)
calib$beta <- beta.fixed[calib$struct]
calib$cens.active <- mapply(cens_marg, calib$struct, calib$beta, calib$cutoff)

ess_cgc <- function(nu, theta0, cens_times) {
  z <- as.numeric(exp(nu) * log(cens_times) - theta0)
  Pz <- pnorm(z); pz <- dnorm(z)
  hz <- ifelse(z > 8, z + 1 / z, pz / pmax(1 - Pz, .Machine$double.xmin))
  sum(pmax(pmin(Pz + pz * (hz - z), 1), 0))
}
mk_par <- function(gam, scale, n)
  list(nu = 0, theta0 = 0, theta = 0, eta = 0, gam = gam, g = scale * gE, n = scale,
       g.theta = scale / n, g.eta = scale / n, nu.alpha = 0.01, nu.beta = 0.01)
compute_cell <- function(task) {
  set.seed(task$seed)
  x <- matrix(rnorm(task$n), task$n, 1)
  bb <- if (task$active) beta_dir(task$struct, task$beta) else c(0, 0)
  Tt <- simGH(seed = task$seed + 1000000L, n = task$n, des_h = x, des_t = x, theta = c(mu0, sig0),
              beta_h = bb[2], beta_t = bb[1], hstr = "GH", baseline = "LN")
  Tt <- pmax(Tt, 1e-8); Tt[!is.finite(Tt)] <- Inf
  st <- as.numeric(Tt <= task$cutoff)
  dat <- list(x = x, xt = x, y = pmin(Tt, task$cutoff), status = st)
  dat$xt2 <- dat$x2 <- crossprod(x); dat$xt2.inv <- dat$x2.inv <- solve(dat$x2)
  scale <- task$n
  if (task$ess) {
    f0 <- tryCatch(mle.gh(mk_par(0, task$n, task$n), dat), error = function(e) NULL)
    if (!is.null(f0)) scale <- ess_cgc(f0$nu, f0$theta0, rep(task$cutoff, task$n))
  }
  lm <- function(gam) {
    par <- mk_par(gam, scale, task$n)
    tryCatch(if (task$prior == "LCM") marg.loglik(par, dat) else lapprox.margin.likelihood(par, dat),
             error = function(e) NA_real_)
  }
  c(nevent = sum(st), N = scale, logBF = lm(struct2gam[[task$struct]]) - lm(0))
}

grid <- expand.grid(seed = seq_len(n.seed), cens = cens.grid, active = c(FALSE, TRUE),
                    struct = struct.grid, prior = priors, scaling = names(scalings),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
grid <- merge(grid, calib[c("struct", "cens", "cutoff", "beta")], by = c("struct", "cens"), sort = FALSE)
grid$n <- n.data; grid$ess <- scalings[grid$scaling]
tasks <- lapply(seq_len(nrow(grid)), function(r) as.list(grid[r, ]))
setup <- function(core) {
  assign("BVS_CORE", core, envir = .GlobalEnv)
  suppressMessages({ library(mvtnorm)
    source(file.path(core, "ADS_sampler.R")); source(file.path(core, "utils", "EB_gprior_fixed.R"))
    source(file.path(core, "utils", "gprior_setup_fixed.R")) })
  invisible(TRUE)
}
t0 <- Sys.time()
if (n.cores > 1L) {
  cl <- makeCluster(n.cores)
  clusterExport(cl, c("beta_dir", "gh.ratio", "struct2gam", "mu0", "sig0", "gE", "ess_cgc", "mk_par", "compute_cell"))
  clusterCall(cl, setup, paths$core)
  out <- parLapplyLB(cl, tasks, compute_cell); stopCluster(cl)
} else { setup(paths$core); out <- lapply(tasks, compute_cell) }
out <- do.call(rbind, out)
grid$nevent <- out[, "nevent"]; grid$N <- out[, "N"]; grid$logBF <- out[, "logBF"]
grid$PIP <- plogis(grid$logBF)
cat(sprintf("%d replicates in %.1f min on %d cores; %d failed\n", nrow(grid),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n.cores, sum(is.na(grid$logBF))))
write.csv(grid, file.path(res.dir, paste0("sim_pm", suffix, ".csv")), row.names = FALSE)
} else grid <- read.csv(file.path(res.dir, paste0("sim_pm", suffix, ".csv")), stringsAsFactors = FALSE)

ess.tab <- aggregate(cbind(neff.over.n = N / n) ~ cens, data = grid[grid$scaling == "calib" & !grid$active, ], FUN = mean)
cat("\nn_eff / n by censoring level (inactive covariate):\n"); print(round(ess.tab, 3))
for (sc in names(scalings)) {
  pdf(file.path(res.dir, paste0("pm_", sc, suffix, ".pdf")), width = 12, height = 6)
  op <- par(mfrow = c(length(priors), length(struct.grid)), mar = c(4, 4, 3, 1))
  cols <- c("FALSE" = "#2a78d6", "TRUE" = "#eb6834")
  for (pr in priors) for (s in struct.grid) {
    plot(NA, xlim = range(cens.grid), ylim = c(0, 1), xlab = "censoring proportion",
         ylab = "posterior inclusion probability", main = sprintf("%s prior - %s", pr, s))
    abline(h = 0.5, lty = 3, col = "grey50")
    for (ac in c(FALSE, TRUE)) {
      sub <- grid[grid$scaling == sc & grid$prior == pr & grid$struct == s & grid$active == ac, ]
      cv <- sort(unique(sub$cens))
      q <- sapply(cv, function(c0) { z <- sub$PIP[sub$cens == c0]
        c(mean(z, na.rm = TRUE), quantile(z, c(0.25, 0.75), na.rm = TRUE)) })
      polygon(c(cv, rev(cv)), c(q[2, ], rev(q[3, ])), col = adjustcolor(cols[as.character(ac)], 0.15), border = NA)
      lines(cv, q[1, ], col = cols[as.character(ac)], lwd = 2); points(cv, q[1, ], col = cols[as.character(ac)], pch = 19)
    }
    legend("topright", bty = "n", cex = 0.8, legend = c("full model (active)", "null model (inactive)"),
           col = cols[c("TRUE", "FALSE")], lwd = 2, pch = 19)
  }
  par(op); dev.off()
}
cat("saved results/sim_pm", suffix, ".csv and results/pm_orig", suffix, ".pdf, pm_calib", suffix, ".pdf\n", sep = "")
