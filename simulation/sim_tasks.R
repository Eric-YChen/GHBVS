sim_sigma <- function(p) 0.7 ^ abs(outer(seq_len(p), seq_len(p), "-"))

true_model <- function(p, hstr) {
  idx <- as.integer(quantile(seq_len(p), c(0.2, 0.4, 0.6, 0.8)))
  g <- rep(0L, p); g[idx] <- c(AH = 1L, PH = 2L, GH = 3L, AFT = 4L)[[hstr]]
  paste0(g, collapse = "")
}

secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

load_core <- function(ctx) {
  suppressMessages({ library(mvtnorm); library(survival); library(glmnet) })
  assign("BVS_CORE", ctx$paths$core, envir = .GlobalEnv)
  if (!exists("BVS.MCMC", envir = .GlobalEnv, inherits = FALSE))
    suppressMessages(source(file.path(ctx$paths$core, "ADS_sampler.R"), local = FALSE))
}

sim_data <- function(i, j, ctx) {
  g <- ctx$str.grid[i, ]
  d <- simdat(n = g$n, p = g$p, sigma = sim_sigma(g$p), hstr = g$hstr, baseline = g$baseline,
              n.simset = 1, cens = g$cens, seeds = j)
  list(dat = d$data[[1]], cens_rat = d$cens.rat[1])
}

model_par <- function(dat, gam, gam.param = NULL) {
  list(inits = list(nu = rnorm(1), theta0 = rnorm(1), theta = rnorm(ncol(dat$xt), 0, 0.5),
                    eta = rnorm(ncol(dat$x), 0, 0.5), g.theta = 1, g.eta = 1, g = 1,
                    a_g = NULL, gam = gam, nu.alpha = 0.01, nu.beta = 0.01),
       gam.param = gam.param, data = dat)
}

# ch: BVS.MCMC output; ntop: number of top models kept
summarise_fit <- function(ch, ntop = 200L) {
  pp <- ch$post.prob
  incl <- do.call(rbind, lapply(strsplit(pp$Vars, "", fixed = TRUE), function(s) s != "0"))
  w <- pp$prob; v <- pp$visits / sum(pp$visits)
  hz <- factor(pp$haz, levels = c("Null", "AH", "PH", "AFT", "GH"))
  by_haz <- function(wt) { h <- tapply(wt, hz, sum); h[is.na(h)] <- 0; h }
  pip_haz <- t(sapply(c("AH", "PH", "AFT", "GH"), function(l) {
    s <- pp$haz == l
    if (!any(s) || sum(w[s]) <= 0) rep(NA_real_, ncol(incl)) else colSums(w[s] * incl[s, , drop = FALSE]) / sum(w[s])
  }))
  list(hpm = pp$Vars[1],
       hpm.visited = pp$Vars[which(pp$visits > 0)[1]],
       PIP.post = colSums(w * incl), PIP.freq = colSums(v * incl),
       haz.post = by_haz(w), haz.freq = by_haz(v), PIP.haz = pip_haz,
       top = head(pp, ntop), n.eval = nrow(pp), n.visited = sum(pp$visits > 0),
       acc = mean(ch$acc[, 1]), opt.log = ch$opt.log)
}

gam_param <- function(ctx, i) {
  switch(ctx$variant,
         default    = NULL,
         unif       = list(layer1.betabin.par = "uniform"),
         restricted = list(layer2.mulnom.par = switch(ctx$str.grid$hstr[i],
                                                      PH = c(0, 1, 0, 0), AFT = c(0, 0, 0, 1))))
}

task_main <- function(i, j, ctx) {
  load_core(ctx)
  s <- sim_data(i, j, ctx)
  t0 <- Sys.time()
  init <- suppressWarnings(init.gam(s$dat, history.keep = FALSE))
  if (ctx$variant == "restricted")
    init[init != 0] <- c(PH = 2L, AFT = 4L)[[ctx$str.grid$hstr[i]]]
  t_init <- secs(t0)
  gp <- gam_param(ctx, i)
  t0 <- Sys.time()
  ch <- BVS.MCMC(model_par(s$dat, init, gp), ctx$mcmcpar, verbose = FALSE,
                 gprior.mat = ctx$prior, gamprior = if (is.null(gp)) "default" else "custom")
  c(list(i = i, j = j, cens_rat = s$cens_rat, init = paste0(init, collapse = "")),
    summarise_fit(ch), list(runtime = c(init = t_init, sampler = secs(t0))))
}

# cache: environment for the MLE of each evaluated model
make_target_bic <- function(cache) {
  function(par, data, ...) {
    if (all(par$gam == 0)) return(-Inf)
    lp <- prior.gam(par, log = TRUE)
    if (!is.finite(lp)) return(-Inf)
    key <- paste0(par$gam, collapse = "")
    ll <- get0(key, envir = cache, inherits = FALSE)
    if (is.null(ll)) {
      fit <- tryCatch(mle.gh(par, data), error = function(e) NULL)
      ll <- if (is.null(fit)) -Inf else tryCatch(lglik(fit, data), error = function(e) -Inf)
      if (!is.finite(ll)) ll <- -Inf
      assign(key, ll, envir = cache)
    }
    d <- sum(par$gam %in% c(1, 2, 4)) + 2 * sum(par$gam == 3)
    lp + ll - 0.5 * (d + 2) * log(nrow(data$x))
  }
}

task_bic <- function(i, j, ctx) {
  load_core(ctx)
  s <- sim_data(i, j, ctx)
  init <- suppressWarnings(init.gam(s$dat, history.keep = FALSE))
  cache <- new.env(hash = TRUE)
  fits <- list(); rt <- c(); err <- c()
  for (m in names(ctx$methods)) {
    t0 <- Sys.time()
    r <- tryCatch({
      ch <- BVS.MCMC(model_par(s$dat, init), ctx$mcmcpar, verbose = FALSE, gprior.mat = "LCM",
                     gamprior = ctx$methods[[m]], target = make_target_bic(cache))
      summarise_fit(ch)
    }, error = function(e) conditionMessage(e))
    if (is.character(r)) err[m] <- r else fits[[m]] <- r
    rt[m] <- secs(t0)
  }
  if (!length(fits)) stop(paste(names(err), err, collapse = "; "))
  list(i = i, j = j, cens_rat = s$cens_rat, init = paste0(init, collapse = ""),
       fits = fits, runtime = rt, method_errors = err)
}

# px: number of covariates; m: number of variables
rand_valid_gam <- function(px, m) {
  m <- max(1L, min(as.integer(m), px)); idx <- sample.int(px, m); gam <- rep(0L, px)
  strc <- sample(c("AH", "PH", "AFT", "GH"), 1L)
  if (strc == "GH") {
    gam[idx] <- sample(1:3, m, replace = TRUE)
    if (all(gam[idx] == 1L) || all(gam[idx] == 2L)) gam[idx[1]] <- 3L
  } else gam[idx] <- c(AH = 1L, PH = 2L, AFT = 4L)[[strc]]
  gam
}

# chains: list of inclusion-indicator matrices, one per chain
compute_rhat <- function(chains) {
  chains <- Filter(function(m) !is.null(m) && nrow(m) >= 20L, chains)
  if (length(chains) < 2L) return(NULL)
  n <- min(vapply(chains, nrow, integer(1))); m <- length(chains)
  vapply(seq_len(ncol(chains[[1]])), function(k) {
    x <- vapply(chains, function(ch) ch[seq_len(n), k], numeric(n))
    W <- mean(apply(x, 2, var)); B <- n / (m - 1) * sum((colMeans(x) - mean(x))^2)
    if (W < .Machine$double.eps) NA_real_ else sqrt(((n - 1) / n * W + B / n) / W)
  }, numeric(1))
}

task_conv <- function(i, j, ctx) {
  load_core(ctx)
  s <- sim_data(i, j, ctx); px <- ncol(s$dat$x)
  t0 <- Sys.time()
  init <- tryCatch(suppressWarnings(init.gam(s$dat, history.keep = FALSE)),
                   error = function(e) { g <- rep(0L, px); g[1] <- 2L; g })
  t_init <- secs(t0)
  m0 <- max(1L, sum(init != 0))
  inits <- lapply(ctx$init_types, function(type) switch(type,
    init.gam = as.integer(init), rand_sparse = rand_valid_gam(px, m0),
    underdispersed = rand_valid_gam(px, 1L), overdispersed = rand_valid_gam(px, min(2L * m0, px))))
  n_ch <- length(inits); incl <- vector("list", n_ch)
  acc <- rep(NA_real_, n_ch); err <- rep(NA_character_, n_ch); rt <- rep(NA_real_, n_ch)
  for (k in seq_len(n_ch)) {
    t0 <- Sys.time()
    out <- tryCatch(BVS.MCMC(model_par(s$dat, inits[[k]]), ctx$mcmcpar, verbose = FALSE,
                             gprior.mat = ctx$prior, gamprior = "default"),
                    error = function(e) list(err = conditionMessage(e)))
    rt[k] <- secs(t0)
    if (!is.null(out$err)) { err[k] <- out$err; next }
    incl[[k]] <- as.matrix(out$chain[[1]]) != 0
    acc[k] <- mean(out$acc[, 1]); rm(out)
  }
  practical <- which(ctx$init_types %in% c("init.gam", "rand_sparse"))
  gr <- lapply(ctx$gr_settings, function(st) {
    keep <- seq(st$burnin + st$thin, st$n.iter, st$thin)
    B <- lapply(incl, function(x) if (is.null(x)) NULL else x[keep, , drop = FALSE] * 1)
    all <- compute_rhat(B); prac <- compute_rhat(B[practical])
    pip <- t(vapply(B, function(x) if (is.null(x)) rep(NA_real_, px) else colMeans(x), numeric(px)))
    list(psrf = all, psrf.practical = prac,
         mean = if (is.null(all)) NA else mean(all, na.rm = TRUE),
         max = if (is.null(all) || all(is.na(all))) NA else max(all, na.rm = TRUE),
         mean.practical = if (is.null(prac)) NA else mean(prac, na.rm = TRUE),
         max.practical = if (is.null(prac) || all(is.na(prac))) NA else max(prac, na.rm = TRUE),
         pip = pip)
  })
  if (sum(is.na(err)) < 2L) stop("fewer than two chains ran: ", paste(na.omit(err), collapse = "; "))
  list(i = i, j = j, cens_rat = s$cens_rat, gr = gr, acc = acc, chain_errors = err,
       init_models = vapply(inits, paste0, character(1), collapse = ""), n_chains_ok = sum(is.na(err)),
       runtime = c(init = t_init, chains = rt))
}

task_freq <- function(i, j, ctx) {
  load_core(ctx)
  if (!exists("step.gh.freq", envir = .GlobalEnv, inherits = FALSE))
    source(file.path(ctx$paths$sim, "utils", "freq_stepwise_multi.R"), local = FALSE)
  s <- sim_data(i, j, ctx); tru <- ctx$str.grid$true.mod[i]
  t0 <- Sys.time()
  f <- step.gh.freq(s$dat, crit = ctx$crit, lambda = "cox", dist = "LogNormal", maxit = 5000L,
                    max.steps = min(2L * ncol(s$dat$x), 40L), max.active = min(ncol(s$dat$x), 20L),
                    gh.strict = FALSE, verbose = FALSE)
  b <- f$best; mm <- metrics(b$gam, tru)
  list(i = i, j = j, cens_rat = s$cens_rat,
       sel = list(gam.str = b$gam.str, structure = b$structure, hstr = b$hstr, score = b$score,
                  loglik = b$loglik, k = b$k, AIC = b$AIC, BIC = b$BIC),
       arm.table = f$table, coef = b$coef,
       metrics = c(Sens = mm$Sens, Spec = mm$Spec, FDR = mm$FDR, phi = mm$phi,
                   haz_correct = as.numeric(b$hstr == gam2haz(tru))),
       runtime = c(selection = secs(t0)))
}

# i: scenario (row of ctx$str.grid); j: replicate; ctx: study settings
sim_task <- function(i, j, ctx) {
  tryCatch(get(ctx$task)(i, j, ctx),
           error = function(e) list(i = i, j = j, error = paste0(conditionMessage(e), " [scenario ", i,
                                                                  ", replicate ", j, "]")))
}
