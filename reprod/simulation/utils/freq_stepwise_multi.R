## ==========================================================================
## 0.  Low-level helpers
## ==========================================================================

## Number of baseline parameters GHMLE() fits for each distribution, i.e. the
## offset before the first regression coefficient in OPT$par.
.n_base <- function(dist) {
  if (dist %in% c("PGW", "EW", "GenGamma")) 3L
  else if (dist %in% c("LogNormal", "LogLogistic", "Gamma", "Weibull")) 2L
  else stop("Unsupported baseline distribution: ", dist)
}

## Moment-based starting values for the baseline parameters.  Starting all of
## them at 0 (as some of the older routines do) is a poor guess whenever the
## survival times are not of order 1.
.start_base <- function(times, dist) {
  lt <- log(pmax(as.vector(times), 1e-8))
  m  <- mean(lt)
  s  <- stats::sd(lt)
  if (!is.finite(s) || s <= 0) s <- 1
  switch(dist,
         LogNormal   = c(m, log(s)),
         LogLogistic = c(m, log(s * sqrt(3) / pi)),
         Weibull     = c(m, 0),                # (log scale, log shape)
         Gamma       = c(m, 0),
         PGW         = c(m, 0, 0),             # (log sigma, log nu, log gamma)
         EW          = c(m, 0, 0),
         GenGamma    = c(m, 0, 0),
         stop("Unsupported baseline distribution: ", dist))
}

## Smoothed |b|, differentiable at 0.
.abs_smooth <- function(b) sqrt(b * b + 1e-10)

## nlminb reports $objective, optim reports $value.
.opt_obj <- function(OPT) if (!is.null(OPT$objective)) OPT$objective else OPT$value

## paste0() recycles zero-length arguments to "", so paste0("t", integer(0))
## is "t", not character(0).  Guard against that: an empty index set must
## contribute no names at all, otherwise phantom parameters creep into `init`.
.nm <- function(prefix, idx) if (length(idx)) paste0(prefix, idx) else character(0)

## GHMLE names parameters positionally.  These are the names we attach, in the
## order GHMLE lays them out: baseline, then des_t coefficients (alpha), then
## des coefficients (beta).  AFT lives on its own arm so "a" never co-occurs
## with "t"/"h".
.par_names <- function(t_idx, h_idx, aft_idx, dist) {
  c(paste0("base", seq_len(.n_base(dist))),
    .nm("t", t_idx), .nm("h", h_idx), .nm("a", aft_idx))
}

## Hazard structure implied by a set of active indices.
.hstr_of <- function(t_idx, h_idx, aft_idx) {
  if (length(aft_idx))                        "AFT"
  else if (length(t_idx) && length(h_idx))    "GH"
  else if (length(t_idx))                     "AH"
  else if (length(h_idx))                     "PH"
  else                                        "baseline"
}


## ==========================================================================
## 1.  Fit a single candidate model
## ==========================================================================
#
#   t_idx    columns of xt carrying a time-level effect   (codes 1 and 3)
#   h_idx    columns of x  carrying a hazard-level effect (codes 2 and 3)
#   aft_idx  columns of x  carrying an AFT effect         (code 4)
#
# Exactly one of {aft_idx} and {t_idx, h_idx} may be non-empty.
#
#   lambda   L1 penalty on the regression coefficients; 0 for AIC/BIC
#   start    optional *named* parameter vector from a previously fitted model,
#            used as a warm start.  Names are matched, so a coefficient keeps
#            its value when a different covariate is added elsewhere and the
#            newly added one starts at 0.  This both speeds the search up and
#            makes deep steps far more stable than restarting from scratch.
#
# Returns a list with ok / loglik / k / par / l1 / hstr.  ok = FALSE (and an
# infinite score downstream) whenever the optimiser fails or the objective is
# not finite, so a failed candidate can never win the argmin.
.fit_one <- function(t_idx, h_idx, aft_idx, data, dist, lambda = 0,
                     maxit = 5000, start = NULL) {

  t_idx   <- sort(as.integer(t_idx))
  h_idx   <- sort(as.integer(h_idx))
  aft_idx <- sort(as.integer(aft_idx))

  nb   <- .n_base(dist)
  nm   <- .par_names(t_idx, h_idx, aft_idx, dist)
  hstr <- .hstr_of(t_idx, h_idx, aft_idx)

  if (hstr == "GH" && dist == "Weibull")
    stop("GHMLE() has no GH branch for the Weibull baseline.")

  ## -- initial value ---------------------------------------------------------
  init        <- c(.start_base(data$y, dist), rep(0, length(nm) - nb))
  names(init) <- nm
  if (!is.null(start)) {
    hit <- intersect(names(start), nm)
    if (length(hit)) init[hit] <- start[hit]
  }

  des   <- if (hstr == "AFT")     data$x[,  aft_idx, drop = FALSE]
           else if (length(h_idx)) data$x[,  h_idx,   drop = FALSE]
           else                    NULL
  des_t <- if (length(t_idx))      data$xt[, t_idx,   drop = FALSE] else NULL

  fail <- list(ok = FALSE, loglik = -Inf, k = length(nm), par = NULL,
               l1 = Inf, hstr = hstr,
               t_idx = t_idx, h_idx = h_idx, aft_idx = aft_idx)

  ## -- unpenalised MLE -------------------------------------------------------
  fit <- try(GHMLE(init = unname(init), times = data$y, status = data$status,
                   hstr = hstr, dist = dist, des = des, des_t = des_t,
                   method = "nlminb", maxit = maxit),
             silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$OPT)) return(fail)

  nll  <- fit$log_lik                    # GHMLE minimises the NEGATIVE loglik
  par  <- fit$OPT$par
  obj  <- .opt_obj(fit$OPT)
  if (!is.finite(obj)) return(fail)

  ## -- L1-penalised refit, warm started at the unpenalised MLE ---------------
  ncoef <- length(nm) - nb
  if (lambda > 0 && ncoef > 0L) {
    ci      <- seq.int(nb + 1L, length(nm))
    pen_obj <- function(p) nll(p) + lambda * sum(.abs_smooth(p[ci]))
    opt     <- try(nlminb(par, pen_obj, control = list(iter.max = maxit)),
                   silent = TRUE)
    if (!inherits(opt, "try-error") && is.finite(.opt_obj(opt))) {
      par <- opt$par
      obj <- nll(par)
      if (!is.finite(obj)) return(fail)
    }
  }

  ## `k` is read straight off the fitted vector, so a mismatch between the
  ## names we expect and the parameters GHMLE actually fitted would silently
  ## corrupt every AIC/BIC.  Fail loudly instead.
  if (length(par) != length(nm))
    stop(sprintf("parameter/name length mismatch for hstr=%s (%d vs %d)",
                 hstr, length(par), length(nm)))

  names(par) <- nm
  list(ok      = TRUE,
       loglik  = -obj,
       k       = length(par),
       par     = par,
       l1      = if (ncoef > 0L) sum(abs(par[seq.int(nb + 1L, length(nm))])) else 0,
       hstr    = hstr,
       t_idx   = t_idx,
       h_idx   = h_idx,
       aft_idx = aft_idx)
}


## Score a fitted candidate on the chosen objective (smaller is better).
.crit_value <- function(f, crit, n_eff, lambda) {
  if (!isTRUE(f$ok) || !is.finite(f$loglik)) return(Inf)
  switch(crit,
         AIC   = -2 * f$loglik + 2 * f$k,
         BIC   = -2 * f$loglik + f$k * log(n_eff),
         LASSO = -2 * f$loglik + 2 * lambda * f$l1,
         stop("crit must be one of \"AIC\", \"BIC\", \"LASSO\"."))
}


## ==========================================================================
## 2.  Forward stepwise within one structure
## ==========================================================================
#
#   struct   "AH", "PH", "AFT" or "GH"
#
# Legal moves per step (each worth exactly one degree of freedom):
#   AH   0 -> 1
#   PH   0 -> 2
#   AFT  0 -> 4
#   GH   0 -> 1,  0 -> 2,  1 -> 3,  2 -> 3
#
# The search stops when no single move improves the score by more than `tol`,
# when `max.steps` is reached, or when `max.active` covariates are already in.
.forward_arm <- function(struct, data, p, crit, dist, lambda, maxit, n_eff,
                         max.steps, max.active, null_fit, tol, verbose) {

  gam       <- rep(0L, p)
  best_fit  <- null_fit
  best_sc   <- .crit_value(null_fit, crit, n_eff, lambda)
  n_fits    <- 0L
  steps     <- 0L

  ## Best genuinely GH-typed model seen anywhere in the search -- accepted or
  ## merely evaluated (see gh.strict in step.gh.freq).  Tracking evaluated
  ## candidates and not just the accepted path matters: on PH-like data the arm
  ## can walk from the null model to a PH-typed optimum without ever *accepting*
  ## a GH-typed state, and would otherwise report nothing at all.
  gh_fit <- NULL; gh_sc <- Inf; gh_gam <- NULL

  repeat {
    if (steps >= max.steps) break

    cand_sc  <- best_sc
    cand_fit <- NULL
    cand_gam <- NULL

    for (j in seq_len(p)) {
      moves <- switch(
        struct,
        AH  = if (gam[j] == 0L) 1L else integer(0),
        PH  = if (gam[j] == 0L) 2L else integer(0),
        AFT = if (gam[j] == 0L) 4L else integer(0),
        GH  = if (gam[j] == 0L) c(1L, 2L)
              else if (gam[j] == 1L || gam[j] == 2L) 3L
              else integer(0))

      for (v in moves) {
        g <- gam; g[j] <- v
        if (sum(g != 0L) > max.active) next

        f <- .fit_one(t_idx   = which(g == 1L | g == 3L),
                      h_idx   = which(g == 2L | g == 3L),
                      aft_idx = which(g == 4L),
                      data = data, dist = dist, lambda = lambda,
                      maxit = maxit, start = best_fit$par)
        n_fits <- n_fits + 1L

        sc <- .crit_value(f, crit, n_eff, lambda)
        if (!is.finite(sc)) next

        if (struct == "GH" && sc < gh_sc && gam2haz(g) == "GH") {
          gh_fit <- f; gh_sc <- sc; gh_gam <- g
        }
        if (sc < cand_sc - tol) {
          cand_sc <- sc; cand_fit <- f; cand_gam <- g
        }
      }
    }

    if (is.null(cand_fit)) break

    gam <- cand_gam; best_fit <- cand_fit; best_sc <- cand_sc
    steps <- steps + 1L

    if (verbose)
      cat(sprintf("  [%-3s] step %2d  %s  %s = %.3f\n",
                  struct, steps, paste0(gam, collapse = ""), crit, best_sc))
  }

  list(structure = struct,
       fit       = best_fit,
       gam       = gam,
       score     = best_sc,
       steps     = steps,
       n.fits    = n_fits,
       gh.fit    = gh_fit,          # NULL outside the GH arm, or if none seen
       gh.gam    = gh_gam,
       gh.score  = gh_sc)
}


## ==========================================================================
## 3.  Main entry point
## ==========================================================================
#
#   data        list with x, xt, y, status (as produced by simdat()).  If xt is
#               absent it defaults to x, matching the simulation set-up.
#   crit        "BIC" (default), "AIC" or "LASSO"
#   lambda      only used when crit == "LASSO":
#                 numeric  a fixed penalty on the sum-scale log-likelihood
#                 "cox"    (default) n * cv.glmnet(...)$lambda.1se from a Cox
#                          lasso on the same data -- cheap and data driven
#                 "cv"     K-fold CV over `lambda.grid`, scored by held-out
#                          log-likelihood.  Statistically the right thing but
#                          it reruns the whole four-arm search nfolds*|grid|
#                          times; only practical at small p.
#   dist        baseline distribution passed to GHMLE.  "LogNormal" matches the
#               baseline the Bayesian comparators assume.
#   n_eff       sample size in the BIC penalty; nrow(x) by default.  Set to
#               sum(status) to penalise by the number of events instead.
#   max.steps   cap on forward steps per arm
#   max.active  cap on the number of active covariates per arm
#   gh.strict   FALSE (default): the GH arm reports its unrestricted optimum,
#               which may itself be AH- or PH-typed because those are nested in
#               GH.  TRUE: the GH arm instead reports the best *genuinely
#               GH-typed* model it evaluated, which keeps the four-way
#               comparison a comparison of four distinct structures.  If the
#               arm never reached a GH-typed model its score is Inf and its row
#               of the table is NA.
#   standardize scale covariates to unit SD before fitting (default TRUE).  The
#               likelihood is invariant to this; it only makes the L1 penalty
#               scale free.  Coefficients are returned on both scales.
#
# Returns an object of class "step.gh.freq" (so print() gives the summary):
#   best        the selected model: structure, gam, gam.str, hstr, score,
#               loglik, k, AIC, BIC, coef (original scale), coef.std, par
#   arms        the four per-structure results, named AH/PH/AFT/GH
#   table       one row per structure, sorted by `score`, with AIC and BIC
#               shown alongside whichever criterion drove the search.  The l1
#               column is the L1 norm on the *standardised* scale, i.e. the
#               quantity the LASSO objective actually penalises.
#   null        the covariate-free baseline fit
#   crit, lambda, dist, n_eff, sd.x, sd.xt, n.fits, runtime
step.gh.freq <- function(data,
                         crit        = c("BIC", "AIC", "LASSO"),
                         lambda      = "cox",
                         dist        = "LogNormal",
                         maxit       = 5000,
                         n_eff       = NULL,
                         max.steps   = NULL,
                         max.active  = NULL,
                         gh.strict   = FALSE,
                         standardize = TRUE,
                         tol         = 1e-6,
                         verbose     = FALSE) {

  crit <- match.arg(crit)
  t0   <- Sys.time()

  if (!exists("GHMLE"))
    stop("GHMLE() not found -- source(\"./HazReg/routines.R\") first.")
  if (!exists("gam2haz"))
    stop("gam2haz() not found -- source(\"./utils/helper_functions.R\") first.")

  ## ---- normalise the data object -----------------------------------------
  if (is.null(data$y) && !is.null(data$time))  data$y <- data$time
  if (is.null(data$y) && !is.null(data$times)) data$y <- data$times
  X  <- as.matrix(data$x)
  Xt <- if (is.null(data$xt)) X else as.matrix(data$xt)
  if (ncol(X) != ncol(Xt))
    stop("x and xt must have the same number of columns.")

  n <- nrow(X)
  p <- ncol(X)
  status <- as.integer(data$status)

  ## ---- standardisation ----------------------------------------------------
  sd_x  <- if (standardize) apply(X,  2, stats::sd) else rep(1, p)
  sd_xt <- if (standardize) apply(Xt, 2, stats::sd) else rep(1, p)
  sd_x[!is.finite(sd_x)   | sd_x  <= 0] <- 1
  sd_xt[!is.finite(sd_xt) | sd_xt <= 0] <- 1
  dat <- list(x      = sweep(X,  2, sd_x,  "/"),
              xt     = sweep(Xt, 2, sd_xt, "/"),
              y      = as.vector(data$y),
              status = status)

  if (is.null(n_eff))      n_eff <- n
  if (is.null(max.steps))  max.steps <- min(2L * p, 40L)
  if (is.null(max.active)) max.active <- min(p, max(1L, floor(sum(status) / 2)))

  ## ---- penalty ------------------------------------------------------------
  lam <- 0
  if (crit == "LASSO") {
    lam <- if (is.numeric(lambda)) as.numeric(lambda)
           else if (identical(lambda, "cox")) lambda.cox(dat)
           else if (identical(lambda, "cv"))
             cv.lambda.gh(dat, dist = dist, maxit = maxit,
                          max.steps = max.steps, max.active = max.active,
                          gh.strict = gh.strict, verbose = verbose)$lambda
           else stop("lambda must be numeric, \"cox\" or \"cv\".")
    if (!is.finite(lam) || lam <= 0)
      warning("lambda <= 0: the LASSO objective reduces to the raw ",
              "log-likelihood and forward selection will not stop early.")
    if (verbose) cat(sprintf("LASSO penalty: lambda = %.4f\n", lam))
  }

  ## ---- null model ---------------------------------------------------------
  null_fit <- .fit_one(integer(0), integer(0), integer(0), dat, dist,
                       lambda = lam, maxit = maxit)
  if (!null_fit$ok)
    stop("The covariate-free baseline model failed to fit; check `dist`.")

  ## ---- Stage 1: one forward search per structure --------------------------
  arms <- list()
  for (s in c("AH", "PH", "AFT", "GH")) {
    if (s == "GH" && dist == "Weibull") next
    if (verbose) cat(sprintf("Structure %s\n", s))
    arms[[s]] <- .forward_arm(s, dat, p, crit, dist, lam, maxit, n_eff,
                              max.steps, max.active, null_fit, tol, verbose)
  }

  ## The GH arm nests AH and PH, so its optimum need not be GH-typed.  Under
  ## gh.strict the arm instead contributes the best GH-typed model it visited.
  if (gh.strict && !is.null(arms$GH)) {
    if (is.null(arms$GH$gh.fit)) {
      arms$GH$fit <- NULL; arms$GH$score <- Inf; arms$GH$gam <- rep(0L, p)
    } else {
      arms$GH$fit <- arms$GH$gh.fit
      arms$GH$gam <- arms$GH$gh.gam
      arms$GH$score <- arms$GH$gh.score
    }
  }

  ## ---- Stage 2: compare the four winners ----------------------------------
  ## Rescale coefficients back to the original covariate scale.
  ## Coefficients only (drop the baseline parameters), on the original scale.
  ## x enters solely through x'b, so b_orig = b_std / sd_j.
  coefs_of <- function(f) {
    if (is.null(f) || is.null(f$par)) return(NULL)
    f$par[!grepl("^base", names(f$par))]
  }
  unstd <- function(f) {
    cf <- coefs_of(f)
    if (is.null(cf) || !length(cf)) return(cf)
    lvl <- substr(names(cf), 1, 1)
    idx <- as.integer(substring(names(cf), 2))
    cf / ifelse(lvl == "t", sd_xt[idx], sd_x[idx])
  }

  rows <- lapply(names(arms), function(s) {
    a <- arms[[s]]
    f <- a$fit
    if (is.null(f) || !isTRUE(f$ok)) {
      return(data.frame(structure = s, gam.str = NA_character_, hstr = NA_character_,
                        n.active = NA_integer_, k = NA_integer_, loglik = NA_real_,
                        l1 = NA_real_, score = a$score, AIC = NA_real_,
                        BIC = NA_real_, steps = a$steps, n.fits = a$n.fits,
                        stringsAsFactors = FALSE))
    }
    data.frame(structure = s,
               gam.str   = paste0(a$gam, collapse = ""),
               hstr      = gam2haz(a$gam),
               n.active  = sum(a$gam != 0L),
               k         = f$k,
               loglik    = f$loglik,
               l1        = f$l1,
               score     = a$score,
               AIC       = -2 * f$loglik + 2 * f$k,
               BIC       = -2 * f$loglik + f$k * log(n_eff),
               steps     = a$steps,
               n.fits    = a$n.fits,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$score), ]
  rownames(tab) <- NULL

  win <- arms[[tab$structure[1]]]
  wf  <- win$fit

  best <- list(structure = win$structure,
               gam       = as.integer(win$gam),
               gam.str   = paste0(win$gam, collapse = ""),
               hstr      = gam2haz(win$gam),
               score     = win$score,
               loglik    = wf$loglik,
               k         = wf$k,
               AIC       = -2 * wf$loglik + 2 * wf$k,
               BIC       = -2 * wf$loglik + wf$k * log(n_eff),
               par       = wf$par,          # standardised scale, incl. baseline
               coef      = unstd(wf),       # original covariate scale
               coef.std  = coefs_of(wf))

  out <- list(
    best    = best,
    arms    = arms,
    table   = tab,
    null    = list(loglik = null_fit$loglik, k = null_fit$k,
                   score  = .crit_value(null_fit, crit, n_eff, lam)),
    crit    = crit,
    lambda  = lam,
    dist    = dist,
    n_eff   = n_eff,
    sd.x    = sd_x,
    sd.xt   = sd_xt,
    n.fits  = sum(vapply(arms, function(a) a$n.fits, integer(1))) + 1L,
    runtime = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  class(out) <- "step.gh.freq"
  out
}


## ==========================================================================
## 4.  Choosing lambda
## ==========================================================================

## Cheap data-driven default: a Cox lasso on the same design.  glmnet minimises
## -(1/n)*partial loglik + lambda*||b||_1, whereas our objective is on the sum
## scale, hence the factor n.
lambda.cox <- function(data, s = c("lambda.1se", "lambda.min"), nfolds = 10) {
  s <- match.arg(s)
  if (!requireNamespace("glmnet", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE))
    stop("lambda = \"cox\" needs the glmnet and survival packages.")
  cvfit <- glmnet::cv.glmnet(x = as.matrix(data$x),
                             y = survival::Surv(data$y, data$status),
                             family = "cox", alpha = 1,
                             nfolds = nfolds, maxit = 10000)
  as.numeric(nrow(data$x) * cvfit[[s]])
}

## Negative log-likelihood of a given active set, as a function of the
## parameters, on an arbitrary data set.  Built by asking GHMLE() for one
## iteration and keeping the closure it returns -- this reuses the likelihood
## code in HazReg/routines.R rather than restating it.
.nll_fun <- function(t_idx, h_idx, aft_idx, data, dist) {
  hstr  <- .hstr_of(t_idx, h_idx, aft_idx)
  des   <- if (hstr == "AFT")      data$x[,  aft_idx, drop = FALSE]
           else if (length(h_idx)) data$x[,  h_idx,   drop = FALSE]
           else                    NULL
  des_t <- if (length(t_idx))      data$xt[, t_idx,   drop = FALSE] else NULL
  fit <- try(GHMLE(init = c(.start_base(data$y, dist),
                            rep(0, length(t_idx) + length(h_idx) + length(aft_idx))),
                   times = data$y, status = data$status, hstr = hstr,
                   dist = dist, des = des, des_t = des_t,
                   method = "nlminb", maxit = 1L), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit$log_lik
}

## K-fold CV for lambda.  For each fold and each candidate lambda the whole
## four-arm search is rerun on the training folds and the selected model is
## scored by its negative log-likelihood on the held-out fold.  Cost is
## nfolds * length(lambda.grid) full searches -- only practical at small p.
cv.lambda.gh <- function(data, lambda.grid = NULL, nfolds = 5L,
                         dist = "LogNormal", maxit = 5000, max.steps = NULL,
                         max.active = NULL, gh.strict = FALSE, seed = 1L,
                         verbose = FALSE) {

  n <- nrow(as.matrix(data$x))
  if (is.null(lambda.grid)) {
    hi <- lambda.cox(data, s = "lambda.min")
    if (!is.finite(hi) || hi <= 0) hi <- n / 10
    lambda.grid <- exp(seq(log(hi * 4), log(hi / 8), length.out = 8L))
  }

  set.seed(seed)
  folds <- sample(rep_len(seq_len(nfolds), n))
  cvnll <- rep(0, length(lambda.grid))

  for (k in seq_len(nfolds)) {
    tr <- which(folds != k)
    te <- which(folds == k)
    dtr <- list(x = data$x[tr, , drop = FALSE], xt = data$xt[tr, , drop = FALSE],
                y = data$y[tr], status = data$status[tr])
    dte <- list(x = data$x[te, , drop = FALSE], xt = data$xt[te, , drop = FALSE],
                y = data$y[te], status = data$status[te])

    for (l in seq_along(lambda.grid)) {
      f <- try(step.gh.freq(dtr, crit = "LASSO", lambda = lambda.grid[l],
                            dist = dist, maxit = maxit, max.steps = max.steps,
                            max.active = max.active, gh.strict = gh.strict,
                            standardize = FALSE, verbose = FALSE),
               silent = TRUE)
      if (inherits(f, "try-error")) { cvnll[l] <- cvnll[l] + Inf; next }

      wf  <- f$arms[[f$best$structure]]$fit
      nll <- .nll_fun(wf$t_idx, wf$h_idx, wf$aft_idx, dte, dist)
      v   <- if (is.null(nll)) Inf else try(nll(wf$par), silent = TRUE)
      cvnll[l] <- cvnll[l] + if (inherits(v, "try-error") || !is.finite(v)) Inf
                             else as.numeric(v)
    }
    if (verbose) cat(sprintf("  CV fold %d/%d done\n", k, nfolds))
  }

  list(lambda = lambda.grid[which.min(cvnll)],
       grid   = lambda.grid,
       cv.nll = cvnll)
}


## ==========================================================================
## 5.  Reporting
## ==========================================================================

print.step.gh.freq <- function(x, ...) {
  cat(sprintf("Forward stepwise per hazard structure  (criterion: %s%s, baseline: %s)\n",
              x$crit,
              if (x$crit == "LASSO") sprintf(", lambda = %.4f", x$lambda) else "",
              x$dist))
  cat(sprintf("%d model fits in %.1f s\n\n", x$n.fits, x$runtime))
  print(x$table, digits = 5)
  cat(sprintf("\nSelected: %s structure, model %s (%s), %s = %.3f\n",
              x$best$structure, x$best$gam.str, x$best$hstr,
              x$crit, x$best$score))
  invisible(x)
}
