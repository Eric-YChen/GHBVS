###############################################################################
# utils/freq_stepwise_bic.R
#
# Frequentist forward-stepwise selection for the general hazard (GH) family,
# scored by BIC. Selects the active covariate set AND the hazard structure
# jointly, over all four structures gam2haz() recognises: AH, PH, AFT and GH.
#
# Why this exists
# ---------------
# The earlier comparator reused init.gam(), which ends in bools2gam(). That
# maps a 2-bit (time-level, hazard-level) indicator onto codes {0,1,2,3}, so
# code 4 (AFT) is unreachable by construction. On AFT scenarios -- 40% of the
# BIC comparison grid -- it therefore scored exactly 0 hazard accuracy no
# matter what the data said. This routine can reach every structure.
#
# Encoding (as elsewhere in the project)
#   0 = excluded
#   1 = time-level effect only      (AH)
#   2 = hazard-level effect only    (PH)
#   3 = both time and hazard level  (GH)
#   4 = AFT
#
# Search
#   Arm A  codes {0,1,2,3}. Every step adds exactly one parameter, via
#          0->1, 0->2, 1->3 or 2->3, so successive BIC values differ by one
#          degree of freedom. The arm starts at the null model, passes through
#          AH / PH states and continues into GH.
#   Arm B  codes {0,4}. Every step adds one AFT covariate. AFT is searched on
#          its own arm because gam2haz() treats a model mixing code 4 with
#          1/2/3 as "Invalid model".
#   The returned model is the global BIC minimum over the two arms.
#
# Scoring
#   BIC = -2*loglik + k*log(n_eff),  k = length(fitted par)
#                                      = 2 baseline params + coefficients
#   This matches marg.loglik.bic() in utils/bic_functions.R with n_type = "n"
#   (k = d_gamma + 2), so the frequentist and Bayesian BIC arms are on the
#   same scale.
#
# Dependencies: HazReg/routines.R (GHMLE), utils/helper_functions.R (gam2haz)
###############################################################################

## Moment-based starting values for the LogNormal baseline: par[1] is the
## location on the log-time scale and par[2] is log(scale). Starting both at 0
## (as the older routines do) is a poor guess when survival times are large.
.lnorm_start <- function(times, k) {
  lt <- log(pmax(as.vector(times), 1e-8))
  s  <- stats::sd(lt)
  if (!is.finite(s) || s <= 0) s <- 1
  c(mean(lt), log(s), rep(0, k))
}

## Fit one model and return its BIC.
##   t_idx   columns with a time-level effect   (des_t)  -- codes 1 and 3
##   h_idx   columns with a hazard-level effect (des)    -- codes 2 and 3
##   aft_idx columns with an AFT effect         (des)    -- code 4
## Exactly one of {aft_idx} and {t_idx, h_idx} may be non-empty.
## Returns Inf if the optimiser fails or the objective is not finite, so a
## failed candidate can never win the argmin.
.bic_fit <- function(t_idx, h_idx, aft_idx, data, n_eff, dist, maxit) {
  X  <- data$x
  Xt <- data$xt
  n_t <- length(t_idx); n_h <- length(h_idx); n_a <- length(aft_idx)

  fit <- try({
    if (n_a > 0L) {
      GHMLE(init = .lnorm_start(data$y, n_a),
            times = data$y, status = data$status, hstr = "AFT", dist = dist,
            des = X[, aft_idx, drop = FALSE],
            method = "nlminb", maxit = maxit)
    } else if (n_t > 0L && n_h > 0L) {
      ## GH: the LogNormal branch reads des_t coefficients first, then des
      GHMLE(init = .lnorm_start(data$y, n_t + n_h),
            times = data$y, status = data$status, hstr = "GH", dist = dist,
            des   = X[,  h_idx, drop = FALSE],
            des_t = Xt[, t_idx, drop = FALSE],
            method = "nlminb", maxit = maxit)
    } else if (n_t > 0L) {
      GHMLE(init = .lnorm_start(data$y, n_t),
            times = data$y, status = data$status, hstr = "AH", dist = dist,
            des_t = Xt[, t_idx, drop = FALSE],
            method = "nlminb", maxit = maxit)
    } else if (n_h > 0L) {
      GHMLE(init = .lnorm_start(data$y, n_h),
            times = data$y, status = data$status, hstr = "PH", dist = dist,
            des = X[, h_idx, drop = FALSE],
            method = "nlminb", maxit = maxit)
    } else {
      GHMLE(init = .lnorm_start(data$y, 0L),
            times = data$y, status = data$status, hstr = "baseline",
            dist = dist, method = "nlminb", maxit = maxit)
    }
  }, silent = TRUE)

  if (inherits(fit, "try-error") || is.null(fit$OPT)) return(Inf)
  obj <- fit$OPT$objective            # GHMLE minimises the NEGATIVE log-lik
  if (!is.finite(obj)) return(Inf)
  2 * obj + length(fit$OPT$par) * log(n_eff)
}

## Convenience: BIC of a gamma vector under this encoding.
bic.of.gam <- function(gam, data, n_eff = nrow(data$x),
                       dist = "LogNormal", maxit = 5000) {
  gam <- as.integer(gam)
  if (any(gam == 4L))
    return(.bic_fit(integer(0), integer(0), which(gam == 4L),
                    data, n_eff, dist, maxit))
  .bic_fit(which(gam == 1L | gam == 3L), which(gam == 2L | gam == 3L),
           integer(0), data, n_eff, dist, maxit)
}

###############################################################################
# step.bic.gh
#
#   data       list with x, xt, y, status (as produced by simdat)
#   dist       baseline distribution passed to GHMLE. "LogNormal" matches the
#              baseline the Bayesian methods assume (hessLLS(..., dist="lnorm"))
#   n_eff      sample size in the BIC penalty; nrow(data$x) by default
#   max.steps  cap on forward steps per arm, guarding against runaway searches
#              at large p. BIC's log(n) penalty normally halts well before this.
#   n.threshold optional cap on model size, mirroring prior.gam()'s constraint
#
# Returns a list:
#   gam         selected gamma vector (integer codes)
#   gam.str     the same as a model string
#   hstr        gam2haz(gam)
#   bic         BIC of the selected model
#   bic.armA    best BIC over codes {0,1,2,3}
#   bic.armB    best BIC over codes {0,4}
#   n.fits      number of GHMLE fits performed (for runtime accounting)
#   steps       steps taken on each arm
###############################################################################
step.bic.gh <- function(data, dist = "LogNormal", maxit = 5000,
                        n_eff = NULL, max.steps = NULL, n.threshold = NULL) {

  X <- as.matrix(data$x)
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(n_eff))     n_eff <- n
  if (is.null(max.steps)) max.steps <- min(2L * p, 40L)
  if (is.null(n.threshold))
    n.threshold <- n * (1 - 0.5 * (1 - mean(data$status)))

  n_fits <- 0L
  score <- function(t_idx, h_idx, aft_idx) {
    n_fits <<- n_fits + 1L
    .bic_fit(t_idx, h_idx, aft_idx, data, n_eff, dist, maxit)
  }

  bic_null <- score(integer(0), integer(0), integer(0))

  ## ---- Arm A: codes {0,1,2,3} -- AH / PH / GH ----------------------------
  gamA <- rep(0L, p)
  bicA <- bic_null
  stepsA <- 0L
  repeat {
    if (stepsA >= max.steps) break
    best_bic <- bicA
    best_gam <- NULL
    for (j in seq_len(p)) {
      moves <- if (gamA[j] == 0L) c(1L, 2L) else
               if (gamA[j] == 1L || gamA[j] == 2L) 3L else integer(0)
      for (v in moves) {
        g <- gamA; g[j] <- v
        if (sum(g != 0L) > n.threshold) next
        b <- score(which(g == 1L | g == 3L), which(g == 2L | g == 3L),
                   integer(0))
        if (is.finite(b) && b < best_bic - 1e-8) {
          best_bic <- b; best_gam <- g
        }
      }
    }
    if (is.null(best_gam)) break
    gamA <- best_gam; bicA <- best_bic; stepsA <- stepsA + 1L
  }

  ## ---- Arm B: codes {0,4} -- AFT -----------------------------------------
  gamB <- rep(0L, p)
  bicB <- bic_null
  stepsB <- 0L
  repeat {
    if (stepsB >= max.steps) break
    best_bic <- bicB
    best_gam <- NULL
    for (j in which(gamB == 0L)) {
      g <- gamB; g[j] <- 4L
      if (sum(g != 0L) > n.threshold) next
      b <- score(integer(0), integer(0), which(g == 4L))
      if (is.finite(b) && b < best_bic - 1e-8) {
        best_bic <- b; best_gam <- g
      }
    }
    if (is.null(best_gam)) break
    gamB <- best_gam; bicB <- best_bic; stepsB <- stepsB + 1L
  }

  ## ---- Global minimum over the two arms ----------------------------------
  if (is.finite(bicB) && bicB < bicA) {
    gam_sel <- gamB; bic_sel <- bicB
  } else {
    gam_sel <- gamA; bic_sel <- bicA
  }

  list(gam      = as.integer(gam_sel),
       gam.str  = paste0(as.integer(gam_sel), collapse = ""),
       hstr     = gam2haz(gam_sel),
       bic      = bic_sel,
       bic.armA = bicA,
       bic.armB = bicB,
       bic.null = bic_null,
       n.fits   = n_fits,
       steps    = c(armA = stepsA, armB = stepsB))
}
