if (!exists("BVS_CORE")) BVS_CORE <- "."
source(file.path(BVS_CORE, "utils", "mle_gh.R"))

marg.loglik <- function(par, data, max.iter = 200) {
  par.orig <- mle.gh(par, data, max.iter = max.iter)
  par <- par.gamma(par.orig)
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  gam <- par$gam

  theta.keep = as.logical(par$gam == 1 | par$gam == 3 | par$gam == 4)
  eta.keep = as.logical(par$gam == 2 | par$gam == 3)
  pt <- sum(theta.keep)
  p <- sum(eta.keep)

  log_lik_map <- lglik(par.orig, data)
  hess_full <- -hessLLS(par.orig, data)
  keep.all <- c(TRUE, TRUE, theta.keep, eta.keep)
  if (!is.finite(log_lik_map) || !is.pd(hess_full[keep.all, keep.all, drop = FALSE])) {
    log.opt.problem(gam, "observed information not positive definite at the MLE")
    return(-Inf)
  }

  idx_nu <- 1
  idx_theta0 <- 2
  idx_theta <- which(c(FALSE, FALSE, theta.keep))
  idx_eta <- which(c(FALSE, FALSE, rep(FALSE, length(theta.keep)), eta.keep))
  idx_beta <- c(idx_theta, idx_eta)

  J_nu_nu <- hess_full[idx_nu, idx_nu]
  J_nu_theta0 <- hess_full[idx_nu, idx_theta0]
  J_theta0_theta0 <- hess_full[idx_theta0, idx_theta0]

  has_beta <- (pt + p) > 0
  if (has_beta) {
    J_theta_eta <- hess_full[idx_beta, idx_beta]
    J_theta_eta_nu <- hess_full[idx_beta, idx_nu]
    J_theta_eta_theta0 <- hess_full[idx_beta, idx_theta0]
    J_theta_eta_inv <- solve(J_theta_eta)
    quad_nu <- as.numeric(t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_nu)
    quad_cross <- as.numeric(t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_theta0)
    quad_theta0 <- as.numeric(t(J_theta_eta_theta0) %*% J_theta_eta_inv %*% J_theta_eta_theta0)

    khat_gamma   <- c(theta, eta)
    khat_Jnu     <- as.numeric(t(khat_gamma) %*% J_theta_eta_nu)
    khat_Jtheta0 <- as.numeric(t(khat_gamma) %*% J_theta_eta_theta0)
    W_gamma      <- as.numeric(t(khat_gamma) %*% J_theta_eta %*% khat_gamma)
  } else {
    quad_nu <- 0
    quad_cross <- 0
    quad_theta0 <- 0
    khat_Jnu <- 0
    khat_Jtheta0 <- 0
    W_gamma <- 0
  }

  mu_nu <- 9.34
  sigma2_nu <- 41.15^2
  K <- 1e6

  base_C <- -mu_nu^2 / (2 * sigma2_nu) - 0.5 * nu^2 * J_nu_nu - nu * theta0 * J_nu_theta0 - 0.5 * theta0^2 * J_theta0_theta0

  k <- pt + p
  n <- par$n

  if (is.null(par$a_g) || par$a_g <= 2 || par$a_g > 4) {
    a <- (1 + n) / (k + 1) - 1
    use_hyper_g <- FALSE
  } else {
    a_g <- par$a_g
    a <- 0
    use_hyper_g <- TRUE
  }

  log_integrand_t <- function(t) {
    g_val <- a + exp(t)
    one_plus_g <- 1 + g_val
    factor_full <- g_val / one_plus_g
    factor_half <- 0.5 * factor_full

    if (has_beta) {
      J_tilde_nu_nu <- J_nu_nu - factor_full * quad_nu
      J_tilde_nu_theta0 <- J_nu_theta0 - factor_full * quad_cross
      J_tilde_theta0_theta0 <- J_theta0_theta0 - factor_full * quad_theta0
    } else {
      len <- length(t)
      J_tilde_nu_nu <- rep(J_nu_nu, len)
      J_tilde_nu_theta0 <- rep(J_nu_theta0, len)
      J_tilde_theta0_theta0 <- rep(J_theta0_theta0, len)
    }

    P11 <- J_tilde_nu_nu + 1 / sigma2_nu
    P12 <- J_tilde_nu_theta0
    P22 <- J_tilde_theta0_theta0 + 1 / K
    det_P <- P11 * P22 - P12 * P12
    valid <- is.finite(det_P) & det_P > 0 & is.finite(one_plus_g) & (one_plus_g > 0)
    res <- rep(-Inf, length(t))
    if (!any(valid)) return(res)

    log_det_P <- log(det_P[valid])

    h1 <- J_tilde_nu_nu[valid] * nu + J_tilde_nu_theta0[valid] * theta0 + mu_nu / sigma2_nu
    h2 <- J_tilde_nu_theta0[valid] * nu + J_tilde_theta0_theta0[valid] * theta0
    if (has_beta) {
      h1 <- h1 + khat_Jnu / one_plus_g[valid]
      h2 <- h2 + khat_Jtheta0 / one_plus_g[valid]
    }
    invP_h1 <- ( P22[valid] * h1 - P12[valid] * h2) / det_P[valid]
    invP_h2 <- (-P12[valid] * h1 + P11[valid] * h2) / det_P[valid]
    quad_h <- 0.5 * (h1 * invP_h1 + h2 * invP_h2)

    Cg <- rep(base_C, sum(valid))
    if (has_beta) {
      Cg <- Cg + factor_half[valid] * (nu^2 * quad_nu + theta0^2 * quad_theta0) + factor_full[valid] * (nu * theta0 * quad_cross)
      Cg <- Cg - W_gamma / (2 * one_plus_g[valid])
      Cg <- Cg - (khat_Jnu * nu + khat_Jtheta0 * theta0) / one_plus_g[valid]
    }

    if (use_hyper_g) {
      log_prior_g <- log(a_g - 2) - log(2 * n) - (a_g / 2) * log(1 + g_val[valid] / n)
    } else {
      log_prior_g <- log(0.5) + 0.5 * log((1 + n) / (k + 1)) - 1.5 * log(one_plus_g[valid])
    }
    log_jac <- t[valid]

    out_valid <- log_lik_map - 0.5 * k * log(one_plus_g[valid]) - 0.5 * log_det_P + quad_h + Cg + log_prior_g + log_jac
    res[valid] <- as.numeric(out_valid)
    res
  }

  grid_t <- seq(-10, 15, length.out = 101)
  vals <- vapply(grid_t, log_integrand_t, numeric(1))
  m_shift <- max(vals[is.finite(vals)])
  if (!is.finite(m_shift)) m_shift <- 0

  f_scaled <- function(t) {
    vals <- log_integrand_t(t) - m_shift
    vals[!is.finite(vals)] <- -Inf
    exp(vals)
  }

  int_res <- try(integrate(f_scaled, lower = -Inf, upper = Inf, rel.tol = 1e-6, abs.tol = 0), silent = TRUE)
  if (inherits(int_res, "try-error")) {
    int_res <- integrate(f_scaled, lower = -20, upper = 20, rel.tol = 1e-6, abs.tol = 0)
  }
  log_int <- log(int_res$value) + m_shift

  return(as.numeric(log_int))
}

targetdist = function(par, data, ...){
  if (sum(par$gam)==0 && !(exists("admit.null") && isTRUE(admit.null))){
    return(-Inf)
  }
  lp <- prior.gam(par, log = TRUE)
  if (!is.finite(lp)) return(-Inf)
  log.approx.marg.llh <- marg.loglik(par, data)
  return(lp + log.approx.marg.llh)
}
