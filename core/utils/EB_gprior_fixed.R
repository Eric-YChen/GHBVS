if (!exists("BVS_CORE")) BVS_CORE <- "."
source(file.path(BVS_CORE, "utils", "mle_gh.R"))

marg.loglik <- function(par, data, max.iter = 200, return.terms = FALSE) {
  par.orig <- mle.gh(par, data, max.iter = max.iter)
  par <- par.gamma(par.orig)
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  g <- par$g; gam <- par$gam

  theta.keep <- as.logical(gam == 1 | gam == 3 | gam == 4)
  eta.keep <- as.logical(gam == 2 | gam == 3)
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

  if (pt + p > 0) {
    J_theta_eta <- hess_full[idx_beta, idx_beta]
    J_theta_eta_nu <- hess_full[idx_beta, idx_nu]
    J_theta_eta_theta0 <- hess_full[idx_beta, idx_theta0]

    J_theta_eta_inv <- solve(J_theta_eta)
    factor <- g / (1 + g)

    J_tilde_nu_nu <- J_nu_nu - factor * t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_nu
    J_tilde_nu_theta0 <- J_nu_theta0 - factor * t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_theta0
    J_tilde_theta0_theta0 <- J_theta0_theta0 - factor * t(J_theta_eta_theta0) %*% J_theta_eta_inv %*% J_theta_eta_theta0

    khat_gamma   <- c(theta, eta)
    khat_Jnu     <- as.numeric(t(khat_gamma) %*% J_theta_eta_nu)
    khat_Jtheta0 <- as.numeric(t(khat_gamma) %*% J_theta_eta_theta0)
    W_gamma      <- as.numeric(t(khat_gamma) %*% J_theta_eta %*% khat_gamma)
  } else {
    J_tilde_nu_nu <- J_nu_nu
    J_tilde_nu_theta0 <- J_nu_theta0
    J_tilde_theta0_theta0 <- J_theta0_theta0
  }

  mu_nu <- 9.34
  sigma2_nu <- 41.15^2
  K <- 1e6

  P <- matrix(0, 2, 2)
  P[1,1] <- J_tilde_nu_nu + 1/sigma2_nu
  P[1,2] <- P[2,1] <- J_tilde_nu_theta0
  P[2,2] <- J_tilde_theta0_theta0 + 1/K

  h <- c(
    J_tilde_nu_nu * nu + J_tilde_nu_theta0 * theta0 + mu_nu/sigma2_nu,
    J_tilde_nu_theta0 * nu + J_tilde_theta0_theta0 * theta0
  )
  if (pt + p > 0) {
    h[1] <- h[1] + khat_Jnu / (1 + g)
    h[2] <- h[2] + khat_Jtheta0 / (1 + g)
  }

  C <- -mu_nu^2 / (2 * sigma2_nu) - 0.5 * nu^2 * J_nu_nu - nu * theta0 * J_nu_theta0 - 0.5 * theta0^2 * J_theta0_theta0
  if (pt + p > 0) {
    factor_half <- g / (2 * (1 + g))
    factor_full <- g / (1 + g)
    C <- C + factor_half * nu^2 * (t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_nu)
    C <- C + factor_full * nu * theta0 * (t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_theta0)
    C <- C + factor_half * theta0^2 * (t(J_theta_eta_theta0) %*% J_theta_eta_inv %*% J_theta_eta_theta0)
    C <- C - W_gamma / (2 * (1 + g))
    C <- C - (khat_Jnu * nu + khat_Jtheta0 * theta0) / (1 + g)
  }

  log_det_P <- determinant(P, logarithm = TRUE)$modulus[1]
  P_inv <- solve(P)

  d_gam        <- pt + p
  term_ell     <- as.numeric(log_lik_map)
  term_dim     <- -0.5 * d_gam * log(1 + g)
  term_logdetP <- -0.5 * log_det_P
  term_quad    <-  0.5 * as.numeric(t(h) %*% P_inv %*% h)
  term_const   <-  as.numeric(C)

  marg_log_lik <- term_ell + term_dim + term_logdetP + term_quad + term_const

  if (isTRUE(return.terms)) {
    n_obs  <- nrow(as.matrix(data$x))
    n_ev   <- sum(as.numeric(data$status) == 1)
    k_gam  <- d_gam + 2
    bic_pen_n  <- -0.5 * k_gam * log(n_obs)
    bic_pen_ev <- -0.5 * k_gam * log(max(n_ev, 1))
    return(list(
      gam          = paste0(par.orig$gam, collapse = ""),
      d_gam        = d_gam,
      g            = g,
      n_obs        = n_obs,
      n_events     = n_ev,
      ell          = term_ell,
      ila_dim      = term_dim,
      ila_logdetP  = term_logdetP,
      ila_quad     = term_quad,
      ila_const    = term_const,
      ila_extra    = term_dim + term_logdetP + term_quad + term_const,
      ila_total    = as.numeric(marg_log_lik),
      bic_pen_n    = bic_pen_n,
      bic_pen_ev   = bic_pen_ev,
      bic_total_n  = term_ell + bic_pen_n,
      bic_total_ev = term_ell + bic_pen_ev
    ))
  }

  return(as.numeric(marg_log_lik))
}

targetdist = function(par, data, max.iter = 200, ...){
  if (sum(par$gam)==0 && !(exists("admit.null") && isTRUE(admit.null))){
    return(-Inf)
  }
  lp <- prior.gam(par, log = TRUE)
  if (!is.finite(lp)) return(-Inf)
  log.approx.marg.llh <- marg.loglik(par, data, max.iter = max.iter)
  return(lp + log.approx.marg.llh)
}
