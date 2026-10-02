###########################################################################################
# Marginal Likelihood Calculation with Normal Prior on nu = log(tau)
###########################################################################################
neg.post <- function(par.vec, otherpar, data2,...){
  pars = otherpar
  pxt = ncol(data2$xt); px = ncol(data2$x)
  theta.keep = as.logical(pars$gam == 1 | pars$gam == 3 | pars$gam == 4)
  eta.keep = as.logical(pars$gam == 2 | pars$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  pars$nu = par.vec[1]; pars$theta0 = par.vec[2]
  if (ptheta == 0){
    pars$theta = rep(0,pxt)
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2+ptheta]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2+ptheta]
    pars$theta = thetas; pars$eta = etas
  }
  return(-lglik(pars,data2))
}

grad.post <- function(par.vec, otherpar, data2, ...){
  pars = otherpar
  pxt = ncol(data2$xt); px = ncol(data2$x)
  theta.keep = as.logical(pars$gam == 1 | pars$gam == 3 | pars$gam == 4)
  eta.keep = as.logical(pars$gam == 2 | pars$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  pars$nu = par.vec[1]; pars$theta0 = par.vec[2]
  if (ptheta == 0){
    pars$theta = rep(0,pxt)
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2+ptheta]
    pars$theta = thetas; pars$eta = etas
  }
  out.keep = c(T,T,theta.keep,eta.keep)
  out <- - gradLLS(pars,data2)[out.keep]
  return(out)
}

hess.post <- function(par.vec, otherpar, data2,...){
  pars = otherpar
  pxt = ncol(data2$xt); px = ncol(data2$x)
  theta.keep = as.logical(pars$gam == 1 | pars$gam == 3 | pars$gam == 4)
  eta.keep = as.logical(pars$gam == 2 | pars$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  pars$nu = par.vec[1]; pars$theta0 = par.vec[2]
  if (ptheta == 0){
    pars$theta = rep(0,pxt)
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2+ptheta]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[seq_len(ptheta)+2]
    etas = eta.keep
    etas[etas==1] = par.vec[seq_len(peta)+2+ptheta]
    pars$theta = thetas; pars$eta = etas
  }
  cols.keep = c(T,T,theta.keep,eta.keep)
  rows.keep = cols.keep
  out <- -hessLLS(pars,data2)[rows.keep,cols.keep]
  return(out)
}

mle.gh <- function(par1, data1, max.iter = 100){
  pxt = ncol(data1$xt); px = ncol(data1$x)
  theta.keep = as.logical(par1$gam == 1 | par1$gam == 3 | par1$gam == 4)
  eta.keep = as.logical(par1$gam == 2 | par1$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  D = 2 + ptheta + peta
  
  # Calculate MAP with nlminb
  result <- nlminb(start=rep(0,D), otherpar=par1, data2=data1, 
                   control = list(iter.max = max.iter),
                   objective=neg.post, gradient=grad.post, hessian=hess.post)
  if (!exists(quote(itervec))){
    itervec <<- c(result$iterations)
  }else{
    itervec <<- c(itervec,result$iterations)  
  }
  MAP = result$par
  MAP0 = list(nu=MAP[1], theta0 = MAP[2], theta = theta.keep, eta = eta.keep)
  if (ptheta == 0){
    MAP0$theta = rep(0,pxt)
    MAP0$eta[eta.keep==1] = MAP[seq_len(peta)+2]
  }else if (peta == 0){
    MAP0$eta = rep(0,px)
    MAP0$theta[theta.keep==1] = MAP[seq_len(ptheta)+2]
  }else{
    MAP0$theta[theta.keep==1] = MAP[seq_len(ptheta)+2]
    MAP0$eta[eta.keep==1] = MAP[seq_len(peta)+2+ptheta]
  }
  MAP = MAP0
  MAP = utils::modifyList(par1, MAP)
  return(MAP)
}

marg.loglik <- function(par, data, max.iter = 100, return.terms = FALSE) {
  # Extract parameters and data
  par <- mle.gh(par, data, max.iter = max.iter)
  par.orig <- par
  par <- par.gamma(par)
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  g <- par$g; gam <- par$gam
  X <- as.matrix(data$x); Xt <- as.matrix(data$xt)
  
  theta.keep = as.logical(par$gam == 1 | par$gam == 3 | par$gam == 4)
  eta.keep = as.logical(par$gam == 2 | par$gam == 3)
  pt <- sum(theta.keep)
  p <- sum(eta.keep)
  
  
  # Calculate log-likelihood at MAP
  log_lik_map <- lglik(par.orig, data)
  
  # Calculate the Hessian at MAP
  hess_full <- -hessLLS(par.orig, data)
  
  # Extract relevant blocks from the Hessian
  # Order is: nu, theta0, theta, eta
  idx_nu <- 1
  idx_theta0 <- 2
  idx_theta <- which(c(FALSE, FALSE, theta.keep))
  idx_eta <- which(c(FALSE, FALSE, rep(FALSE, length(theta.keep)), eta.keep))
  idx_beta <- c(idx_theta, idx_eta)  # Combined theta and eta indices
  
  # Extract Hessian blocks
  J_nu_nu <- hess_full[idx_nu, idx_nu]
  J_nu_theta0 <- hess_full[idx_nu, idx_theta0]
  J_theta0_theta0 <- hess_full[idx_theta0, idx_theta0]
  
  if (pt + p > 0) {
    J_theta_eta <- hess_full[idx_beta, idx_beta]
    J_theta_eta_nu <- hess_full[idx_beta, idx_nu]
    J_theta_eta_theta0 <- hess_full[idx_beta, idx_theta0]
    
    # Calculate adjusted Hessian components
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
    # No theta or eta parameters
    J_tilde_nu_nu <- J_nu_nu
    J_tilde_nu_theta0 <- J_nu_theta0
    J_tilde_theta0_theta0 <- J_theta0_theta0
  }
  
  # Prior parameters
  mu_nu <- 9.34
  sigma2_nu <- 41.15^2
  K <- 1e6  # Assuming unit variance prior for theta0
  
  # Construct P matrix
  P <- matrix(0, 2, 2)
  P[1,1] <- J_tilde_nu_nu + 1/sigma2_nu
  P[1,2] <- P[2,1] <- J_tilde_nu_theta0
  P[2,2] <- J_tilde_theta0_theta0 + 1/K
  
  # Construct h vector
  h <- c(
    J_tilde_nu_nu * nu + J_tilde_nu_theta0 * theta0 + mu_nu/sigma2_nu,
    J_tilde_nu_theta0 * nu + J_tilde_theta0_theta0 * theta0
  )

  if (pt + p > 0) {
    h[1] <- h[1] + khat_Jnu / (1 + g)
    h[2] <- h[2] + khat_Jtheta0 / (1 + g)
  }
  
  # Calculate constant term C
  C <- -mu_nu^2 / (2 * sigma2_nu) - 0.5 * nu^2 * J_nu_nu - nu * theta0 * J_nu_theta0 - 0.5 * theta0^2 * J_theta0_theta0
  
  if (pt + p > 0) {
    # Add the remaining terms involving J_theta_eta
    factor_half <- g / (2 * (1 + g))
    factor_full <- g / (1 + g)
    
    C <- C + factor_half * nu^2 * (t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_nu)
    C <- C + factor_full * nu * theta0 * (t(J_theta_eta_nu) %*% J_theta_eta_inv %*% J_theta_eta_theta0)
    C <- C + factor_half * theta0^2 * (t(J_theta_eta_theta0) %*% J_theta_eta_inv %*% J_theta_eta_theta0)
    C <- C - W_gamma / (2 * (1 + g))
    C <- C - (khat_Jnu * nu + khat_Jtheta0 * theta0) / (1 + g)
  }
  
  # Calculate marginal log-likelihood
  log_det_P <- determinant(P, logarithm = TRUE)$modulus[1]
  P_inv <- solve(P)

  d_gam        <- pt + p
  term_ell     <- as.numeric(log_lik_map)                          # ell(Psi_hat)
  term_dim     <- -0.5 * d_gam * log(1 + g)                        # (A)
  term_logdetP <- -0.5 * log_det_P                                 # (B)
  term_quad    <-  0.5 * as.numeric(t(h) %*% P_inv %*% h)          # (C1)
  term_const   <-  as.numeric(C)                                   # (C2)

  marg_log_lik <- term_ell + term_dim + term_logdetP + term_quad + term_const

  if (isTRUE(return.terms)) {
    n_obs  <- nrow(as.matrix(data$x))
    n_ev   <- sum(as.numeric(data$status) == 1)
    k_gam  <- d_gam + 2
    bic_pen_n  <- -0.5 * k_gam * log(n_obs)
    bic_pen_ev <- -0.5 * k_gam * log(max(n_ev, 1))

    return(list(
      gam          = paste0(par.orig$gam, collapse = ""),
      d_gam        = d_gam,          # number of regression coefficients (pt + p)
      g            = g,              # = n * g_E in this parametrisation
      n_obs        = n_obs,
      n_events     = n_ev,
      ## --- ILA terms (Eq. 11) ---
      ell          = term_ell,       # shared maximised log-likelihood
      ila_dim      = term_dim,       # (A) refined dimension penalty
      ila_logdetP  = term_logdetP,   # (B) common-parameter curvature
      ila_quad     = term_quad,      # (C1) prior-location quadratic
      ila_const    = term_const,     # (C2) constant
      ila_extra    = term_dim + term_logdetP + term_quad + term_const,
      ila_total    = as.numeric(marg_log_lik),
      ## --- BIC penalties on the same model ---
      bic_pen_n    = bic_pen_n,
      bic_pen_ev   = bic_pen_ev,
      bic_total_n  = term_ell + bic_pen_n,
      bic_total_ev = term_ell + bic_pen_ev
    ))
  }

  return(as.numeric(marg_log_lik))
}

targetdist = function(par, data, max.iter = 100, ...){
  if (sum(par$gam)==0 && !(exists("admit.null") && isTRUE(admit.null))){
    return(-Inf)
  }
  log.approx.marg.llh <- marg.loglik(par, data, max.iter = max.iter)
  return(prior.gam(par,log=T) + log.approx.marg.llh)
}