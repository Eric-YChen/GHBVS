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
    etas[etas==1] = par.vec[(3+ptheta):(2+ptheta+peta)]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    etas = eta.keep
    etas[etas==1] = par.vec[(3+ptheta):(2+ptheta+peta)]
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
    etas[etas==1] = par.vec[3:(2+peta)]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    etas = eta.keep
    etas[etas==1] = par.vec[(3+ptheta):(2+ptheta+peta)]
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
    etas[etas==1] = par.vec[(3+ptheta):(2+ptheta+peta)]
    pars$eta = etas
  }else if (peta == 0){
    pars$eta = rep(0,px)
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    pars$theta = thetas
  }else{
    thetas = theta.keep
    thetas[thetas==1] = par.vec[3:(2+ptheta)]
    etas = eta.keep
    etas[etas==1] = par.vec[(3+ptheta):(2+ptheta+peta)]
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
  
  MAP = result$par
  MAP0 = list(nu=MAP[1], theta0 = MAP[2], theta = theta.keep, eta = eta.keep)
  if (ptheta == 0){
    MAP0$theta = rep(0,pxt)
    MAP0$eta[eta.keep==1] = MAP[3:(2+peta)]
  }else if (peta == 0){
    MAP0$eta = rep(0,px)
    MAP0$theta[theta.keep==1] = MAP[3:(2+ptheta)]
  }else{
    MAP0$theta[theta.keep==1] = MAP[3:(2+ptheta)]
    MAP0$eta[eta.keep==1] = MAP[(3+ptheta):(2+ptheta+peta)]
  }
  MAP = MAP0
  MAP = utils::modifyList(par1, MAP)
  return(MAP)
}

marg.loglik <- function(par, data, max.iter = 100) {
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
  
  # Calculate log-likelihood at MAP (independent of g)
  log_lik_map <- lglik(par.orig, data)
  
  # Calculate the Hessian at MAP (independent of g)
  hess_full <- hessLLS(par.orig, data)
  
  # Extract relevant blocks from the Hessian
  # Order is: nu, theta0, theta, eta
  idx_nu <- 1
  idx_theta0 <- 2
  idx_theta <- which(c(FALSE, FALSE, theta.keep))
  idx_eta <- which(c(FALSE, FALSE, rep(FALSE, length(theta.keep)), eta.keep))
  idx_beta <- c(idx_theta, idx_eta)  # Combined theta and eta indices
  
  # Extract Hessian blocks (reused inside integrand)
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
  } else {
    # Placeholders to avoid R CMD check notes
    quad_nu <- 0
    quad_cross <- 0
    quad_theta0 <- 0
  }
  
  # Prior parameters (for nu and theta0)
  mu_nu <- 9.34 # 9.34
  sigma2_nu <- 41.15^2 # 3.72
  K <- 1e6
  
  # Constants independent of g
  base_C <- -mu_nu^2 / (2 * sigma2_nu) - 0.5 * nu^2 * J_nu_nu - nu * theta0 * J_nu_theta0 - 0.5 * theta0^2 * J_theta0_theta0
  
  # Set up prior for g
  k <- pt + p
  n <- par$n
  
  # Check which prior to use based on par$a_g
  if (is.null(par$a_g) || par$a_g <= 2 || par$a_g > 4) {
    # Original robust prior: pi(g) = 0.5 * sqrt((1+n)/(k+1)) * (1+g)^(-3/2), support g > a
    a <- (1 + n) / (k + 1) - 1
    use_hyper_g <- FALSE
  } else {
    # New hyper-g prior: pi(g) = (a_g-2)/(2*n) * (1/(1+g/n))^(a_g/2), support g > 0
    a_g <- par$a_g
    a <- 0  # support is g > 0
    use_hyper_g <- TRUE
  }
  #message(sprintf("Using hyper-g: %s (a_g=%.3f)", use_hyper_g, par$a_g))
  
  # Log integrand as a function of t where g = a + exp(t)
  log_integrand_t <- function(t) {
    # Vectorized in t
    g_val <- a + exp(t)
    one_plus_g <- 1 + g_val
    factor_full <- g_val / one_plus_g
    factor_half <- 0.5 * factor_full
    
    if (has_beta) {
      J_tilde_nu_nu <- J_nu_nu - factor_full * quad_nu
      J_tilde_nu_theta0 <- J_nu_theta0 - factor_full * quad_cross
      J_tilde_theta0_theta0 <- J_theta0_theta0 - factor_full * quad_theta0
    } else {
      # recycle scalars to vector length
      len <- length(t)
      J_tilde_nu_nu <- rep(J_nu_nu, len)
      J_tilde_nu_theta0 <- rep(J_nu_theta0, len)
      J_tilde_theta0_theta0 <- rep(J_theta0_theta0, len)
    }
    
    # P matrix and its log-det (vectorized)
    P11 <- J_tilde_nu_nu + 1 / sigma2_nu
    P12 <- J_tilde_nu_theta0
    P22 <- J_tilde_theta0_theta0 + 1 / K
    det_P <- P11 * P22 - P12 * P12
    valid <- is.finite(det_P) & det_P > 0 & is.finite(one_plus_g) & (one_plus_g > 0)
    res <- rep(-Inf, length(t))
    if (!any(valid)) return(res)
    
    log_det_P <- log(det_P[valid])
    
    # Solve P^{-1} h via 2x2 explicit formulas (vectorized on valid)
    h1 <- J_tilde_nu_nu[valid] * nu + J_tilde_nu_theta0[valid] * theta0 + mu_nu / sigma2_nu
    h2 <- J_tilde_nu_theta0[valid] * nu + J_tilde_theta0_theta0[valid] * theta0
    invP_h1 <- ( P22[valid] * h1 - P12[valid] * h2) / det_P[valid]
    invP_h2 <- (-P12[valid] * h1 + P11[valid] * h2) / det_P[valid]
    quad_h <- 0.5 * (h1 * invP_h1 + h2 * invP_h2)
    
    # Constant term depending on g via beta-block adjustments (vectorized)
    Cg <- rep(base_C, sum(valid))
    if (has_beta) {
      Cg <- Cg + factor_half[valid] * (nu^2 * quad_nu + theta0^2 * quad_theta0) + factor_full[valid] * (nu * theta0 * quad_cross)
    }
    
    # Prior on g (log) and Jacobian from transformation (vectorized)
    if (use_hyper_g) {
      # Hyper-g prior: pi(g) = (a_g-2)/(2*n) * (1/(1+g/n))^(a_g/2)
      log_prior_g <- log(a_g - 2) - log(2 * n) - (a_g / 2) * log(1 + g_val[valid] / n)
    } else {
      # Original robust prior: pi(g) = 0.5 * sqrt((1+n)/(k+1)) * (1+g)^(-3/2)
      log_prior_g <- log(0.5) + 0.5 * log((1 + n) / (k + 1)) - 1.5 * log(one_plus_g[valid])
    }
    log_jac <- t[valid]
    
    # Combine (vectorized)
    out_valid <- log_lik_map - 0.5 * k * log(one_plus_g[valid]) - 0.5 * log_det_P + quad_h + Cg + log_prior_g + log_jac
    res[valid] <- as.numeric(out_valid)
    res
  }
  
  # Find a stable shift for numerical integration
  grid_t <- seq(-6, 6, length.out = 49)
  vals <- vapply(grid_t, log_integrand_t, numeric(1))
  m_shift <- max(vals[is.finite(vals)])
  if (!is.finite(m_shift)) m_shift <- 0
  
  # Scaled integrand
  f_scaled <- function(t) {
    vals <- log_integrand_t(t) - m_shift
    vals[!is.finite(vals)] <- -Inf
    exp(vals)
  }
  
  # Perform 1D numerical integration over t in (-Inf, Inf)
  int_res <- try(integrate(f_scaled, lower = -Inf, upper = Inf, rel.tol = 1e-6, abs.tol = 0), silent = TRUE)
  if (inherits(int_res, "try-error")) {
    # Fallback: integrate over a wide finite range
    int_res <- integrate(f_scaled, lower = -20, upper = 20, rel.tol = 1e-6, abs.tol = 0)
  }
  log_int <- log(int_res$value) + m_shift
  
  return(as.numeric(log_int))
}

targetdist = function(par, data, ...){
  if (sum(par$gam)==0){
    return(-Inf)
  }
  log.approx.marg.llh <- marg.loglik(par, data)
  return(prior.gam(par,log=T) + log.approx.marg.llh)
}