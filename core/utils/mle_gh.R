gam.blocks <- function(gam){
  list(theta = as.logical(gam == 1 | gam == 3 | gam == 4),
       eta   = as.logical(gam == 2 | gam == 3))
}

# par.vec: free parameters (nu, theta0, theta_gamma, eta_gamma); otherpar: full parameter list
vec2par <- function(par.vec, otherpar, pxt, px){
  b <- gam.blocks(otherpar$gam); pt <- sum(b$theta); pe <- sum(b$eta)
  theta <- rep(0, pxt); eta <- rep(0, px)
  theta[b$theta] <- par.vec[seq_len(pt) + 2]
  eta[b$eta]     <- par.vec[seq_len(pe) + 2 + pt]
  otherpar$nu <- par.vec[1]; otherpar$theta0 <- par.vec[2]
  otherpar$theta <- theta; otherpar$eta <- eta
  otherpar
}

opt.log <- character(0)
log.opt.problem <- function(gam, msg){
  opt.log[paste0(gam, collapse = "")] <<- msg
  invisible(NULL)
}

# start, objective, gradient, hessian, ...: passed to nlminb; iter.max: iteration budget; gam: model logged on failure
safe.nlminb <- function(start, objective, gradient = NULL, hessian = NULL, ...,
                        iter.max = 200, gam = NULL){
  ctrl <- function(k) list(iter.max = k, eval.max = 2 * k)
  res <- nlminb(start, objective, gradient, hessian, ..., control = ctrl(iter.max))
  if (res$convergence != 0 || !is.finite(res$objective)){
    res2 <- nlminb(res$par, objective, gradient, hessian, ..., control = ctrl(10 * iter.max))
    if (is.finite(res2$objective) &&
        (!is.finite(res$objective) || res2$objective <= res$objective)) res <- res2
    if ((res$convergence != 0 || !is.finite(res$objective)) && !is.null(gam))
      log.opt.problem(gam, paste("nlminb did not converge:", res$message))
  }
  res
}

is.pd <- function(M) all(is.finite(M)) && !inherits(try(chol(M), silent = TRUE), "try-error")

mle.negloglik <- function(par.vec, otherpar, data2){
  -lglik(vec2par(par.vec, otherpar, ncol(data2$xt), ncol(data2$x)), data2)
}
mle.grad <- function(par.vec, otherpar, data2){
  pars <- vec2par(par.vec, otherpar, ncol(data2$xt), ncol(data2$x))
  b <- gam.blocks(pars$gam)
  -gradLLS(pars, data2)[c(TRUE, TRUE, b$theta, b$eta)]
}
mle.hess <- function(par.vec, otherpar, data2){
  pars <- vec2par(par.vec, otherpar, ncol(data2$xt), ncol(data2$x))
  b <- gam.blocks(pars$gam); keep <- c(TRUE, TRUE, b$theta, b$eta)
  -hessLLS(pars, data2)[keep, keep, drop = FALSE]
}

# par1: parameter list with model gam; data1: data list
mle.gh <- function(par1, data1, max.iter = 200){
  b <- gam.blocks(par1$gam)
  D <- 2 + sum(b$theta) + sum(b$eta)
  res <- safe.nlminb(rep(0, D), mle.negloglik, mle.grad, mle.hess,
                     otherpar = par1, data2 = data1, iter.max = max.iter, gam = par1$gam)
  vec2par(res$par, par1, ncol(data1$xt), ncol(data1$x))
}
