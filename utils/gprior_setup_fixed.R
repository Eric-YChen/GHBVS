################################################################################
# prior settings for fixed-g
################################################################################

# Efficient parameter parsing function (no caching - simpler and more robust)
parse_par_vec_fast <- function(par.vec, otherpar, data2) {
  pars <- otherpar
  pxt <- ncol(data2$xt)
  px <- ncol(data2$x)

  # Pre-compute logical vectors once
  theta.keep <- (pars$gam == 1L) | (pars$gam == 3L) | (pars$gam == 4L)
  eta.keep <- (pars$gam == 2L) | (pars$gam == 3L)
  ptheta <- sum(theta.keep)
  peta <- sum(eta.keep)

  pars$nu <- par.vec[1L]
  pars$theta0 <- par.vec[2L]

  if (ptheta == 0L) {
    pars$theta <- rep(0, pxt)
    pars$eta <- replace(as.numeric(eta.keep), eta.keep, par.vec[3L:(2L + peta)])
  } else if (peta == 0L) {
    pars$eta <- rep(0, px)
    pars$theta <- replace(as.numeric(theta.keep), theta.keep, par.vec[3L:(2L + ptheta)])
  } else {
    pars$theta <- replace(as.numeric(theta.keep), theta.keep, par.vec[3L:(2L + ptheta)])
    pars$eta <- replace(as.numeric(eta.keep), eta.keep, par.vec[(3L + ptheta):(2L + ptheta + peta)])
  }

  list(
    pars = pars,
    theta.keep = theta.keep,
    eta.keep = eta.keep,
    out.keep = c(TRUE, TRUE, theta.keep, eta.keep)
  )
}

lgprior = function(par, data, data.subbed = F){
  if (!data.subbed){
    data = data.gamma(par,data)
  }
  par = par.gamma(par)
  theta = par$theta; eta = par$eta; nu = par$nu
  x2 = data$x2; xt2 = data$xt2
  pxt = ncol(data$xt); px = ncol(data$x)
  g.eta = par$g.eta; g.theta = par$g.theta
  nu.alpha = par$nu.alpha; nu.beta = par$nu.beta
  
  if (any(par$gam == 4)){
    g.aft = par$g.aft
    n = nrow(data$xt)
    out = dmvnorm(x=theta, mean=rep(0,pxt), sigma=g.aft*n*data$xt2.inv, log=T) + 
      dgamma(exp(nu), nu.alpha, nu.beta, log=T) + nu
    return(out)
  }
  
  if (is.null(pxt)){cp1 = 0}else{
    n = nrow(data$xt)
    cp1 = dmvnorm(x=theta, mean=rep(0,pxt), sigma=g.theta*n*data$xt2.inv, log=T)
  }
  if (is.null(px)){cp2 = 0}else{
    n = nrow(data$x)
    cp2 = dmvnorm(x=eta, mean=rep(0,px), sigma=g.eta*n*data$x2.inv, log=T)
  }
  out <- cp1 + cp2 + dgamma(exp(nu), nu.alpha, nu.beta, log=T) + nu
  return(out)
}

gradprior = function(par, data, data.subbed = F){
  if (!data.subbed){
    data = data.gamma(par,data)
  }
  par = par.gamma(par)
  n = ifelse(is.null(nrow(data$xt)), nrow(data$x), nrow(data$xt))
  x2 = data$x2; xt2 = data$xt2
  pxt = ncol(data$xt); px = ncol(data$x)
  theta = par$theta; eta = par$eta
  
  grad.nu = par$nu.alpha - par$nu.beta*exp(par$nu)
  grad.theta0 = 0
  if (any(par$gam == 4)){
    n = nrow(data$xt)
    g.aft = par$g.aft
    grad.theta = -1/n*(xt2%*%theta)/par$g.aft
    grad.eta = numeric(0)
  }else{
    grad.theta = if(is.null(pxt)){numeric(0)}else{-1/n*(xt2%*%theta)/par$g.theta}
    grad.eta = if(is.null(px)){numeric(0)}else{-1/n*(x2%*%eta)/par$g.eta}
  }
  
  return(as.vector(c(grad.nu, grad.theta0, grad.theta, grad.eta)))
}

hessprior = function(par, data, data.subbed = F){
  if (!data.subbed){
    data = data.gamma(par,data)
  }
  par = par.gamma(par)
  n = ifelse(is.null(nrow(data$xt)), nrow(data$x), nrow(data$xt))
  x2 = data$x2; xt2 = data$xt2
  pxt = ncol(data$xt); px = ncol(data$x)
  
  hess.nu = -par$nu.beta*exp(par$nu)
  hess.theta0 = 0
  if (any(par$gam == 4)){
    n = nrow(data$xt)
    hess.theta = -xt2/(1*n*par$g.aft)
    hess.eta = numeric(0)
  }else{
    hess.theta = if(is.null(pxt)){numeric(0)}else{-xt2/(1*n*par$g.theta)}
    hess.eta = if(is.null(px)){numeric(0)}else{-x2/(1*n*par$g.eta)}
  }
  pxt = ifelse(is.null(pxt), 0, pxt)
  px = ifelse(is.null(px), 0, px)
  
  HESS <- as.matrix(rbind(rbind(c(hess.nu, 0, rep(0,pxt), rep(0,px))),
                           rbind(c(0, hess.theta0, rep(0,pxt), rep(0,px))),
                           if(pxt==0){numeric(0)}else{cbind(cbind(rep(0,pxt)), cbind(rep(0,pxt)), hess.theta, matrix(0,nrow=pxt,ncol=px))},
                           if(px==0){numeric(0)}else{cbind(cbind(rep(0,px)), cbind(rep(0,px)), matrix(0,nrow=px,ncol=pxt), hess.eta)}))
  return(HESS)
}

neg.post = function(par.vec,otherpar,data2,data2.sub,
                    lglik2=lglik,lgprior2=lgprior,...){
  parsed <- parse_par_vec_fast(par.vec, otherpar, data2)
  pars <- parsed$pars
  return(-lglik2(pars,data2) - lgprior2(pars, data2.sub, data.subbed = T))
}

grad.post = function(par.vec,otherpar,data2,data2.sub,
                     gradlglik=gradLLS,gradlgprior=gradprior,...){
  parsed <- parse_par_vec_fast(par.vec, otherpar, data2)
  pars <- parsed$pars
  out.keep <- parsed$out.keep
  out <- -gradlglik(pars,data2)[out.keep] - gradlgprior(pars, data2.sub, data.subbed = T)
  return(out)
}

hess.post = function(par.vec,otherpar,data2,data2.sub,
                     hesslglik=hessLLS,hesslgprior=hessprior,...){
  parsed <- parse_par_vec_fast(par.vec, otherpar, data2)
  pars <- parsed$pars
  idx.keep <- parsed$out.keep
  out <- -hesslglik(pars,data2)[idx.keep,idx.keep] - hesslgprior(pars, data2.sub, data.subbed = T)
  return(out)
}

# Laplace approximation with nlminb for MAP
lapprox.margin.likelihood <- function(par1, data1,
                                      loglik=lglik, logprior=lgprior, 
                                      grad.prior=gradprior, hess.prior=hessprior,
                                      grad.lik=gradLLS, hess.lik=hessLLS, 
                                      laplace.iter.max=50, reinit=F){
  # Intialize
  pxt = ncol(data1$xt); px = ncol(data1$x)
  theta.keep = as.logical(par1$gam == 1 | par1$gam == 3 | par1$gam == 4)
  eta.keep = as.logical(par1$gam == 2 | par1$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  D = 2 + ptheta + peta
  data1.sub = data.gamma(par1,data1)
  
  # Calculate MAP with nlminb
  result <- nlminb(start=rep(0,D), otherpar=par1, data2.sub = data1.sub,
                  data2=data1, control = list(iter.max = laplace.iter.max),
                  objective=neg.post, gradient=grad.post, hessian=hess.post)
  
  # store number of GD iterations
  if (!exists(quote(itervec))){
    itervec <<- c(result$iterations)
  }else{
    itervec <<- c(itervec,result$iterations)  
  }
  
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
  
  MAP.pars <<- MAP  # store current MAP as initial values for next iteration
  MAP = utils::modifyList(par1, MAP)
  # Calculate log marginal posterior
  hess.MAP <<- hess.lik(MAP,data1)[c(T,T,theta.keep,eta.keep),c(T,T,theta.keep,eta.keep)] + 
    hess.prior(MAP,data1.sub, data.subbed = T)
  log.marg.post = loglik(MAP,data1)+logprior(MAP,data1.sub, data.subbed = T) + sum(par1$gam != 0)/2*log(2*pi) -
    as.vector(0.5*determinant(hess.MAP,logarithm = T)$modulus)
  return(log.marg.post)
}

targetdist = function(par, data, loglik=lglik, 
                      prior=prior.gam, prior.type = "gam",...){
  if (sum(par$gam)==0){
    return(-Inf)
  }
  log.approx.margin.llh = lapprox.margin.likelihood(par, data, laplace.iter.max = 50)
  return(prior(par,log=T) + log.approx.margin.llh)
}
