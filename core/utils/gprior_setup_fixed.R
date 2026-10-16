if (!exists("BVS_CORE")) BVS_CORE <- "."
source(file.path(BVS_CORE, "utils", "mle_gh.R"))
K.theta0 <- 1e6

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

  lp.theta0 = dnorm(par$theta0, 0, sqrt(K.theta0), log=TRUE)
  if (any(par$gam == 4)){
    g.aft = if (is.null(par$g.aft)) par$g.theta else par$g.aft
    n = nrow(data$xt)
    out = dmvnorm(x=theta, mean=rep(0,pxt), sigma=g.aft*n*data$xt2.inv, log=T) +
      dgamma(exp(nu), nu.alpha, nu.beta, log=T) + nu + lp.theta0
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
  out <- cp1 + cp2 + dgamma(exp(nu), nu.alpha, nu.beta, log=T) + nu + lp.theta0
  return(out)
}

gradprior = function(par, data, data.subbed = F){
  if (!data.subbed){
    data = data.gamma(par,data)
  }
  par = par.gamma(par)
  n = length(data$y)
  x2 = data$x2; xt2 = data$xt2
  pxt = ncol(data$xt); px = ncol(data$x)
  theta = par$theta; eta = par$eta

  grad.nu = par$nu.alpha - par$nu.beta*exp(par$nu)
  grad.theta0 = -par$theta0/K.theta0
  if (any(par$gam == 4)){
    n = nrow(data$xt)
    g.aft = if (is.null(par$g.aft)) par$g.theta else par$g.aft
    grad.theta = -1/n*(xt2%*%theta)/g.aft
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
  n = length(data$y)
  x2 = data$x2; xt2 = data$xt2
  pxt = ncol(data$xt); px = ncol(data$x)

  hess.nu = -par$nu.beta*exp(par$nu)
  hess.theta0 = -1/K.theta0
  if (any(par$gam == 4)){
    n = nrow(data$xt)
    g.aft = if (is.null(par$g.aft)) par$g.theta else par$g.aft
    hess.theta = -xt2/(1*n*g.aft)
    hess.eta = numeric(0)
  }else{
    hess.theta = if(is.null(pxt)){numeric(0)}else{-xt2/(1*n*par$g.theta)}
    hess.eta = if(is.null(px)){numeric(0)}else{-x2/(1*n*par$g.eta)}
  }
  pxt = ifelse(is.null(pxt), 0, pxt)
  px = ifelse(is.null(px), 0, px)

  HESS <-  as.matrix(rbind(rbind(c(hess.nu, 0, rep(0,pxt), rep(0,px))),
                           rbind(c(0, hess.theta0, rep(0,pxt), rep(0,px))),
                           if(pxt==0){numeric(0)}else{cbind(cbind(rep(0,pxt)), cbind(rep(0,pxt)), hess.theta, matrix(0,nrow=pxt,ncol=px))},
                           if(px==0){numeric(0)}else{cbind(cbind(rep(0,px)), cbind(rep(0,px)), matrix(0,nrow=px,ncol=pxt), hess.eta)}))
  return(HESS)
}

neg.post = function(par.vec,otherpar,data2,data2.sub,
                    lglik2=lglik,lgprior2=lgprior,...){
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

  return(-lglik2(pars,data2) - lgprior2(pars, data2.sub, data.subbed = T))
}

grad.post = function(par.vec,otherpar,data2,data2.sub,
                     gradlglik=gradLLS,gradlgprior=gradprior,...){
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
  out <- -gradlglik(pars,data2)[out.keep] - gradlgprior(pars, data2.sub, data.subbed = T)
  return(out)
}

hess.post = function(par.vec,otherpar,data2,data2.sub,
                     hesslglik=hessLLS,hesslgprior=hessprior,...){
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
  out <- -hesslglik(pars,data2)[rows.keep,cols.keep] - hesslgprior(pars, data2.sub, data.subbed = T)
  return(out)
}

lapprox.margin.likelihood <- function(par1, data1,
                                      loglik=lglik, logprior=lgprior,
                                      grad.prior=gradprior, hess.prior=hessprior,
                                      grad.lik=gradLLS, hess.lik=hessLLS,
                                      laplace.iter.max=200, reinit=F){
  pxt = ncol(data1$xt); px = ncol(data1$x)
  theta.keep = as.logical(par1$gam == 1 | par1$gam == 3 | par1$gam == 4)
  eta.keep = as.logical(par1$gam == 2 | par1$gam == 3)
  ptheta = sum(theta.keep)
  peta = sum(eta.keep)
  D = 2 + ptheta + peta
  data1.sub = data.gamma(par1,data1)

  result <- safe.nlminb(rep(0,D), neg.post, grad.post, hess.post,
                        otherpar=par1, data2.sub = data1.sub, data2=data1,
                        iter.max = laplace.iter.max, gam = par1$gam)

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

  MAP.pars <<- MAP
  MAP = utils::modifyList(par1, MAP)
  hess.MAP <<- hess.lik(MAP,data1)[c(T,T,theta.keep,eta.keep),c(T,T,theta.keep,eta.keep)] +
    hess.prior(MAP,data1.sub, data.subbed = T)
  if (!is.pd(-hess.MAP)){
    log.opt.problem(par1$gam, "negative Hessian not positive definite at the MAP")
    return(-Inf)
  }
  log.marg.post = loglik(MAP,data1)+logprior(MAP,data1.sub, data.subbed = T) + D/2*log(2*pi) -
    as.vector(0.5*determinant(-hess.MAP,logarithm = T)$modulus)
  return(log.marg.post)
}

targetdist = function(par, data, loglik=lglik,
                      prior=prior.gam, prior.type = "gam",...){
  if (sum(par$gam)==0 && !(exists("admit.null") && isTRUE(admit.null))){
    return(-Inf)
  }
  lp <- prior(par, log = TRUE)
  if (!is.finite(lp)) return(-Inf)
  log.approx.margin.llh = lapprox.margin.likelihood(par, data, laplace.iter.max = 200)
  return(lp + log.approx.margin.llh)
}
