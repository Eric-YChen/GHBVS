# Probability Density Function: location parameter mu and scale parameter sigma
dlap = function(x,mu=0,sigma=1,log=FALSE){
  log.den = - log(2) - log(sigma) - abs(x-mu)/sigma
  if(log) return(log.den)
  else return(exp(log.den))
}

# Cumulative Distribution Function
plap = function(x,mu=0,sigma=1,log.p=FALSE){
  val <- rep(0, length(x))
  val[x<mu]<- -log(2) + (x[x<mu]-mu)/sigma
  val[x>=mu] <- -log(2) + (-x[x>=mu]-mu)/sigma
  if(log.p) val[x<mu] = val[x<mu] else val = exp(val[x<mu])
  if(log.p) val[x>=mu] <- log(1 - exp(val[x>=mu])) else 1 - exp(val[x>=mu])
  return(as.vector(val))
}

# Quantile Function
qlap = function(p,mu=0,sigma=1){
  val <- rep(0,length(p))
val[p<0.5] <-  mu + sigma*log(2*p[p<0.5]) 
val[p>=0.5] <-  mu - sigma*log(2*(1-p[p>=0.5])) 
return(as.vector(val))
}

#--------------------------------------------------------------------------------------------------------------------------
# The Normal Inverse Mills Ratio (IMR) function
#--------------------------------------------------------------------------------------------------------------------------
NIMR <- Vectorize( function(x) exp( dnorm(x,log=T) - pnorm(x,log.p = TRUE, lower.tail = FALSE) ) )


lglik = function(par, data,
                 pdf="dnorm", cdf="pnorm"){
  ### Parameters:
  # nu = -log(sigma)
  # theta0 = mu/sigma
  # theta = -alpha/sigma
  # eta = -beta
  ###
  X <- data$x; Xt <- data$xt; times <- data$y; status <- data$status
  X <- as.matrix(X);  Xt <- as.matrix(Xt); times <- as.vector(times)
  p <- ncol(X);  pt <- ncol(Xt)
  status <- as.logical(status); delta <- as.numeric(status);
  yi = log(times); log.t <- sum(yi[status]); n.obs <- sum(status)
  # Parameters
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  
  # useful quantities
  enu <- exp(nu)
  e_nu <- exp(-nu)
  X.theta <- as.vector(Xt%*%theta)
  if (any(par$gam == 4)){  # AFT model
    v.i <- 1
  }else{
    X.eta <- as.vector(X%*%eta)
    v.i <- exp(X.theta*e_nu - X.eta)
  }
  u.i <- as.vector(enu*yi - X.theta - theta0)
  vi.di <- as.vector(v.i - delta)
  
  if(pdf=="dnorm" & cdf=="pnorm"){
    lgf_t = dnorm(u.i,log = T)  # log(f(tau*log(t_i) - x_til^T*theta - theta0))
    lgF_t = pnorm(-u.i,log.p = T)  # log(F(-...))
  }else{
    lgf_t = eval(parse(text = paste0("log(",pdf,"(u.i))")))  # log(f(tau*log(t_i) - x_til^T*theta - theta0))
    lgF_t = eval(parse(text = paste0("log(",cdf,"(-u.i))")))  # log(F(-...))
  }
  
  out = n.obs*nu - sum(yi[status]) + ifelse(any(par$gam == 4),0,sum(log(v.i)[status])) + 
    sum(lgf_t[status]) + sum(lgF_t*vi.di)
  return(as.numeric(out))
}

#--------------------------------------------------------------------------------------------------------------------------
# Gradient function
#--------------------------------------------------------------------------------------------------------------------------
gradLLS <- function(par, data, dist = "lnorm"){
  # Common variables
  X <- data$x; Xt <- data$xt; times <- data$y; status <- data$status
  X <- as.matrix(X);  Xt <- as.matrix(Xt); times <- as.vector(times)
  p <- ncol(X);  pt <- ncol(Xt)
  status <- as.logical(status); delta <- as.numeric(status);
  yi = log(times); log.t <- sum(yi[status]); n.obs <- sum(status)
  # Parameters
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta

  # Common quantities
  enu <- exp(nu)
  e_nu <- exp(-nu)
  X.theta <- as.vector(Xt%*%theta)
  if (any(par$gam == 4)){  # AFT model
    v.i <- 1
  }else{
    X.eta <- as.vector(X%*%eta)
    v.i <- exp(X.theta*e_nu - X.eta)
  }
  u.i <- as.vector(enu*yi - X.theta - theta0)
  vi.di <- as.vector(v.i - delta)
  
  if(dist == "lnorm"){
    fp.fu.i.obs <- -u.i[status]  # -ui*di
    IMR.ui_vi.di <- as.vector(NIMR(u.i)*vi.di)  # f(ui)/F(-ui)*[vi - di]
    log.Fv <- as.vector(pnorm(u.i, log.p = TRUE, lower.tail = FALSE)*v.i)  # F(-ui)*vi
    # Gradient
    grad.theta0 <- -sum(fp.fu.i.obs ) + sum(IMR.ui_vi.di)
    if (any(par$gam == 4)){  # AFT model
      grad.theta <- - colSums(fp.fu.i.obs*Xt[status,,drop=FALSE]) +
        colSums(IMR.ui_vi.di*Xt)
      grad.eta <- rep(0,p)
      grad.nu <- n.obs + sum(fp.fu.i.obs*yi[status])*enu - sum(IMR.ui_vi.di*yi*enu)
    }else{  # Otherwise
      grad.theta <- colSums(Xt[status,,drop=FALSE])*e_nu - colSums(fp.fu.i.obs*Xt[status,,drop=FALSE]) +
        colSums(IMR.ui_vi.di*Xt) + colSums(log.Fv*Xt)*e_nu
      grad.eta <- -colSums(X[status,,drop=FALSE]) - colSums(log.Fv*X)
      grad.nu <- n.obs - sum(X.theta[status])*e_nu + sum(fp.fu.i.obs*yi[status])*enu -
        sum(IMR.ui_vi.di*yi*enu) - sum(log.Fv*X.theta)*e_nu
    }
  }
  if(dist == "llogis"){  ######## AFT model not implemented ########
    fp.fu.i.obs <- - tanh(0.5*u.i[status])
    IMR.ui_vi.di <- as.vector(plogis(u.i)*vi.di)
    log.Fv <- as.vector(plogis(u.i, log.p = TRUE, lower.tail = FALSE)*v.i)
    # Gradient
    grad.theta0 <- -sum(fp.fu.i.obs ) + sum(IMR.ui_vi.di)
    grad.theta <- colSums(Xt[status,,drop=FALSE])*e_nu - colSums(fp.fu.i.obs*Xt[status,,drop=FALSE]) +
      colSums(IMR.ui_vi.di*Xt) + colSums(log.Fv*Xt)*e_nu
    grad.eta <- -colSums(X[status,,drop=FALSE]) - colSums(log.Fv*X)
    grad.nu <- n.obs - sum(X.theta[status])*e_nu + sum(fp.fu.i.obs*yi[status])*enu -
      sum(IMR.ui_vi.di*yi*enu) - sum(log.Fv*X.theta)*e_nu
  } 
  return(as.vector(c(grad.nu,grad.theta0,grad.theta,grad.eta)))
}

#--------------------------------------------------------------------------------------------------------------------------
# Hessian matrix function
#--------------------------------------------------------------------------------------------------------------------------
hessLLS <- function(par, data, dist="lnorm"){
  # Common variables
  X <- data$x; Xt <- data$xt; times <- data$y; status <- data$status
  X <- as.matrix(X);  Xt <- as.matrix(Xt); times <- as.vector(times)
  p <- dim(X)[2];  pt <- dim(Xt)[2]
  status <- as.logical(status); delta <- as.numeric(status);
  yi = log(times); log.t <- sum(yi[status]); n.obs <- sum(status)
  # Parameters
  nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  
  # Common quantities
  enu <- exp(nu)
  e_nu <- exp(-nu)
  e2nu <- enu^2
  e_2nu <- e_nu^2
  Xt.obs <- as.matrix(Xt[status,,drop=FALSE])
  X.theta <- as.vector(Xt%*%theta)
  X.theta.e_nu <- as.vector(X.theta*e_nu)
  if (any(par$gam == 4)){  # AFT model
    v.i <- 1
  }else{
    X.eta <- as.vector(X%*%eta)
    v.i <- exp(X.theta*e_nu - X.eta)
    X.X <- lapply(1:nrow(X), function(i) tcrossprod(X[i,]))
    Xt.X <- lapply(1:nrow(Xt), function(i) tcrossprod(Xt[i,],X[i,]))
  }
  u.i <- as.vector(enu*yi - X.theta - theta0)
  u.i.obs <- u.i[status]
  u.i2 <- u.i^2
  u.i2.obs <- u.i2[status]
  Xt.Xt <- lapply(1:nrow(Xt), function(i) tcrossprod(Xt[i,]))
  yi2 <- yi^2
  yi.enu <- yi*enu
  yi.enu2 <- yi.enu^2
  
  if(dist == "lnorm"){
    # Common quantities
    fpp.f.obs <- as.vector(-1 + u.i2.obs)
    fp.f.obs2 <- as.vector(u.i2.obs)
    f.F <- as.vector(NIMR(u.i))
    f.F2 <- as.vector(f.F^2)
    fp.F <- as.vector(-u.i*f.F)
    vi.di <- as.vector(v.i - delta)
    log.Fv <- as.vector(pnorm(u.i, log.p = TRUE, lower.tail = FALSE)*v.i)
    # Entries of the Hessian matrix
    if (any(par$gam == 4)){
      hess.theta0 <- sum(fpp.f.obs - fp.f.obs2)- sum( (fp.F + f.F2)*vi.di )
      w.theta.obs <- (fpp.f.obs - fp.f.obs2)
      w.theta <-  -(fp.F + f.F2)*vi.di 
      w.theta <- as.vector(w.theta); w.theta.obs <- as.vector(w.theta.obs)
      hess.theta <- Reduce('+', Map(`*`, Xt.Xt[status], w.theta.obs)) + Reduce('+', Map(`*`, Xt.Xt, w.theta)) 
      hess.eta <- matrix(0, nrow = p, ncol = p)
      hess.nu <-  sum((fpp.f.obs - fp.f.obs2)*yi.enu2[status]) + sum( -u.i.obs*yi.enu[status] ) - 
        sum((fp.F + f.F2)*yi.enu2*vi.di) - sum( f.F*vi.di*yi.enu)
      hess.theta0.theta <- colSums( (fpp.f.obs - fp.f.obs2)*Xt.obs ) - 
        colSums( (fp.F + f.F2)*vi.di*Xt )
      hess.theta0.eta <- rep(0,p)
      hess.theta0.nu <- -sum( (fpp.f.obs - fp.f.obs2)*yi.enu[status] ) + 
        sum( (fp.F + f.F2)*vi.di*yi.enu )
      hess.theta.eta <- matrix(0, nrow = pt, ncol = p)
      hess.eta.nu <- rep(0,p)
      hess.theta.nu <- - colSums((fpp.f.obs - fp.f.obs2)*Xt.obs*yi.enu[status]) + colSums( (fp.F + f.F2)*yi.enu*vi.di*Xt )
      
    }else{
      hess.theta0 <- sum(fpp.f.obs - fp.f.obs2)- sum( (fp.F + f.F2)*vi.di )
      w.theta.obs <- (fpp.f.obs - fp.f.obs2)
      w.theta <-  -(fp.F + f.F2)*vi.di + 2*f.F*v.i*e_nu + log.Fv * e_2nu
      w.theta <- as.vector(w.theta); w.theta.obs <- as.vector(w.theta.obs)
      hess.theta <- Reduce('+', Map(`*`, Xt.Xt[status], w.theta.obs)) + Reduce('+', Map(`*`, Xt.Xt, w.theta))       
      hess.eta <- Reduce('+', Map(`*`, X.X, log.Fv)) 
      # fp.fobs*u.i.obs = -u.i2.obs
      hess.nu <-  sum(X.theta.e_nu[status]) + sum((fpp.f.obs - fp.f.obs2)*yi.enu2[status]) + 
        sum( -u.i.obs*yi.enu[status] ) - sum((fp.F + f.F2)*yi.enu2*vi.di) + sum(log.Fv*X.theta.e_nu*(1+X.theta.e_nu)) + 
        2*sum( f.F*yi*v.i*X.theta) - sum( f.F*vi.di*yi.enu)
      
      hess.theta0.theta <- colSums( (fpp.f.obs - fp.f.obs2)*Xt.obs ) + colSums( f.F*v.i*Xt )*e_nu - 
        colSums( (fp.F + f.F2)*vi.di*Xt )
      hess.theta0.eta <- - colSums( f.F*v.i*X )
      # fp.f = -u.i here
      hess.theta0.nu <- -sum( (fpp.f.obs - fp.f.obs2)*yi.enu[status] ) + 
        sum( (fp.F + f.F2)*vi.di*yi.enu ) - sum(f.F*v.i*X.theta.e_nu)
      
      w.theta.eta <-  as.vector(-(f.F)*v.i - log.Fv*e_nu)
      hess.theta.eta <- Reduce('+', Map(`*`, Xt.X, w.theta.eta))    
      hess.theta.nu <- colSums(-Xt.obs)*e_nu - colSums((fpp.f.obs - fp.f.obs2)*Xt.obs*yi.enu[status]) - 
        colSums(f.F*v.i*Xt*X.theta.e_nu) - colSums(log.Fv*Xt*(1+X.theta.e_nu))*e_nu +
        colSums( (fp.F + f.F2)*yi.enu*vi.di*Xt ) - colSums(f.F*yi*v.i*Xt)
      
      hess.eta.nu <- colSums(f.F*v.i*X*yi.enu) + colSums(log.Fv*X*X.theta.e_nu)
    }
  }
  if(dist == "llogis"){  ######## AFT model not implemented ########
    # Common quantities
    f.F <- as.vector(plogis(u.i))
    f.F2 <- as.vector(f.F^2)
    fp.f <- -as.vector(tanh(0.5*u.i))
    fp.fobs <- as.vector(fp.f[status])
    fp.f.obs2 <- as.vector(fp.fobs^2)
    fp.F <- -as.vector( dlogis(u.i)*( exp(u.i)-1 ) )
    fpp.f <- as.vector(1 - 3/( cosh(u.i) + 1) )
    fpp.f.obs <- as.vector(fpp.f[status])
    vi.di <- as.vector(v.i - delta)
    log.Fv <- as.vector(plogis(u.i, log.p = TRUE, lower.tail = FALSE)*v.i)
    # Entries of the Hessian matrix
    hess.theta0 <- sum(fpp.f.obs - fp.f.obs2)- sum( (fp.F + f.F2)*vi.di )
    w.theta.obs <- (fpp.f.obs - fp.f.obs2)
    w.theta <-  -(fp.F + f.F2)*vi.di + 2*f.F*v.i*e_nu + log.Fv * e_2nu
    w.theta <- as.vector(w.theta); w.theta.obs <- as.vector(w.theta.obs)
    hess.theta <- Reduce('+', Map(`*`, Xt.Xt[status], w.theta.obs)) + Reduce('+', Map(`*`, Xt.Xt, w.theta))       
    hess.eta <- Reduce('+', Map(`*`, X.X, log.Fv)) 
    # fp.fobs*u.i.obs = -u.i2.obs
    hess.nu <-  sum(X.theta.e_nu[status]) + sum((fpp.f.obs - fp.f.obs2)*yi.enu2[status]) + 
      sum( fp.fobs*yi.enu[status] ) - sum((fp.F + f.F2)*yi.enu2*vi.di) + sum(log.Fv*X.theta.e_nu*(1+X.theta.e_nu)) + 
      2*sum( f.F*yi*v.i*X.theta) - sum( f.F*vi.di*yi.enu)
    hess.theta0.theta <- colSums( (fpp.f.obs - fp.f.obs2)*Xt.obs ) + colSums( f.F*v.i*Xt )*e_nu - 
      colSums( (fp.F + f.F2)*vi.di*Xt )
    hess.theta0.eta <- - colSums( f.F*v.i*X )
    # fp.f = -u.i here
    hess.theta0.nu <- -sum( (fpp.f.obs - fp.f.obs2)*yi.enu[status] ) + 
      sum( (fp.F + f.F2)*vi.di*yi.enu ) - sum(f.F*v.i*X.theta.e_nu)
    
    w.theta.eta <-  as.vector(-(f.F)*v.i - log.Fv*e_nu)
    hess.theta.eta <- Reduce('+', Map(`*`, Xt.X, w.theta.eta))    
    hess.theta.nu <- colSums(-Xt.obs)*e_nu - colSums((fpp.f.obs - fp.f.obs2)*Xt.obs*yi.enu[status]) - 
      colSums(f.F*v.i*Xt*X.theta.e_nu) - colSums(log.Fv*Xt*(1+X.theta.e_nu))*e_nu +
      colSums( (fp.F + f.F2)*yi.enu*vi.di*Xt ) - colSums(f.F*yi*v.i*Xt)
    
    hess.eta.nu <- colSums(f.F*v.i*X*yi.enu) + colSums(log.Fv*X*X.theta.e_nu)
  }
  HESS <-  as.matrix(rbind(c(hess.nu, hess.theta0.nu, hess.theta.nu, hess.eta.nu),
                           c(hess.theta0.nu, hess.theta0, hess.theta0.theta, hess.theta0.eta),
                           cbind(cbind(hess.theta.nu), cbind(hess.theta0.theta), hess.theta, hess.theta.eta),
                           cbind(cbind(hess.eta.nu), cbind(hess.theta0.eta), t(hess.theta.eta), hess.eta)))
  
  colnames(HESS) <- NULL;   rownames(HESS) <- NULL
  return(HESS)
}


#--------------------------------------------------------------------------------------------------------------------------
# Norm functions
#--------------------------------------------------------------------------------------------------------------------------
norm_vec <- function(x) sqrt(sum(x^2))
norm2_vec <- function(x) sum(x^2)
#--------------------------------------------------------------------------------------------------------------------------
# MLE using the Gradient Descent Method 
# method: "fixed" for gamma fixed, 
#         "BB" for Barzilai-Borwein adaptive method
#         "MM" for Malitsky-Mishchenko adaptive method  
# gamma0:  Initial value of the step size for "BB" and "MM" methods
# gamma:   Step size for "fixed" method
#--------------------------------------------------------------------------------------------------------------------------
GDMLE <- function(init, dist, data, method = "BB", 
                  in.original.par = F, out.original.par = F,
                  gamma = NULL, gamma0 = NULL, maxit = 100, tol = 1e-3, lambda = 1e-6){
  # Common variables
  X <- data$x; Xt <- data$xt; times <- data$y; status <- data$status
  X <- as.matrix(X);  Xt <- as.matrix(Xt); times <- as.vector(times)
  p <- ncol(X);  pt <- ncol(Xt)
  status <- as.logical(status); delta <- as.numeric(status);
  yi = log(times); log.t <- sum(yi[status]); n.obs <- sum(status)
  par = init
  if (in.original.par){
    sigma <- par[1]
    nu <- -log(sigma); theta0 <- par[2]/sigma; theta <- -par[3:(pt+2)]/sigma; eta <- -par[(pt+3):(p+pt+2)]
  }else{
    nu <- par$nu; theta0 <- par$theta0; theta <- par$theta; eta <- par$eta
  }
  par <- c(nu, theta0, theta, eta)
  
  if(method == "fixed"){
    iter <- 1
    grad.par <- gradLLS(init, dist = dist, times, status, X, Xt)
    while(iter < maxit & norm_vec(grad.par) > tol){
      # Parameters
      #sigma <- exp(par[1]); theta0 <- par[2]; theta <- par[3:(pt+2)]; eta <- par[(pt+3):(p+pt+2)]
      par <- par - gamma*grad.par
      grad.par <- gradLLS(par, dist = dist, times, status, X, Xt)
      iter <- iter + 1
    }
    grad.par1 <- grad.par
  }
  if(method == "BB"){
    iter <- 1
    grad.par1 <- gradLLS(init, dist = dist, times, status, X, Xt)
    while(iter < maxit & norm_vec(grad.par1) > tol){
      # GD cycle
      grad.par0 <- grad.par1
      par0 <- par
      par <- as.vector(par + gamma0*grad.par0)
      grad.par1 <- gradLLS(par, dist = dist, times, status, X, Xt)
      # adaptive value of gamma
      gamma0 <- as.numeric(abs( (par - par0)%*%(grad.par1 - grad.par0))/norm2_vec(grad.par1 - grad.par0))
      iter <- iter + 1
    }
  }
  if(method == "MM"){
    iter <- 1
    grad.par1 <- gradLLS(init, dist = dist, times, status, X, Xt)
    theta = Inf
    while(iter < maxit & norm_vec(grad.par1) > tol){
      # GD cycle
      grad.par0 <- grad.par1
      par0 <- par
      par <- as.vector(par + lambda*grad.par0)
      grad.par1 <- gradLLS(par, dist = dist, times, status, X, Xt)
      # adaptive value of lambda
      lambda0 <- lambda
      lambda <- as.numeric(min( sqrt(1+theta)*lambda, 0.5*norm_vec( par - par0 )/norm_vec(grad.par1 - grad.par0)))
      theta <- lambda/lambda0
      iter <- iter + 1
    }
  }
  # output
  if(out.original.par){
    par <- c(exp(-par[1]),pargdbb$par[2:(pt+2)]*exp(-par[1]),pargdbb$par[(pt+3):(p+pt+2)])
  }
  out <- list(par = par, maxiter = iter, gradient = grad.par1)
  return(out)
}
