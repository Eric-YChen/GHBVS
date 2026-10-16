# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log: log scale (TRUE or FALSE)
dpgw <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(nu) - log(gamma) - nu*log(sigma) + (nu-1)*log(t) +
    (1/gamma - 1)*log( 1 + (t/sigma)^nu ) +
    ( 1 - ( 1 + (t/sigma)^nu )^(1/gamma) )
  if(log) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log.p: log scale (TRUE or FALSE)
spgw <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- 1 - ( 1 + (t/sigma)^nu )^(1/gamma)
  if(log.p) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log: log scale (TRUE or FALSE)
hpgw <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(nu) - log(gamma) - nu*log(sigma) + (nu-1)*log(t) +
    (1/gamma - 1)*log( 1 + (t/sigma)^nu )
  if(log) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument
chpgw <- function(t, sigma, nu, gamma){
  val <- -1 + ( 1 + (t/sigma)^nu )^(1/gamma)
  return(val)
}

# sigma: scale parameter; nu, gamma: shape parameters; n: number of observations
rpgw <- function(n, sigma, nu, gamma){
  p <- runif(n)
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(as.vector(out))
}

# sigma: scale parameter; nu, gamma: shape parameters; p: probability in (0,1)
qpgw <- function(p, sigma, nu, gamma){
  out <- sigma*(  ( 1 - log(1-p) )^gamma - 1 )^(1/nu)
  return(out)
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log: log scale (TRUE or FALSE)
hew <- function(t, sigma, nu, gamma, log=FALSE){
  log.pdf <-  log(gamma) + (gamma-1)*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) +
    dweibull(t,scale=sigma,shape=nu,log=TRUE)
  cdf <- exp(gamma*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) )
  log.h <- log.pdf - log(1-cdf)
  ifelse(log, return(log.h), return(exp(log.h)))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument
chew <- function(t, sigma, nu, gamma){
  cdf <- exp(gamma*pweibull(t,scale=sigma,shape=nu,log.p=TRUE) )
  return(-log(1-cdf))
}

# p: probability in (0,1); sigma: scale parameter; nu, gamma: shape parameters
qew <- function(p, sigma, nu, gamma){
  quant <-  qweibull(p^(1/gamma),scale=sigma,shape=nu)
  return(quant)
}

# sigma: scale parameter; nu: shape parameter; t: positive argument; log: log scale (TRUE or FALSE)
hweibull <- function(t, sigma, nu, log=FALSE){
  log.pdf <-  dweibull(t,scale=sigma,shape=nu, log = TRUE)
  log.s <- pweibull(t,scale=sigma,shape=nu, lower.tail = FALSE, log.p = TRUE)
  log.h <- log.pdf - log.s
  ifelse(log, return(log.h), return(exp(log.h)))
}

# sigma: scale parameter; nu: shape parameter; t: positive argument
chweibull <- function(t, sigma, nu){
  H0 <- -pweibull(t,scale=sigma,shape=nu, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

# mu: mean parameter in the log scale; sigma: scale parameter in the log scale; t: positive argument
# log: log scale (TRUE or FALSE)
hlnorm <- function(t, mu, sigma, log = FALSE){
  lpdf0 <-  dlnorm(t, mu, sigma, log = TRUE)
  ls0 <- plnorm(t, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

# mu: mean parameter in the log scale; sigma: scale parameter in the log scale; t: positive argument
chlnorm <- function(t, mu, sigma){
  H0 <- -plnorm(t, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

# mu: mean parameter in the log scale; sigma: scale parameter in the log scale; t: positive argument
# log: log scale (TRUE or FALSE)
hllogis <- function(t, mu, sigma, log = FALSE){
  lpdf0 <-  dlogis(log(t),mu,sigma, log = TRUE) - log(t)
  ls0 <- plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

# mu: mean parameter in the log scale; sigma: scale parameter in the log scale; t: positive argument
chllogis <- function(t, mu, sigma){
  H0 <- -plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

# mu: mean parameter in the log scale; sigma: scale parameter in the log scale
# p: probability in (0,1)
qllogis <- function(p, mu, sigma){
  qq <- exp(qlogis(p, mu, sigma))
  return(qq)
}

# scale: scale parameter; shape: shape parameter; t: positive argument; log: log scale (TRUE or FALSE)
hgamma <- function(t, shape, scale, log = FALSE){
  lpdf0 <-  dgamma(t, shape = shape, scale = scale, log = TRUE)
  ls0 <- pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

# scale: scale parameter; shape: shape parameter; t: positive argument; log: log scale (TRUE or FALSE)
chgamma <- function(t, shape, scale){
  H0 <- -pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log: log scale (TRUE or FALSE)
dggamma <- function(t, sigma, nu, gamma, log = FALSE){
  val <- log(gamma) - nu*log(sigma) - lgamma(nu/gamma) + (nu - 1)*log(t) -
    (t/sigma)^gamma
  if(log) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log.p: log scale (TRUE or FALSE)
pggamma <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE)
  if(log.p) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log.p: log scale (TRUE or FALSE)
sggamma <- function(t, sigma, nu, gamma, log.p = FALSE){
  val <- pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE, lower.tail =  FALSE)
  if(log.p) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument; log: log scale (TRUE or FALSE)
hggamma <- function(t, sigma, nu, gamma, log = FALSE){
  log.pdf <- dggamma(t, sigma, nu, gamma, log = TRUE)
  log.s <- sggamma(t, sigma, nu, gamma, log.p = TRUE)
  val <- log.pdf - log.s
  if(log) return(val) else return(exp(val))
}

# sigma: scale parameter; nu, gamma: shape parameters; t: positive argument
chggamma <- function(t, sigma, nu, gamma){
  val <- -pgamma( t^gamma, shape = nu/gamma, scale = sigma^gamma, log.p = TRUE, lower.tail =  FALSE)
  return(val)
}

# sigma: scale parameter; nu, gamma: shape parameters; p: probability in (0,1)
qggamma <- function(p, sigma, nu, gamma){
  out <- qgamma(p, shape = nu/gamma, scale = sigma^gamma)^(1/gamma)
  return(out)
}

# sigma: scale parameter; nu, gamma: shape parameters; n: number of observations
rggamma <- function(n, sigma, nu, gamma){
  p <- runif(n)
  out <- qgamma(p, shape = nu/gamma, scale = sigma^gamma)^(1/gamma)
  return(as.vector(out))
}

GHMLE <-
  function(init,
           times,
           status,
           hstr = NULL,
           dist = NULL,
           des = NULL,
           des_t = NULL,
           method = "Nelder-Mead",
           maxit = 100) {
    times <- as.vector(times)
    status <- as.vector(as.logical(status))
    times.obs <- times[status]
    if (!is.null(des)){
      des <- as.matrix(des)
      des.obs <- des[status, ]
    }
    if (!is.null(des_t)){
      des_t <- as.matrix(des_t)
      des_t.obs <- des_t[status, ]
    }

    if (hstr == "baseline") {
      if (dist == "PGW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "EW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chew(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "GenGamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) + sum(chggamma(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "LogNormal") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "Gamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chgamma(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "Weibull") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) + sum(chweibull(times, ae0, be0))
          return(val)
        }
      }
    }

    if (hstr == "PH") {
      if (dist == "PGW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chpgw(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chew(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chggamma(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE)  + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chweibull(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }
    }

    if (hstr == "AFT") {
      if (dist == "PGW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hpgw(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hew(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hggamma(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hlnorm(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chlnorm(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hllogis(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hgamma(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]
          lhaz0 <- hweibull(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }
    }

    if (hstr == "AH") {
      if (dist == "PGW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]
          lhaz0 <- hweibull(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }
    }

    if (hstr == "GH") {
      if (dist == "PGW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "EW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "GenGamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "LogNormal") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "LogLogistic") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "Gamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) + x.beta.obs

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }
    }

    if (method != "nlminb") {
      OPT <- optim(init,
                   log.lik,
                   control = list(maxit = maxit),
                   method = method)
    }
    if (method == "nlminb") {
      OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
    }

    OUT <- list(log_lik = log.lik, OPT = OPT)
    return(OUT)
  }

# FUN: minus log-likelihood function to be used to calculate the confidence intervals
# MLE: maximum likelihood estimator of the parameters of interest; level: confidence level
# index: position of the positive parameters under the original parameterisation
Conf_Int <- function(FUN,MLE,level=0.95,index=NULL){
  sd.int <- abs(qnorm(0.5*(1-level)))
  tempf <- function(par){
    par[index] = exp(par[index])
    return(FUN( par ))
  }
  r.MLE <- MLE
  r.MLE[index] <- log(MLE[index])
  HESS <- hessian(tempf,x=r.MLE)
  Fisher.Info <- solve(HESS)
  Sigma <- sqrt(diag(Fisher.Info))
  U<- r.MLE + sd.int*Sigma
  L<- r.MLE - sd.int*Sigma
  C.I <- cbind(L,U,r.MLE, Sigma)
  names.row <- paste0("par", seq_along(1:length(MLE)))
  names.row[index] <- paste0("log.par", seq_along(index))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Lower","Upper","Transf MLE", "Std. Error")
  return(C.I)
}

# seed: seed for simulation; n: sample size (number of individuals); theta: parameters of the baseline hazard
# beta_h: regression parameters multiplying the hazard for GH model
# beta_t: regression parameters multiplying the time scale for GH model
# beta: regression parameters for AFT, PH, and AH models; des_h: Design matrix for GH model (hazard scale)
# des_t: Design matrix for GH model (time scale); des: Design matrix for AFT, PH, and AH models
# hstr: hazard structure (AH, AFT, PH, GH); baseline: baseline hazard distribution
simGH <- function(seed, n, des = NULL, des_h = NULL, des_t = NULL, theta,
                  beta_h = NULL, beta_t = NULL, beta = NULL, hstr, baseline){
  if(!is.null(des))   des <- as.matrix(des)
  if(!is.null(des_h)) des_h <- as.matrix(des_h)
  if(!is.null(des_t)) des_t <- as.matrix(des_t)

  if(baseline == "LN")      quantf <- function(p) qlnorm(p, theta[1], theta[2])
  if(baseline == "LL")      quantf <- function(p) qllogis(p, theta[1], theta[2])
  if(baseline == "G")       quantf <- function(p) qgamma(p, theta[1], theta[2])
  if(baseline == "W")       quantf <- function(p) qweibull(p, theta[1], theta[2])
  if(baseline == "PGW")     quantf <- function(p) qpgw(p, theta[1], theta[2], theta[3])
  if(baseline == "EW")      quantf <- function(p) qew(p, theta[1], theta[2], theta[3])
  if(baseline == "GG")      quantf <- function(p) qggamma(p, theta[1], theta[2], theta[3])

  set.seed(seed)
  u = runif(n)

  if( hstr == "GH" ){
    exp.xbeta_t  <- exp(des_t%*%beta_t)
    exp.dif <- exp(des_t%*%beta_t - des_h%*%beta_h )

    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  if( hstr == "PH" ){
    exp.xbeta_t  <- 1
    exp.dif <- exp( - des%*%beta )

    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  if( hstr == "AFT" ){
    exp.xbeta_t  <- exp(des%*%beta)
    exp.dif <- 1

    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  if( hstr == "AH" ){
    exp.xbeta_t  <- exp(des%*%beta)
    exp.dif <- exp(des%*%beta)

    p0 <- as.vector(1 - exp(log(1-u)*exp.dif))
    times <- as.vector(quantf(p0)/exp.xbeta_t)
  }

  return(as.vector(times))
}

GEHMLE <-
  function(init,
           times,
           status,
           hp,
           hstr = NULL,
           dist = NULL,
           des = NULL,
           des_t = NULL,
           method = "Nelder-Mead",
           maxit = 100) {
    times <- as.vector(times)
    status <- as.vector(as.logical(status))
    times.obs <- times[status]
    if (!is.null(des)){
      des <- as.matrix(des)
      des.obs <- des[status, ]
    }
    if (!is.null(des_t)){
      des_t <- as.matrix(des_t)
      des_t.obs <- des_t[status, ]
    }
    hp <- as.vector(hp)
    hp.obs <- as.vector(hp[status])
    lhp.obs <- log(hp.obs)

    if (hstr == "baseline") {
      if (dist == "PGW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT = cbind(lhp.obs , hpgw(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "EW") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT <- cbind(lhp.obs , hew(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chew(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "GenGamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])

          MAT <- cbind(lhp.obs , hggamma(times.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chggamma(times, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "LogNormal") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hlnorm(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hllogis(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "Gamma") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hgamma(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chgamma(times, ae0, be0))
          return(val)
        }
      }

      if (dist == "Weibull") {
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])

          MAT <- cbind(lhp.obs , hweibull(times.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chweibull(times, ae0, be0))
          return(val)
        }
      }
    }

    if (hstr == "PH") {
      if (dist == "PGW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs , hpgw(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chpgw(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hew(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hggamma(times.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chggamma(times, ae0, be0, ce0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hlnorm(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chlnorm(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hllogis(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) + sum(chllogis(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hgamma(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))

          MAT <- cbind(lhp.obs , hweibull(times.obs, ae0, be0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times, ae0, be0) * exp.x.beta)
          return(val)
        }
      }
    }

    if (hstr == "AFT") {
      if (dist == "PGW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          beta <- par[4:(3 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.beta.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.beta, ae0, be0, ce0))
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chlnorm(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          beta <- par[3:(2 + p)]

          x.beta <- des %*% beta
          x.beta.obs <- x.beta[status]
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.obs <- exp.x.beta[status]

          MAT <- cbind(lhp.obs ,
                         hweibull(times.obs * exp.x.beta.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.beta, ae0, be0))
          return(val)
        }
      }
    }

    if (hstr == "AH") {
      if (dist == "PGW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "EW") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "GenGamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "LogNormal") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "LogLogistic") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "Gamma") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }

      if (dist == "Weibull") {
        p <- ncol(des_t)
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p)]

          exp.x.alpha <- as.vector(exp(des_t %*% alpha))
          exp.x.alpha.obs <- exp.x.alpha[status]

          MAT <- cbind(lhp.obs ,
                         hweibull(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE))
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chweibull(times * exp.x.alpha, ae0, be0) / exp.x.alpha)
          return(val)
        }
      }
    }

    if (hstr == "GH") {
      if (dist == "PGW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hpgw(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chpgw(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "EW") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hew(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) + x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chew(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "GenGamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          ce0 <- exp(par[3])
          alpha <- par[4:(3 + p0)]
          beta <- par[(4 + p0):(3 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hggamma(times.obs * exp.x.alpha.obs, ae0, be0, ce0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <-  -sum(lhaz0) +
            sum(chggamma(times * exp.x.alpha, ae0, be0, ce0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "LogNormal") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hlnorm(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chlnorm(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "LogLogistic") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- par[1]
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hllogis(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chllogis(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }

      if (dist == "Gamma") {
        p0 <- dim(des_t)[2]
        p1 <- dim(des)[2]
        log.lik <- function(par) {
          ae0 <- exp(par[1])
          be0 <- exp(par[2])
          alpha <- par[3:(2 + p0)]
          beta <- par[(3 + p0):(2 + p0 + p1)]

          x.alpha <- des_t %*% alpha
          x.beta <- des %*% beta
          x.alpha.obs <- x.alpha[status]
          x.beta.obs <- x.beta[status]
          exp.x.alpha <- as.vector(exp(x.alpha))
          exp.x.beta <- as.vector(exp(x.beta))
          exp.x.beta.dif <-  as.vector(exp(x.beta - x.alpha))
          exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
          exp.x.beta.obs <- as.vector(exp.x.beta[status])

          MAT <- cbind(lhp.obs ,
                         hgamma(times.obs * exp.x.alpha.obs, ae0, be0, log = TRUE) +
                         x.beta.obs)
          lhaz0 <- apply(MAT,1,logSumExp)

          val <- -sum(lhaz0) +
            sum(chgamma(times * exp.x.alpha, ae0, be0) * exp.x.beta.dif)
          return(sum(val))
        }
      }
    }

    if (method != "nlminb") {
      OPT <- optim(init,
                   log.lik,
                   control = list(maxit = maxit),
                   method = method)
    }
    if (method == "nlminb") {
      OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
    }

    OUT <- list(log_lik = log.lik, OPT = OPT)
    return(OUT)
  }
