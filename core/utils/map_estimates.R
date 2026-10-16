if (!exists("BVS_CORE")) BVS_CORE <- "."
source(file.path(BVS_CORE, "utils", "mle_gh.R"))

map.base.par <- function(dat, gam, g.E = 1){
  list(nu = 0, theta0 = 0, theta = rep(0, ncol(dat$xt)), eta = rep(0, ncol(dat$x)),
       n = nrow(dat$x), g = nrow(dat$x) * g.E, g.theta = 1, g.eta = 1, g.aft = 1,
       a_g = NULL, gam = gam, nu.alpha = 0.01, nu.beta = 0.01)
}

map.to.orig <- function(M){
  b <- gam.blocks(M$gam)
  eta <- M$eta
  eta[M$gam == 4] <- exp(-M$nu) * M$theta[M$gam == 4]
  sigma <- exp(-M$nu)
  alpha <- -M$theta * sigma; beta <- -eta
  names(alpha) <- names(beta) <- colnames(M$x.names)
  list(gam = M$gam, mu = M$theta0 * sigma, sigma = sigma, alpha = alpha, beta = beta)
}

# gam: model vector; dat: data list; hyper: integrate a robust prior on g_E over g_E
map.lcm <- function(gam, dat, hyper = FALSE, g.E = 1, K = 1e6, max.iter = 1000){
  par1 <- mle.gh(map.base.par(dat, gam, g.E), dat, max.iter = max.iter)
  b <- gam.blocks(gam); k <- sum(b$theta) + sum(b$eta)
  keep <- c(FALSE, FALSE, b$theta, b$eta)
  n <- nrow(dat$x)
  if (k > 0){
    J <- -hessLLS(par1, dat)[keep, keep, drop = FALSE]
    ldJ <- as.numeric(determinant(J, logarithm = TRUE)$modulus)
  }
  log.prior.kappa <- function(kap){
    if (k == 0) return(0)
    Q <- as.numeric(t(kap) %*% J %*% kap)
    if (!hyper){
      g <- n * g.E
      return(-0.5 * k * log(2 * pi * g) + 0.5 * ldJ - Q / (2 * g))
    }
    a <- (1 + n) / (k + 1) - 1
    lf <- function(t){
      g <- a + exp(t)
      -0.5 * k * log(2 * pi * g) + 0.5 * ldJ - Q / (2 * g) +
        log(0.5) + 0.5 * log((1 + n) / (k + 1)) - 1.5 * log(1 + g) + t
    }
    m <- optimize(lf, c(-30, 30), maximum = TRUE)$objective
    log(integrate(function(t) exp(lf(t) - m), -Inf, Inf, rel.tol = 1e-8)$value) + m
  }
  negpost <- function(v){
    pr <- vec2par(v, par1, ncol(dat$xt), ncol(dat$x))
    -(lglik(pr, dat) + log.prior.kappa(v[-(1:2)]) +
        dgamma(exp(v[1]), 0.01, 0.01, log = TRUE) + v[1] +
        dnorm(v[2], 0, sqrt(K), log = TRUE))
  }
  start <- c(par1$nu, par1$theta0, par1$theta[b$theta], par1$eta[b$eta])
  opt <- safe.nlminb(start, negpost, iter.max = max.iter, gam = gam)
  M <- vec2par(opt$par, par1, ncol(dat$xt), ncol(dat$x))
  M$x.names <- dat$x[0, , drop = FALSE]
  map.to.orig(M)
}

# gam: model vector; dat: data list
map.product <- function(gam, dat, max.iter = 1000){
  source(file.path(BVS_CORE, "utils", "gprior_setup_fixed.R"), local = TRUE)
  dat$xt2 <- t(dat$xt) %*% dat$xt; dat$x2 <- t(dat$x) %*% dat$x
  dat$xt2.inv <- solve(dat$xt2); dat$x2.inv <- solve(dat$x2)
  par1 <- map.base.par(dat, gam)
  lapprox.margin.likelihood(par1, dat, laplace.iter.max = max.iter)
  M <- utils::modifyList(par1, MAP.pars)
  M$x.names <- dat$x[0, , drop = FALSE]
  map.to.orig(M)
}

# t: times; x: covariate vector; mu, sigma: log-normal baseline; alpha, beta: time- and hazard-level coefficients
S.gh.lnorm <- function(t, x, mu, sigma, alpha, beta){
  la <- sum(x * alpha); lb <- sum(x * beta)
  H0 <- -plnorm(t * exp(la), mu, sigma, lower.tail = FALSE, log.p = TRUE)
  exp(-H0 * exp(lb - la))
}
