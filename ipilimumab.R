library(mvtnorm)
library(survival)
library(glmnet)
library(HazReg)
set.seed(1)

df0 <- read.csv("./data/CA184043_2A.csv", header = TRUE)
trt0 <- vector()
trt0 <- ifelse(as.factor(df0$Arm) == "ipilimumab", 1 , 0)
df <- data.frame(time = df0$Time, status = df0$Event, trt = trt0)
df <- df[order(df$time),]

status <- as.logical(df$status)
survLTimes <- df$time
trt <- df$trt

#######################################################
#   Original study: single treatment variable only    #
#######################################################
dat <- list(x=cbind(trt), xt=cbind(trt), status=status, y=survLTimes)

# Standard approaches
cox.fit <- coxph(Surv(time, status) ~ trt)
aft.fit <- survreg(Surv(time, status) ~ trt, dist = "lognormal")
summary(cox.fit)
summary(aft.fit)

# Joint selection
source("ADS_sampler.R")
mcmcpar <- list(n.iter = 20000, n.burnin = 10000, n.thin = 2, 
                n.chain = 1, printevery = 2000)
modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1),
                             theta=rnorm(ncol(dat$xt),0,0.5), eta=rnorm(ncol(dat$x),0,0.5),
                             g = 1, a_g = NULL, nu.alpha = 0.01, nu.beta = 0.01, 
                             gam = c(3)),
                data = dat)

chains = BVS.MCMC(modelpar, mcmcpar, verbose = T, null.model = T,
                  gprior.mat = "LCM", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "default" # one of ("default", "custom", "uniform")
)
print(chains$PIP)
print(chains$post.prob[1:10,])  # Model PP estimated from the visited models



#######################################################
#  Additional study: an extra normal noise covariate  #
#######################################################
spur <- scale(rnorm(length(status)))
dat <- list(x=cbind(trt,spur), xt=cbind(trt,spur), status=status, y=survLTimes)

# Standard approaches
cox.fit <- coxph(Surv(time, status) ~ trt+spur)
aft.fit <- survreg(Surv(time, status) ~ trt+spur, dist = "lognormal")
summary(cox.fit)
summary(aft.fit)

# Joint selection
source("ADS_sampler.R")
mcmcpar <- list(n.iter = 20000, n.burnin = 10000, n.thin = 2, 
                n.chain = 1, printevery = 2000)
modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1),
                             theta=rnorm(ncol(dat$xt),0,0.5), eta=rnorm(ncol(dat$x),0,0.5),
                             g = 1, a_g = NULL,nu.alpha = 0.01, nu.beta = 0.01,
                             gam = rep(3,ncol(dat$xt))),
                data = dat)

chains = BVS.MCMC(modelpar, mcmcpar, verbose = T, null.model = T,
                  gprior.mat = "LCM", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "default" # one of ("default", "custom", "uniform")
)
print(chains$PIP)
print(chains$post.prob[1:10,])  # Model PP estimated from the visited models