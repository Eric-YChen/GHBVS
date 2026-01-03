rm(list=ls())
library(HazReg)
library(mvtnorm)
library(tictoc)
library(survival)
library(glmnet)
source("ADS_sampler.R")

set.seed(1234)
n = 250
p = 10
sigma = 0.7**abs(matrix(rep(c(1:p),p),nrow=p,byrow = T) - c(1:p))
hstr = "GH"
baseline = "LN"

dat = simdat(n = n, p = p, sigma = sigma,
             cens = "0.5default", # could be a number in (0,1]
             hstr = hstr, baseline = baseline, n.simset = 1)[[1]][[1]]

mcmcpar = list(n.iter = 20000, n.burnin = 10000, n.thin = 2, 
               n.chain = 1, printevery = 5000)
modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1), theta=rnorm(p,0,0.5), eta=rnorm(p,0,0.5),
                             g.theta = 1, g.eta = 1, g.aft = 1, # for Product prior
                             g = 1, # for LCM prior
                             a_g = NULL, # for hyper-g LCM prior, NULL for robust prior, 2 < a_g <= 4 for hyper-g/n prior
                             gam = init.gam(dat,history.keep = F), # initialise
                             nu.alpha = 0.01, nu.beta = 0.01),
                gam.param = list(layer1.betabin.par = c(1,1),
                                 layer2.mulnom.par = c(0.25,0.25,0.25,0.25), # (AH,PH,GH,AFT)
                                 layer3.beta.par = 1/3),
                data = dat)

print("Starting model:")
print(modelpar$inits$gam)

tic()
chains = BVS.MCMC(modelpar, mcmcpar, verbose = T,
                  gprior.mat = "LCM", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "default" # one of ("default", "custom", "uniform")
                  )
time_cpp = toc()

print(chains$PIP)
print(chains$model.freq[1:10])  # Model PP based on proportion of MCMC visits
print(chains$post.prob[1:10,])  # Model PP estimated from the visited models
