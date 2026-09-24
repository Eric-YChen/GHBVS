library(mvtnorm)
library(survival)
library(glmnet)
library(HazReg)
library(parallel)

set.seed(1234)
nki70 <- read.csv("./data/nki70.csv", header = TRUE, stringsAsFactors = T)
X <- nki70[,-c(1:2)]
X$Diam <- as.numeric(X$Diam) - 1
X$N <- as.numeric(X$N)
X$ER <- as.numeric(X$ER) - 1
X$Grade <- as.numeric(X$Grade)
X <- apply(X,2,scale)
dat <- list(y = nki70$time, status = nki70$event, x=X, xt=X)

########################################################
##### Proposed Method #####
########################################################
source("ADS_sampler.R")
mcmcpar <- list(n.iter = 110000, n.burnin = 10000, n.thin = 4, 
                n.chain = 1, printevery = 2000)
modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1),
                             theta=rnorm(ncol(dat$xt),0,0.5), eta=rnorm(ncol(dat$x),0,0.5),
                             g.theta = 1, g.eta = 1, g.aft = 1, # for Product prior
                             g = 1, # for LCM prior
                             a_g = NULL, # for hyper-g LCM prior, NULL for robust prior, 2 < a_g <= 4 for hyper-g/n prior
                             gam = init.gam(dat,history.keep = F), # initialise
                             nu.alpha = 0.01, nu.beta = 0.01),
                gam.param = list(layer1.betabin.par = c(1,1),
                                 layer2.mulnom.par = c(0.25,0.25,0.25,0.25), # (AH,PH,GH,AFT)
                                 layer3.beta.par = 1/3),
                data = dat)

chains = BVS.MCMC(modelpar, mcmcpar, verbose = T,
                  gprior.mat = "LCM", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "custom" # one of ("default", "custom", "uniform")
)

print(chains$PIP)
print(chains$model.freq[1:10])  # Model PP based on proportion of MCMC visits
print(chains$post.prob[1:10,])  # Model PP estimated from the visited models



################ Comparison 1 ##########################
##### mombf: BVS for AFT model with non-local prior #####
########################################################
library(mombf)
ys <- Surv(log(nki70$time), event=nki70$event)
df <- data.frame(ys = ys, as.data.frame(X), time = nki70$time, status = nki70$event)

f1y <- formula(paste('ys ~ ',paste(colnames(X), collapse='+'),sep=''))
priorCoef = zellnerprior(tau=nrow(X))
priorDelta = modelbbprior(1,1)
aft1 <- modelSelection(f1y, data = df,
                      priorCoef = priorCoef, priorDelta = priorDelta,
                      enumerate = FALSE, method = 'Laplace', niter = 110000, 
                      thinning = 4, burnin = 10000,
                      center = FALSE, scale = FALSE)
print(postProb(aft1,10))



################ Comparison 2 ##########################
##### BVSNLP: Cox PiMOM prior #####
########################################################
library(BVSNLP)
y2 <- Surv(time = df$time, df$status)

nlpcox <- bvs(X = data.frame(X), resp = as.matrix(cbind(df$time, df$status)), 
              family = "survival", nlptype = "pMOM", mod_prior = "unif",
              niter = 1000, ncpu = 10)
bvs.post <- function(bvsfit,len.out=10){data.frame(Freq=exp(bvsfit$max_prob_vec-mclust::logsumexp(bvsfit$max_prob_vec))[1:len.out],Vars=unlist(lapply(bvsfit$max_models[1:len.out],function(x){paste0(colnames(X)[t(x)],collapse = ", ")})))}
print(bvs.post(nlpcox,10))

