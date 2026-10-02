library(mvtnorm)
library(survival)
library(glmnet)
library(HazReg)

set.seed(1234)
flchain <- read.csv("./data/flchain.csv", header=TRUE, stringsAsFactors = T)
data.c <- flchain[,c(1,2,4,5,7,8,9,10)]
data.c <- data.c[complete.cases(data.c),]
ind0 <- which(data.c$futime==0)
data.c <- data.c[-ind0,]

X <- as.matrix(cbind(scale(data.c$age),scale(data.c$kappa), scale(data.c$lambda), as.numeric(data.c$sex)-1, scale(data.c$creatinine), scale(data.c$mgus) ))
colnames(X) <- colnames(data.c)[c(1,3,4,2,5,6)]
Xt <- X
status <- as.vector(data.c$death)
times <- as.vector(data.c$futime)/365.24 # in years
dat <- list(x=X, xt=X, status=status, y=times)

########################################################
##### Proposed Method #####
########################################################
source("ADS_sampler.R")
mcmcpar <- list(n.iter = 50000, n.burnin = 10000, n.thin = 10, 
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
                  gamprior = "default" # one of ("default", "custom", "uniform")
                  )

print(chains$PIP)
print(chains$post.prob[1:10,])  # Model PP estimated from the visited models
print(chains$hazPIP)

chains.pr = BVS.MCMC(modelpar, mcmcpar, verbose = T,
                  gprior.mat = "Product", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "default" # one of ("default", "custom", "uniform")
)

print(chains.pr$PIP)
print(chains.pr$post.prob[1:10,])  # Model PP estimated from the visited models
print(chains.pr$hazPIP)

chains.hyper = BVS.MCMC(modelpar, mcmcpar, verbose = T,
                  gprior.mat = "LCM-hyper", # one of ("LCM", "LCM-hyper", "Product")
                  gamprior = "default" # one of ("default", "custom", "uniform")
)

print(chains.hyper$PIP)
print(chains.hyper$post.prob[1:10,])  # Model PP estimated from the visited models
print(chains.hyper$hazPIP)

