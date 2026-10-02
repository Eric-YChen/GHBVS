library(mvtnorm)
library(survival)
library(glmnet)
library(HazReg)
library(parallel)

set.seed(1234)
nki70 <- read.csv("nki70.csv", header = TRUE, stringsAsFactors = T)
X <- nki70[,-c(1:2)]
X$Diam <- as.numeric(X$Diam) - 1
X$N <- as.numeric(X$N)
X$ER <- as.numeric(X$ER) - 1
X$Grade <- as.numeric(X$Grade)
X <- apply(X,2,scale)
dat <- list(y = nki70$time, status = nki70$event, x=X, xt=X)


################ Comparison 1 ##########################
##### mombf: BVS for AFT model with non-local prior #####
########################################################
library(mombf)
library(foreach)
library(doParallel)
library(doRNG)

ys <- Surv(log(nki70$time), event=nki70$event)
df <- data.frame(ys = ys, as.data.frame(X), time = nki70$time, status = nki70$event)

f1y <- formula(paste('ys ~ ',paste(colnames(X), collapse='+'),sep=''))
priorCoef= zellnerprior(tau=nrow(X))
priorDelta= modelbbprior(1,1)

set.seed(1234)
ms1.list <- list()
time.ms1 <- Sys.time()
for(i in 1:10){
  ms1 <- modelSelection(f1y, data = df,
                        priorCoef = priorCoef, priorDelta = priorDelta,
                        enumerate = FALSE, method = 'Laplace', niter = 110000, 
                        thinning = 4, burnin = 10000,
                        center = FALSE, scale = FALSE)
  ms1.list[[i]] <- ms1
}
time.ms1 <- Sys.time() - time.ms1
# combine chains and extract top models
a <- Reduce(rbind,lapply(ms1.list, function(x){cbind(x$postSample, prob=x$postProb)}))
b <- data.frame(var=apply(a,1,function(x){paste0(x[1:76],collapse = "")}),prob=a[,77])
cc <- b[!duplicated(b[,1]),]
cc <- cc[order(cc$prob,decreasing = T),]
cc$postprob <- exp(cc$prob - mclust::logsumexp(cc$prob))
cc <- cc[order(cc$postprob,decreasing = T),]
head(cc,10)

################ Comparison 2 ##########################
##### BVSNLP: Cox PiMOM prior #####
########################################################
library(BVSNLP)
y2 <- Surv(time = df$time, df$status)

# Parallelize bvs across 10 cores
numWorkers <- min(10, max(1, parallel::detectCores()))
cl <- parallel::makeCluster(numWorkers, type = "PSOCK")
parallel::clusterEvalQ(cl, {
  library(BVSNLP)
  NULL
})
parallel::clusterExport(cl, varlist = c("X", "df"), envir = environment())

time.bvsnlp <- Sys.time()
bvs.list <- parallel::parLapply(cl, 1:10, function(i) {
  bvs(X = data.frame(X), resp = as.matrix(cbind(df$time, df$status)), 
      family = "survival", nlptype = "pMOM", mod_prior = "unif",
      niter = 11000, inseed = i, ncpu = 1)
})
time.bvsnlp <- Sys.time() - time.bvsnlp
parallel::stopCluster(cl)

# combine chains and extract top models
b <- Reduce(rbind,lapply(bvs.list, function(x){
  mods <- Reduce(c,lapply(x$max_models, function(y){paste0(t(y),collapse = ",")}))
  mods <- data.frame(var = mods, prob = x$max_prob_vec)
  return(mods)
}))
cc <- b[!duplicated(b[,1]),]
cc$postprob <- exp(cc$prob - mclust::logsumexp(cc$prob))
cc <- cc[order(cc$postprob,decreasing = T),]
head(cc,10)
########################################################

# Save the results
save(ms1, bvsfit, time.ms1, time.bvsnlp, 
     file = "nki70_results.RData")
