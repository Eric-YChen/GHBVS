print("current wd:")
print(getwd())
library(HazReg)
library(mvtnorm)
library(survival)
library(glmnet)
library(tictoc)
library(Rmpi)
library(doMPI)
library(doRNG)
library(Rcpp)
library(RcppArmadillo)

################################################################################
source("./utils/helper_functions.R")
n.vec = c(250,500,1000)
p = 10
n.simset = 250
censrat = c("0.25default","0.5default")
priorstr = "LCM"
filename = "mpi_simulation_raw.RData"

sigma = 0.7**abs(matrix(rep(c(1:p),p),nrow=p,byrow = T) - c(1:p))
hstr = c("PH","AFT","GH","AH")
baseline = c("LN","PGW")
#baseline = "PGW"
#hstr = "AH"


str.grid = expand.grid(n.vec, p, hstr, baseline, censrat)
#str.grid = str.grid[(str.grid[,3]!="GH" | str.grid[,4]!="PGW"),]
row.names(str.grid) = NULL
names(str.grid) = c("n","p","hstr","baseline")

true.models = apply(str.grid, 1, function(x){p = x[["p"]]; x = x[["hstr"]]
if (x == "AH"){tmod <- rep(0,p); tmod[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- 1; paste0(tmod,collapse = "")
}else if (x == "PH"){tmod <- rep(0,p); tmod[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- 2; paste0(tmod,collapse = "")
}else if (x == "AFT"){tmod <- rep(0,p); tmod[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- 4; paste0(tmod,collapse = "")
}else if (x == "GH"){tmod <- rep(0,p); tmod[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- 3; paste0(tmod,collapse = "")}})

str.grid = cbind(str.grid, cbind(true.models))
names(str.grid) = c("n","p","hstr","baseline","cens","true.mod")
################################################################################

#rm(list=ls()[!ls() %in% c("n.simset")])
rm(list=ls()[!ls() %in% c("n.simset","str.grid","sigma","filename","priorstr")])
mcmcpar = list(n.iter = 20000, n.burnin = 10000, n.thin = 2, 
               n.chain = 1, printevery = 10000)

# Create indices for all simulation runs
sim_indices <- expand.grid(i = 1:nrow(str.grid), j = 1:n.simset)

# Set up the MPI cluster
cat("Starting MPI cluster\n")
cl <- startMPIcluster()
# Get the number of workers in the MPI cluster
n_workers <- cl$workerCount
cat("Using", n_workers, "MPI workers for parallel processing\n")
registerDoMPI(cl)

set.seed(1234)
tic(msg="###########################\nTotal time elapsed")
# Use foreach with doRNG for proper parallel random number generation
results_parallel <- foreach(idx = 1:nrow(sim_indices), 
                            .packages = c("mvtnorm", "survival", "glmnet", "tictoc","RcppArmadillo","Rcpp"),
                            .combine = 'list', 
                            .options.RNG = 1234) %dorng% {
                              # Extract indices
                              i <- sim_indices$i[idx]
                              j <- sim_indices$j[idx]
                              
                              # Source required files - each worker needs these
                              source("./ADS_sampler.R")
                              source("./HazReg/routines.R")
                              
                              # Generate data on this worker instead of using pre-generated data
                              dat_full = simdat(n = str.grid[i,1], p = str.grid[i,2], sigma = sigma,
                                                hstr = str.grid[i,3], baseline = str.grid[i,4],
                                                n.simset = 1, cens = str.grid[i,5], seeds = j)
                              dat = dat_full$data[[1]]
                              cens_rat = dat_full$cens.rat[1]
                              
                              n = nrow(dat$x)
                              pxt = ncol(dat$xt)
                              px = ncol(dat$x)
                              
                              # Initialize (timing is done per worker)
                              start_time_init <- Sys.time()
                              initgam = init.gam(dat, history.keep = F)
                              end_time_init <- Sys.time()
                              init.runtime <- as.numeric(difftime(end_time_init, start_time_init, units = "secs"))
                              
                              modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1),
                                                           theta=rnorm(pxt,0,0.5), eta=rnorm(px,0,0.5),
                                                           g.theta = 1, g.eta = 1, g.aft = 1, g = 1,
                                                           gam = initgam, a_g = NULL,
                                                           nu.alpha = 0.01, nu.beta = 0.01), 
                                              # gam.param = list(layer1.betabin.par = c(1,1),
                                              #                  layer2.mulnom.par = if (str.grid[i,3]=="PH"){c(0,1,0,0)}else{c(0,0,0,1)},
                                              #                  layer3.beta.par = 1/3),
                                              data = dat)
                              
                              # Run sampler
                              start_time_sampler <- Sys.time()
                              tryCatch({chains = BVS.MCMC(modelpar, mcmcpar, verbose = F,
                                                          gprior.mat = priorstr, gamprior = "default")
                              end_time_sampler <- Sys.time()
                              sampler_runtime <- as.numeric(difftime(end_time_sampler, start_time_sampler, units = "secs"))
                              
                              # Process chains
                              fullchain = data.frame(gam=apply(chains$chain[[1]], 1, function(x){paste0(x, collapse="")}))
                              fullchain = cbind(fullchain, 
                                                acc.gam=chains$acc[seq(mcmcpar$n.burnin+mcmcpar$n.thin, 
                                                                       mcmcpar$n.iter, mcmcpar$n.thin),])
                              fullchain = cbind(fullchain, post.prob = chains$post.prob)
                              
                              # Return results for this simulation
                              list(
                                i = i,
                                j = j,
                                #pip = chains$PIP,
                                #freq = chains$model.freq,
                                #hazpip = chains$hazPIP,
                                fullchain = fullchain,
                                cens_rat = cens_rat,
                                runtime = c(init = init.runtime, sampler = sampler_runtime)
                              )},
                              error = function(e) {
                                cat(paste0(e$message,"\nn=",str.grid[i,1], ", p=",str.grid[i,2]," ", 
                                           str.grid[i,3],"_",str.grid[i,4]))
                                return(paste0(e$message,"\nn=",str.grid[i,1], ", p=",str.grid[i,2]," ", 
                                              str.grid[i,3],"_",str.grid[i,4]))
                              })
                            }

# Stop MPI cluster
closeCluster(cl)

time.elaps = toc()
time.elaps = time.elaps$toc-time.elaps$tic

# Save simulation results directly after foreach
save(str.grid, results_parallel, n.simset, time.elaps, file = filename)

mpi.quit()

