source("utils/routines2.R")
source("utils/helper_functions.R")

ADS.proposal.aft <- function(gam.old, gam.candid=NULL,...){
  # 0: (0,0), 1: (1,0)/t, 2: (0,1)/p, 3: (1,1)/tp, 4: AFT
  # (theta,eta) / (alpha,beta) / (time,PH) / (xt,x)
  
  gam.candid = gam.old
  haz = gam2haz(gam.old)
  if (haz == "Null"){  ## Null model ##
    move = sample(c("A/D","A/D_aft"),1,prob=c(3/4, 1/4))
  }else if (haz == "AH" | haz == "PH"){
    if (all(gam.old == 1) | all(gam.old == 2)){  ## Full AH/PH ##
      move = sample(c("A/D","C","C*"),1,prob=c(0.6, 0.2, 0.2))
    }else{
      move = sample(c("A/D","S","C","C*"),1,prob=c(0.4, 0.2, 0.2, 0.2))
    }
  }else if (haz == "AFT"){
    if (all(gam.old == 4)){  ## Full AFT ##
      move = sample(c("A/D","C*"),1,prob=c(2/3, 1/3))
    }else{
      move = sample(c("A/D","S","C*"),1,prob=c(0.4, 0.2, 0.2))
    }
  }else{
    if (all(gam.old %in% c(0,3))){ ## GH with only 0 & 3 ##
      if (all(gam.old == 3)){
        move = sample(c("A/D_GH","C","C*"),1,prob=c(0.4, 0.4, 0.2))
      }else{
        move = sample(c("A/D_GH","S","C","C*"),1,prob=c(0.4, 0.2, 0.2, 0.2))
      }
    }else{  ## mixed GH ##
      move = sample(c("A/D_GH","S","C"),1,prob=c(0.5, 0.25, 0.25))
    }
  }
  
  current_move <<- move
  # propose move
  if (move == "A/D_GH"){  # Add or Delete Move (between 0 and >0, model size +1 or -1)
    idx = sample(1:length(gam.old), 1)
    if (gam.candid[idx] > 0){
      gam.candid[idx] = 0
      current_move <<- "D_GH"
      log.q = -log(3)
    }else{
      gam.candid[idx] = sample(c(1:3),1)
      current_move <<- "A_GH"
      log.q = log(3)
    }
  }else if (move == "A/D"){  # Add or Delete Move for non-GH model (Naive ADS)
    idx = sample(1:length(gam.old), 1)
    if (gam.candid[idx] > 0){
      gam.candid[idx] = 0
      current_move <<- "D"
      log.q = 0
    }else{
      gam.candid[idx] = c(1,2,4)[haz == c("AH","PH","AFT")]
      current_move <<- "A"
      log.q = 0
    }
  }else if (move == "C"){  # Change Move (between 1,2,3, model size unchanged)
    if (length(which(gam.candid != 0)) == 1){
      idx = which(gam.candid != 0)
      gam.candid[idx] = sample(c(1:3)[-gam.candid[idx]],1)
    }else{
      idx = sample(which(gam.candid != 0), 1)
      gam.candid[idx] = sample(c(1:3)[-gam.candid[idx]],1)
    }
    log.q = 0
  }else if (move == "S"){  # Swap Move (between two variables, might change model size)
    idx1 = sample(c(1:length(gam.candid)),1)
    gam.1 = gam.candid[idx1]
    q_prop = sum(gam.candid != gam.1)
    idx2 = sample(c(1:length(gam.candid))[gam.candid!= gam.1], 1)
    gam.2 = gam.candid[idx2]
    if (gam.1+gam.2 != 3){
      gam.candid[idx1] = gam.2
      gam.candid[idx2] = gam.1
    }else if (gam.1 == 1 | gam.1 == 2){
      a = sample(c(T,F),1)
      if (a){
        gam.candid[idx1] = 0
        gam.candid[idx2] = 3
      }else{
        gam.candid[idx1] = 3
        gam.candid[idx2] = 0
      }
    }else if (gam.1 == 0 | gam.1 == 3){
      a = sample(c(T,F),1)
      if (a){
        gam.candid[idx1] = 1
        gam.candid[idx2] = 2
      }else{
        gam.candid[idx1] = 2
        gam.candid[idx2] = 1
      }
    }
    q_reverse = sum(gam.candid != gam.candid[idx1])
    log.q = log(q_reverse) - log(q_prop)
  }else if (move == "C*"){  # Change Move: changing model type
    if (haz == "AFT"){
      new.mod <- sample(c(1:3),1)
      gam.candid[gam.candid == 4] = new.mod
      if (gam2haz(gam.candid) == "GH"){
        log.q = log(3)
        current_move <<- "AFT_GH"
      }else{
        log.q = log(1.5)
        current_move <<- ifelse(new.mod == 1, "AFT_AH", "AFT_PH")
      }
    }else if (haz == "GH"){
      gam.candid[gam.candid == 3] <- 4
      log.q = -log(3)
      current_move <<- "GH_AFT"
    }else{
      idx = which(gam.candid != 0)
      gam.candid[idx] = sample(c(1,2,4)[-gam.candid[idx]],1)
      if (gam2haz(gam.candid) == "AFT"){
        log.q = -log(1.5)
        current_move <<- ifelse(haz == "PH", "PH_AFT", "AH_AFT")
      }else{
        log.q = 0
        current_move <<- ifelse(haz == "PH", "PH_AH", "AH_PH")
      }
    }
  }
  return(list(gam.proposal = gam.candid, lograt = log.q))
}

BVS.MCMC <- function(model.par, mcmc.par, 
                     gprior.mat = "LCM",gamprior = "default", 
                     verbose = T){
  ##### Useful Arguments #####
  # gamprior = c("unif", "default", "custom")
  # gprior.mat = c("EB", "EB-hyper", "FB"), [requires g.sampler = NULL & g.fixed = T]
  #     "EB-hyper": robust prior OR hyper g/n for g. To use the hyper g/n prior,
  #                 set 2 < a_g <= 4 in model.par$inits. Otherwise robust prior is used.
  
  ##### Initialization #####
  inits <- model.par$inits
  data <- model.par$data
  inits$g <- nrow(data$x)*inits$g
  data$xt2 = t(data$xt)%*%data$xt; data$x2 = t(data$x)%*%data$x
  data$xt2.inv = solve(data$xt2); data$x2.inv = solve(data$x2)
  gam.sampler = "ADS.aft"; g.sampler = "MH"; gprior = "robust"; g.fixed = T 
  if (gprior.mat == "LCM"){gprior.mat = "EB"
  }else if(gprior.mat == "Product"){gprior.mat = "FB"
  }else if(gprior.mat == "LCM-hyper"){gprior.mat = "EB-hyper"
  }else{stop("Invalid gprior type")}
  if (g.fixed){g.sampler = NULL}
  if (is.null(inits$n)){inits$n <- nrow(data$x)}
  if (is.null(inits$n.threshold)){inits$n.threshold <- inits$n*(1 - 0.5*(1-mean(data$status)))}
  #target.dist <- model.par$target.dist
  
  # mcmc.par(MCMC parameters): list(n.iter, n.burn, n.thin, n.save, n.print, n.chain)
  n.iter <- mcmc.par$n.iter; n.burn <- mcmc.par$n.burn; n.thin <- mcmc.par$n.thin
  if (n.thin < 1){stop("n.thin must be integer >= 1")}
  if (n.iter < n.burn){stop("n.iter must be greater than n.burn")}
  if (n.iter%%n.thin != 0){stop("n.iter must be divisible by n.thin")}
  printevery <- mcmc.par$printevery; n.chain <- mcmc.par$n.chain; n.save = floor((n.iter-n.burn)/n.thin)
  
  move.list <<- data.frame(matrix(NA, nrow=0, ncol=2))
  colnames(move.list) <<- c("move","accept")
  chain.list = list()
  haz.list = matrix(NA, nrow=n.save, ncol=n.chain)
  modelsize.df = data.frame(matrix(NA, nrow=n.iter, ncol=n.chain))
  postprob.list = matrix(NA, nrow=n.save, ncol=n.chain)
  gprior <<- gprior
  gamprior <<- gamprior
  gam.param <<- model.par$gam.param
  #marg.liks <<- matrix(NA, nrow=0, ncol=3)
  #colnames(marg.liks) <<- c("model","MAP","marg.lik")

  ##### Choosing sampler for gamma #####
  gam.proposal = ADS.proposal.aft
  
  ##### Choosing sampler for g #####
  if (is.null(g.sampler)){  # if g is not updated iteratively
    source("utils/laplace_g_sampler.R")
    if (g.fixed){
      model.dict <<- data.frame(matrix(NA, nrow=0, ncol=2))
      colnames(model.dict) <<- c("gam","postprob")
      MH.step = MH.step.fixed
      if (gprior.mat == "EB"){
        source("utils/EB_gprior_fixed.R")
      }else if (gprior.mat == "EB-hyper"){
        source("utils/EB_gprior_robust.R")
      }else{
        source("utils/gprior_setup_fixed.R")
      }
    }else{
      model.dict <<- data.frame(matrix(NA, nrow=0, ncol=4))
      colnames(model.dict) <<- c("gam","g.theta","g.eta","postprob")
      MH.step = MH.step.laplace
      source("utils/gprior_setup_laplace.R")
      g.list = list()
    }
    g.proposal = NULL
    acc.list = data.frame(matrix(NA, nrow=n.iter, ncol=n.chain))
    
  }else{  # if updating g, OBSOLETE #
    model.dict <<- data.frame(matrix(NA, nrow=0, ncol=6))
    colnames(model.dict) <<- c("gam","g.theta","g.eta","g.aft","postprob","name")
    lr.t <<- inits$lr.t; lr.h <<- inits$lr.h; lr.a <<- inits$lr.a
    t.steq.seq <<- 0; h.steq.seq <<- 0; a.steq.seq <<- 0
    lr.mat <<- matrix(NA, nrow=0, ncol=3)
    colnames(lr.mat) <<- c("g.theta","g.eta","g.aft")
    acc.list = array(NA, dim=c(n.iter,n.chain,4), dimnames = list(NULL,NULL,c("gam","gt","ge","ga")))
    g.list = list()
    
    source("utils/gibbs_g_sampler.R")
    source("utils/gprior_setup_gibbs.R")
    if (g.sampler == "MALA"){
      hess.list <<- list()
      g.proposal = MALA.sampler
      MH.step = MH.step.MALA
      g.acc.current <<- c(0.574,0.574,0.574)
    }else if (g.sampler == "MH"){
      MH.step = MH.step.ASM
      g.proposal = MH.sampler
      g.acc.current <<- c(0.234,0.234,0.234)
    }
  }
  target.dist <- targetdist
  #fisher.mat <<- compute_fisher_info(inits,data)
  
  ##### MCMC #####
  for(chain in 1:n.chain){
    par = inits
    current.chain = data.frame(matrix(NA, nrow=n.save, ncol=length(par$gam)))
    current.haz = rep(NA, n.save)
    if (!(g.fixed)){
      current.gs = data.frame(matrix(NA, nrow=n.save, ncol=3))
      colnames(current.gs) = c("g.theta","g.eta","g.aft")
    }
    iter.save = 0
    for(iter in 1:n.iter){
      # update gamma #
      update = MH.step(target.dist, proposal.dist=gam.proposal, par, data)
      par.new = update$par
      acc.gam = update$acc
      
      if (!is.null(g.proposal)){ ### OBSOLETE, not in use ###
        if (g.sampler == "MALA"){
          # update g.theta #
          if (any(par.new$gam==1) | any(par.new$gam==3)){
            update = g.step(target.dist, proposal.dist=g.proposal, par.new, data, type = "t")
            par.new = update$par
            acc.gt = update$acc
            g.acc.current[1] <<- acc.gt
          }else{
            g.prop = abs(rt(1, df = par.new$gt.scale))
            par.new$g.theta = g.prop
            acc.gt = 1
          }
          # update g.eta #
          if (any(par.new$gam==2) | any(par.new$gam==3)){
            update = g.step(target.dist, proposal.dist=g.proposal, par.new, data, type = "h")
            par.new = update$par
            acc.ge = update$acc
            g.acc.current[2] <<- acc.ge
          }else{
            g.prop = abs(rt(1, df = par.new$ge.scale))
            par.new$g.eta = g.prop
            acc.ge = 1
          }
        }else if (g.sampler == "MH"){
          # update g.theta #
          if (any(par.new$gam==1) | any(par.new$gam==3)){
            update = g.proposal(target.dist, par.new, data, type = "t")
            par.new = update$par
            acc.gt = update$acc
            g.acc.current[1] <<- acc.gt
          }else{
            # NEED CHANGE HERE
            if (gprior == "t"){
              g.prop = abs(rt(1, df = par.new$gt.scale))
              par.new$g.theta = g.prop
              acc.gt = NA
            }else if (gprior == "robust"){
              g.prop = robust.g(n = nrow(data$xt), p = ncol(data$xt))
              par.new$g.theta = g.prop
              acc.gt = NA
            }
          }
          # update g.eta #
          if (any(par.new$gam==2) | any(par.new$gam==3)){
            update = g.proposal(target.dist, par.new, data, type = "h")
            par.new = update$par
            acc.ge = update$acc
            g.acc.current[2] <<- acc.ge
          }else{
            if (gprior == "t"){
              g.prop = abs(rt(1, df = par.new$ge.scale))
              par.new$g.eta = g.prop
              acc.ge = NA
            }else if (gprior == "robust"){
              g.prop = robust.g(n = nrow(data$x), p = ncol(data$x))
              par.new$g.eta = g.prop
              acc.ge = NA
            }
          }
          # update g.aft #
          if (any(par.new$gam == 4)){
            update = g.proposal(target.dist, par.new, data, type = "a")
            par.new = update$par
            acc.ga = update$acc
            g.acc.current[3] <<- acc.ga
          }else{
            if (gprior == "t"){
              g.prop = abs(rt(1, df = par.new$ga.scale))
              par.new$g.aft = g.prop
              acc.ga = NA
            }else if (gprior == "robust"){
              g.prop = robust.g(n = nrow(data$xt), p = ncol(data$xt))
              par.new$g.aft = g.prop
              acc.ga = NA
            }
          }
        }
        lr.mat <<- rbind(lr.mat, c(lr.t, lr.h, lr.a))
        #print(c(par.new$g.theta,par.new$g.eta))
      }
      # save acceptance rates
      if (!is.null(g.proposal)){
        acc.list[iter,chain,] = c(acc.gam, acc.gt, acc.ge, acc.ga)
      }else{
        acc.list[iter,chain] = acc.gam
      }
      modelsize.df[iter,chain] = sum(par.new$gam!=0)
      postprob = update$postprob

      # Save chain
      if(iter > n.burn & iter %% n.thin == 0){
        iter.save = iter.save + 1
        current.chain[iter.save,] = par.new$gam
        current.haz[iter.save] = gam2haz(par.new$gam)
        postprob.list[iter.save,chain] = postprob
        if (!is.null(g.proposal)){
          current.gs[iter.save,] = c(par.new$g.theta, par.new$g.eta, par.new$g.aft)
        }else{
          if (!g.fixed){
            if (gam2haz(par.new$gam) == "AFT"){
              g.aft.new = as.numeric(model.dict[which(model.dict$gam == paste0(par.new$gam, collapse="")), 2])
              current.gs[iter.save,] = c(NA, NA, g.aft.new)
            }else{
              g.theta.new = as.numeric(model.dict[which(model.dict$gam == paste0(par.new$gam, collapse="")), 2])
              g.eta.new = as.numeric(model.dict[which(model.dict$gam == paste0(par.new$gam, collapse="")), 3])
              current.gs[iter.save,] = c(g.theta.new, g.eta.new, NA)
            }
          }
        }
      }
      # Print progress
      if(verbose & iter %% printevery == 0){
        if (!is.null(g.proposal)){
          print(paste0("Chain:", chain, 
                       ", Iter:", iter, 
                       ", Acc:(", 
                       round(mean(acc.list[(iter-printevery+1):iter,chain,1]),2), " ", 
                       round(mean(acc.list[(iter-printevery+1):iter,chain,2],na.rm=T),2), " ",
                       round(mean(acc.list[(iter-printevery+1):iter,chain,3],na.rm=T),2), " ",
                       round(mean(acc.list[(iter-printevery+1):iter,chain,4],na.rm=T),2),")",
                       ", g: (", round(par.new$g.theta,3), " ", round(par.new$g.eta,3), " ", round(par.new$g.aft,3), 
                       "), lr: (", round(lr.t,3), " ", round(lr.h,3), " ", round(lr.a,3),
                       "), modl.size: ", round(mean(modelsize.df[(iter-printevery+1):iter,chain]),3) ))
        }else{
          print(paste0("Chain:", chain, ", Iter:", iter, 
                       ", Acc:", round(mean(acc.list[(iter-printevery+1):iter,chain]),3),
                       ", modl.size: ", round(mean(modelsize.df[(iter-printevery+1):iter,chain]),3) ))
        }
      }
      par = par.new
    }
    chain.list[[chain]] = current.chain
    haz.list[,chain] = current.haz
    if (!(g.fixed)){
      g.list[[chain]] = current.gs
    }
  }
  modelfreq = model.freq(chain.list)
  # modelprob = as.data.frame(unlist(model.dict))
  # colnames(modelprob) = "log-post"
  # model.accepted = names(unlist(model.dict)) %in% names(modelfreq)
  # modelprob$accept = 1*model.accepted
  
  # removing MAP.pars and itervec from global environment
  # rm(MAP.pars, itervec, model.dict)
  
  if (g.fixed){
    return(list(chain = chain.list,
                acc = acc.list, modelsize = modelsize.df, 
                PIP = PIP(chain.list), hazPIP = haz.freq(haz.list), 
                model.freq = modelfreq, 
                model.prob = model.dict,
                post.prob = model.postprob(chain.list,postprob.list)))
  }
  else{
    return(list(chain = chain.list, g.chain = g.list, 
                acc = acc.list, modelsize = modelsize.df, 
                PIP = PIP(chain.list), hazPIP = haz.freq(haz.list),
                model.freq = modelfreq, 
                model.prob = model.dict,
                post.prob = model.postprob(chain.list,postprob.list)))
  }
}

prior.gam = function(par4, log=T, type = gamprior, 
                     parameters = gam.param){
  # Default parameters:
  # list(layer1.betabin.par = c(1,1), 
  #      layer2.mulnom.par = c(0.25,0.25,0.25,0.25), # (AH,PH,GH,AFT)
  #      layer3.beta.par = 1/3) # closer to 1 means more sparsity
  
  # exclude models with p > n.threshold
  if (sum(par4$gam!=0) > par4$n.threshold){
    if (log) return(-Inf)
    else return(0)
  }
  
  if (type == "unif"){
    ifelse(log, return(0), return(1))
  }else{
    gam = par4$gam; p = length(gam)
    p1 = sum(gam==1); p2 = sum(gam==2); p3 = sum(gam==3); p4 = sum(gam==4); p0 = sum(gam==0)
    curr.model = gam2haz(gam); m = p-p0
    
    if (type == "default"){
      if (curr.model == "GH"){
        pi.gam.layer2 = - log(3^m - 2)
        #pi.gam.layer2 = 0
        if (m==1){
          pi.gam.layer3 = 0
        }else{
          layer3.sum = sum(sapply(c(0:m), function(x){
            1/(dbinom(x, m, 1/3) * ifelse(x==0, 1-2^{1-m}, 1))
          }))
          pi.gam.layer3 = -dbinom(p3, m, 1/3, log=T) - ifelse(p3==0,log(1-2^(1-m)),0) - log(layer3.sum)
        }
        
        pi.gam = pi.gam.layer2 + pi.gam.layer3
      }else{
        pi.gam = 0
      }
      pi.gam = pi.gam + lbeta(1+m, 1+p0)
    }else if (type == "custom"){
      layer1.betabin.par = parameters$layer1.betabin.par
      layer2.mulnom.par = parameters$layer2.mulnom.par
      layer3.beta.par = parameters$layer3.beta.par
      
      if (curr.model == "GH"){
        pi.gam.layer2 = - log(3^m - 2) + log(layer2.mulnom.par[3])
        
        if (m==1){
          pi.gam.layer3 = 0
        }else{
          layer3.sum = sum(sapply(c(0:m), function(x){
            1/(dbinom(x, m, 1/3) * ifelse(x==0, 1-2^{1-m}, 1))
          }))
          pi.gam.layer3 = -dbinom(p3, m, 1/3, log=T) - ifelse(p3==0,log(1-2^(1-m)),0) - log(layer3.sum)
        }
        
        pi.gam = pi.gam.layer2 + pi.gam.layer3
      }else{
        if (curr.model == "AFT"){
          pi.gam.layer2 = log(layer2.mulnom.par[4])
        }else if (curr.model == "PH"){
          pi.gam.layer2 = log(layer2.mulnom.par[2])
        }else if (curr.model == "AH"){
          pi.gam.layer2 = log(layer2.mulnom.par[1])
        }
        pi.gam = pi.gam.layer2
      }
      pi.gam = pi.gam + lbeta(layer1.betabin.par[1] + m, layer1.betabin.par[2] + p0)
    }
    ifelse(log, return(pi.gam), return(exp(pi.gam)))
  }
}

# gh.multpl <- function(p, k){
#   if (!is.numeric(p) || !is.numeric(k) || p <= 0 || k < 0) {
#     stop("p must be a positive integer and k must be non-negative")
#   }
#   
#   # Convert to integers if they aren't already
#   p <- as.integer(p)
#   k <- as.integer(k)
#   
#   # Check if k is achievable (sum can be at most 2p)
#   if (k > 2 * p) {
#     return(0)
#   }
#   
#   # Calculate range for n_d (number of 'd' letters)
#   n_d_min <- max(0, k - p)
#   n_d_max <- min(p, floor(k / 2))
#   
#   # If no valid n_d values exist
#   if (n_d_min > n_d_max) {
#     return(0)
#   }
#   
#   # Calculate the sum over all valid n_d values
#   total_combinations <- 0
#   
#   for (n_d in n_d_min:n_d_max) {
#     # Calculate n_a, and the constraint k - 2*n_d for n_b + n_c
#     n_a <- p - k + n_d
#     sum_bc <- k - 2 * n_d
#     
#     # Calculate the multinomial coefficient term
#     # Formula: (p! * 2^(k-2*n_d)) / ((p-k+n_d)! * n_d! * (k-2*n_d)!)
#     term <- (factorial(p) * (2^sum_bc)) / 
#       (factorial(n_a) * factorial(n_d) * factorial(sum_bc))
#     
#     total_combinations <- total_combinations + term
#   }
#   
#   return(as.integer(total_combinations))
# }

