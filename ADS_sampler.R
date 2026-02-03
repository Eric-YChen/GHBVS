source("utils/routines2.R")
source("utils/helper_functions.R")

get_move_prob <- function(gam, move_type) {
  haz <- gam2haz(gam)
  
  if (haz == "Null") {
    probs <- c("A_null" = 1.0)
  } else if (haz == "AH" || haz == "PH") {
    if (all(gam == 1) || all(gam == 2)) {
      probs <- c("A/D" = 0.5, "C" = 0.25, "C*" = 0.25)
    } else {
      probs <- c("A/D" = 0.5, "S" = 0.1, "C" = 0.2, "C*" = 0.2)
    }
  } else if (haz == "AFT") {
    if (all(gam == 4)) {
      probs <- c("A/D" = 0.6, "C*" = 0.4)
    } else {
      probs <- c("A/D" = 0.6, "S" = 0.2, "C*" = 0.2)
    }
  } else {  # GH
    if (all(gam %in% c(0, 3))) {
      if (all(gam == 3)) {
        probs <- c("A/D_GH" = 0.5, "C" = 0.25, "C*" = 0.25)
      } else {
        probs <- c("A/D_GH" = 0.5, "S" = 0.15, "C" = 0.15, "C*" = 0.2)
      }
    } else {
      probs <- c("A/D_GH" = 0.5, "S" = 0.25, "C" = 0.25)
    }
  }
  
  if (move_type %in% names(probs)) return(probs[move_type])
  return(0)
}

# Helper: get valid indices for A/D_GH
get_valid_idx <- function(gam) {
  p <- length(gam)
  tabgam <- table(factor(gam, levels = 0:4))
  
  if (tabgam["3"] == 1) {
    remaining <- gam[gam != 3]
    becomes_AH <- all(remaining %in% c(0,1)) && any(remaining == 1)
    becomes_PH <- all(remaining %in% c(0,2)) && any(remaining == 2)
    if (becomes_AH || becomes_PH) {
      return((1:p)[-which(gam == 3)])
    }
  } else if (tabgam["3"]==0 && tabgam["1"]>1 && tabgam["2"]==1) {
    return((1:p)[-which(gam == 2)])
  } else if (tabgam["3"]==0 && tabgam["2"]>1 && tabgam["1"]==1) {
    return((1:p)[-which(gam == 1)])
  } else if (tabgam["3"]==0 && tabgam["2"]==1 && tabgam["1"]==1) {
    return((1:p)[gam == 0])
  }
  return(1:p)
}

ADS.proposal.aft <- function(gam.old, gam.candid=NULL, ...){
  
  gam.candid <- gam.old
  haz <- gam2haz(gam.old)
  p <- length(gam.old)
  
  # MOVE SELECTION
  if (haz == "Null") {
    move <- "A_null"
  } else if (haz == "AH" || haz == "PH") {
    if (all(gam.old == 1) || all(gam.old == 2)) {
      move <- sample(c("A/D","C","C*"), 1, prob=c(0.5, 0.25, 0.25))
    } else {
      move <- sample(c("A/D","S","C","C*"), 1, prob=c(0.5, 0.1, 0.2, 0.2))
    }
  } else if (haz == "AFT") {
    if (all(gam.old == 4)) {
      move <- sample(c("A/D","C*"), 1, prob=c(0.6, 0.4))
    } else {
      move <- sample(c("A/D","S","C*"), 1, prob=c(0.6, 0.2, 0.2))
    }
  } else {  # GH
    if (all(gam.old %in% c(0,3))) {
      if (all(gam.old == 3)) {
        move <- sample(c("A/D_GH","C","C*"), 1, prob=c(0.5, 0.25, 0.25))
      } else {
        move <- sample(c("A/D_GH","S","C","C*"), 1, prob=c(0.5, 0.15, 0.15, 0.2))
      }
    } else {
      move <- sample(c("A/D_GH","S","C"), 1, prob=c(0.5, 0.25, 0.25))
    }
  }
  
  current_move <<- move

  # EXECUTE MOVE AND COMPUTE HASTINGS RATIO
  if (move == "A/D_GH") {
    valid_idx <- get_valid_idx(gam.old)
    n_valid <- length(valid_idx)
    idx <- sample(valid_idx, 1)
    
    if (gam.candid[idx] > 0) {
      # DELETE
      gam.candid[idx] <- 0
      current_move <<- "D_GH"
      new_haz <- gam2haz(gam.candid)
      
      # q_fwd = P(A/D_GH|old) × (1/n_valid)
      # q_rev depends on new model type
      if (new_haz == "Null") {
        # Reverse is A_null: q_rev = 1 × (1/p) × (1/4)
        q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid
        q_rev <- 1.0 / p / 4
      } else {
        # Reverse is A/D_GH add: q_rev = P(A/D_GH|new) × (1/n_valid_new) × (1/3)
        n_valid_new <- length(get_valid_idx(gam.candid))
        q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid
        q_rev <- get_move_prob(gam.candid, "A/D_GH") / n_valid_new / 3
      }
      log.q <- log(q_rev) - log(q_fwd)
      
    } else {
      # ADD
      new_val <- sample(1:3, 1)
      gam.candid[idx] <- new_val
      current_move <<- "A_GH"
      
      # q_fwd = P(A/D_GH|old) × (1/n_valid) × (1/3)
      # q_rev = P(A/D_GH|new) × (1/n_valid_new)
      n_valid_new <- length(get_valid_idx(gam.candid))
      q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid / 3
      q_rev <- get_move_prob(gam.candid, "A/D_GH") / n_valid_new
      log.q <- log(q_rev) - log(q_fwd)
    }
    
  } else if (move == "A/D") {
    idx <- sample(1:p, 1)
    
    if (gam.candid[idx] > 0) {
      # DELETE
      gam.candid[idx] <- 0
      current_move <<- "D"
      new_haz <- gam2haz(gam.candid)
      
      if (new_haz == "Null") {
        # Reverse is A_null
        q_fwd <- get_move_prob(gam.old, "A/D") / p
        q_rev <- 1.0 / p / 4
      } else {
        # Reverse is A/D
        q_fwd <- get_move_prob(gam.old, "A/D") / p
        q_rev <- get_move_prob(gam.candid, "A/D") / p
      }
      log.q <- log(q_rev) - log(q_fwd)
      
    } else {
      # ADD
      gam.candid[idx] <- c(1,2,4)[haz == c("AH","PH","AFT")]
      current_move <<- "A"
      
      # Reverse is A/D
      q_fwd <- get_move_prob(gam.old, "A/D") / p
      q_rev <- get_move_prob(gam.candid, "A/D") / p
      log.q <- log(q_rev) - log(q_fwd)
    }
    
  } else if (move == "C") {
    nonzero_idx <- which(gam.candid != 0)
    n_nonzero <- length(nonzero_idx)
    idx <- ifelse(n_nonzero == 1, nonzero_idx, sample(nonzero_idx, 1))
    gam.candid[idx] <- sample(setdiff(1:3, gam.candid[idx]), 1)
    
    # q_fwd = P(C|old) × (1/n_nonzero) × (1/2)
    # q_rev = P(C|new) × (1/n_nonzero) × (1/2)  [same n_nonzero]
    q_fwd <- get_move_prob(gam.old, "C") / n_nonzero / 2
    q_rev <- get_move_prob(gam.candid, "C") / n_nonzero / 2
    log.q <- log(q_rev) - log(q_fwd)
    
  } else if (move == "S") {
    idx1 <- sample(1:p, 1)
    gam.1 <- gam.candid[idx1]
    other_idx <- which(gam.candid != gam.1)
    n_other <- length(other_idx)
    idx2 <- sample(other_idx, 1)
    gam.2 <- gam.candid[idx2]
    
    # swap
    if (gam.1 + gam.2 != 3) {
      gam.candid[idx1] <- gam.2
      gam.candid[idx2] <- gam.1
    } else if (gam.1 %in% c(1, 2)) {
      if (sample(c(TRUE,FALSE), 1)) {
        gam.candid[idx1] <- 0
        gam.candid[idx2] <- 3
      } else {
        gam.candid[idx1] <- 3
        gam.candid[idx2] <- 0
      }
    } else {
      if (sample(c(TRUE,FALSE), 1)) {
        gam.candid[idx1] <- 1
        gam.candid[idx2] <- 2
      } else {
        gam.candid[idx1] <- 2
        gam.candid[idx2] <- 1
      }
    }
    
    n_other_new <- sum(gam.candid != gam.candid[idx1])
    q_fwd <- get_move_prob(gam.old, "S") / p / n_other
    q_rev <- get_move_prob(gam.candid, "S") / p / n_other_new
    
    # Account for the 1/2 factor in special cases
    if (gam.1 + gam.2 == 3) {
      q_fwd <- q_fwd / 2
    }
    gam.new.1 <- gam.candid[idx1]
    gam.new.2 <- gam.candid[idx2]
    if (gam.new.1 + gam.new.2 == 3) {
      q_rev <- q_rev / 2
    }
    
    log.q <- log(q_rev) - log(q_fwd)
    
  } else if (move == "C*") {
    
    if (haz == "AFT") {
      new_val <- sample(1:3, 1)
      gam.candid[gam.candid == 4] <- new_val
      new_haz <- gam2haz(gam.candid)
      
      # q_fwd = P(C*|AFT) × (1/3)
      q_fwd <- get_move_prob(gam.old, "C*") / 3
      
      if (new_haz == "GH") {
        # q_rev = P(C*|GH) × 1
        q_rev <- get_move_prob(gam.candid, "C*")
        current_move <<- "AFT_GH"
      } else {
        # q_rev = P(C*|AH or PH) × (1/2)
        q_rev <- get_move_prob(gam.candid, "C*") / 2
        current_move <<- ifelse(new_val == 1, "AFT_AH", "AFT_PH")
      }
      log.q <- log(q_rev) - log(q_fwd)
      
    } else if (haz == "GH") {
      gam.candid[gam.candid == 3] <- 4
      current_move <<- "GH_AFT"
      
      # q_fwd = P(C*|GH) × 1
      # q_rev = P(C*|AFT) × (1/3)
      q_fwd <- get_move_prob(gam.old, "C*")
      q_rev <- get_move_prob(gam.candid, "C*") / 3
      log.q <- log(q_rev) - log(q_fwd)
      
    } else {
      # AH or PH
      nonzero_idx <- which(gam.candid != 0)
      current_val <- gam.candid[nonzero_idx[1]]
      new_val <- sample(setdiff(c(1,2,4), current_val), 1)
      gam.candid[nonzero_idx] <- new_val
      new_haz <- gam2haz(gam.candid)
      
      # q_fwd = P(C*|AH/PH) × (1/2)
      q_fwd <- get_move_prob(gam.old, "C*") / 2
      
      if (new_haz == "AFT") {
        # q_rev = P(C*|AFT) × (1/3)
        q_rev <- get_move_prob(gam.candid, "C*") / 3
        current_move <<- ifelse(haz == "PH", "PH_AFT", "AH_AFT")
      } else {
        # q_rev = P(C*|AH/PH) × (1/2)
        q_rev <- get_move_prob(gam.candid, "C*") / 2
        current_move <<- ifelse(haz == "PH", "PH_AH", "AH_PH")
      }
      log.q <- log(q_rev) - log(q_fwd)
    }
    
  } else if (move == "A_null") {
    idx <- sample(1:p, 1)
    new_val <- sample(1:4, 1)
    gam.candid[idx] <- new_val
    new_haz <- gam2haz(gam.candid)
    
    # q_fwd = 1 × (1/p) × (1/4)
    q_fwd <- 1.0 / p / 4
    
    # q_rev depends on target structure
    if (new_haz == "GH") {
      # Reverse is A/D_GH delete
      n_valid_new <- length(get_valid_idx(gam.candid))
      q_rev <- get_move_prob(gam.candid, "A/D_GH") / n_valid_new
    } else {
      # Reverse is A/D delete
      q_rev <- get_move_prob(gam.candid, "A/D") / p
    }
    log.q <- log(q_rev) - log(q_fwd)
  }
  
  return(list(gam.proposal = gam.candid, lograt = log.q))
}

MH.step.fixed <- function(target.dist, proposal.dist, par, data,...){
  # initialize
  par.new = par.old = par
  postprob = 0

  # ensure the old model is in the dictionary (optimized: use environment for O(1) lookup)
  old.model.string = paste0(par$gam, collapse="")
  if (!exists(old.model.string, envir = model.dict, inherits = FALSE)){
    post.old = target.dist(par.old,data=data)
    assign(old.model.string, post.old, envir = model.dict)
  }else{
    post.old = get(old.model.string, envir = model.dict)
  }

  # propose new gamma
  step.new = proposal.dist(par$gam)
  gam.new = step.new$gam.proposal
  log.q = step.new$lograt
  par.new$gam = gam.new

  # check if new model is in dictionary (optimized: use environment for O(1) lookup)
  new.model.string = paste0(gam.new, collapse="")
  if (!exists(new.model.string, envir = model.dict, inherits = FALSE)){
    post.new = target.dist(par.new,data=data)
    assign(new.model.string, post.new, envir = model.dict)
  }else{
    post.new = get(new.model.string, envir = model.dict)
  }
  
  # accept probability
  log.acc.prob = post.new - post.old + log.q
  log.acc.prob = min(0, log.acc.prob)
  acc.prob = exp(log.acc.prob)
  if (is.nan(acc.prob)){
    acc.prob = 0
  }
  
  par.out = par
  if (runif(1) > acc.prob){  # reject
    par.out$gam = par.old$gam
    move.list[nrow(move.list)+1,] <<- list(current_move, 0)
    postprob = post.old
  }else{  # accept
    par.out$gam = par.new$gam
    move.list[nrow(move.list)+1,] <<- list(current_move, 1)
    postprob = post.new
  }
  return(list(par = par.out, acc = acc.prob, postprob = postprob))
}

BVS.MCMC <- function(model.par, mcmc.par,
                     gprior.mat = "LCM", gamprior = "default",
                     verbose = T, use_cpp = FALSE){
  ##### Useful Arguments #####
  # gamprior = c("unif", "default", "custom")
  # gprior.mat = c("LCM", "LCM-hyper", "product")
  #     "LCM-hyper": robust prior OR hyper g/n for g. To use the hyper g/n prior,
  #                 set 2 < a_g <= 4 in model.par$inits. Otherwise robust prior is used.
  # use_cpp = TRUE to use C++ optimizer (requires Rcpp, RcppArmadillo, RcppEigen, RcppNumerical)

  ##### C++ Setup #####
  USE_CPP <<- use_cpp
  if (use_cpp) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      warning("Rcpp not available, falling back to pure R")
      USE_CPP <<- FALSE
    } else {
      # Compile C++ code if not already loaded
      if (!exists("mle_gh_cpp")) {
        if (verbose) cat("Compiling C++ optimizer...\n")
        Rcpp::sourceCpp("src/optimizer_cpp.cpp")
      }
      if (!exists("lglik_cpp")) {
        Rcpp::sourceCpp("src/likelihood_functions.cpp")
      }
      if (verbose) cat("Using C++ optimizer\n")
    }
  }

  ##### Initialization #####
  inits <- model.par$inits
  data <- model.par$data
  inits$g <- nrow(data$x)*inits$g
  data$xt2 = t(data$xt)%*%data$xt; data$x2 = t(data$x)%*%data$x
  data$xt2.inv = solve(data$xt2); data$x2.inv = solve(data$x2)
  if (is.null(inits$n)){inits$n <- nrow(data$x)}
  if (is.null(inits$n.threshold)){inits$n.threshold <- inits$n*(1 - 0.5*(1-mean(data$status)))}
  
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
  gamprior <<- gamprior
  gam.param <<- model.par$gam.param
  #marg.liks <<- matrix(NA, nrow=0, ncol=3)
  #colnames(marg.liks) <<- c("model","MAP","marg.lik")
  gam.proposal = ADS.proposal.aft
  
  # Optimized: pre-allocate model.dict as an environment for O(1) lookup
  # Using environment instead of data.frame for faster key-value access
  model.dict <<- new.env(hash = TRUE, size = 10000L)
  model.dict.list <<- list()  # For final output conversion
  MH.step = MH.step.fixed
  if (gprior.mat == "LCM"){
    source("utils/EB_gprior_fixed.R")
  }else if (gprior.mat == "LCM-hyper"){
    source("utils/EB_gprior_robust.R")
  }else{
    source("utils/gprior_setup_fixed.R")
  }
  g.proposal = NULL
  acc.list = data.frame(matrix(NA, nrow=n.iter, ncol=n.chain))
  target.dist <- targetdist
  
  ##### MCMC #####
  for(chain in 1:n.chain){
    par = inits
    current.chain = data.frame(matrix(NA, nrow=n.save, ncol=length(par$gam)))
    current.haz = rep(NA, n.save)
    iter.save = 0
    for(iter in 1:n.iter){
      # update gamma #
      update = MH.step(target.dist, proposal.dist=gam.proposal, par, data)
      par.new = update$par
      acc.gam = update$acc

      # save acceptance rates
      acc.list[iter,chain] = acc.gam
      modelsize.df[iter,chain] = sum(par.new$gam!=0)
      postprob = update$postprob

      # Save chain
      if(iter > n.burn & iter %% n.thin == 0){
        iter.save = iter.save + 1
        current.chain[iter.save,] = par.new$gam
        current.haz[iter.save] = gam2haz(par.new$gam)
        postprob.list[iter.save,chain] = postprob
      }
      # Print progress
      if(verbose & iter %% printevery == 0){
        print(paste0("Chain:", chain, ", Iter:", iter, 
        ", Acc:", round(mean(acc.list[(iter-printevery+1):iter,chain]),3),
        ", modl.size: ", round(mean(modelsize.df[(iter-printevery+1):iter,chain]),3) ))
      }
      par = par.new
    }
    chain.list[[chain]] = current.chain
    haz.list[,chain] = current.haz
  }
  modelfreq = model.freq(chain.list)
  # modelprob = as.data.frame(unlist(model.dict))
  # colnames(modelprob) = "log-post"
  # model.accepted = names(unlist(model.dict)) %in% names(modelfreq)
  # modelprob$accept = 1*model.accepted
  
  # Convert environment back to data.frame for output
  model.dict.df <- data.frame(
    gam = ls(model.dict),
    postprob = sapply(ls(model.dict), function(x) get(x, envir = model.dict)),
    stringsAsFactors = FALSE
  )

  return(list(chain = chain.list,
              acc = acc.list, modelsize = modelsize.df,
              PIP = PIP(chain.list), hazPIP = haz.freq(haz.list),
              model.freq = modelfreq,
              model.prob = model.dict.df,
              post.prob = model.postprob(chain.list,postprob.list)))
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
