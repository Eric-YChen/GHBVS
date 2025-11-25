MH.step.laplace <- function(target.dist, proposal.dist, par, data,...){
  # initialize
  par.new = par.old = par
  postprob = 0
  
  # ensure the old model is in the dictionary
  old.model.string = paste0(par$gam, collapse="")
  if (!(old.model.string %in% model.dict[,1])){
    post.old = target.dist(par.old,data=data)
    if (post.old == -Inf){
      old.gt = NA
      old.ge = NA
    }else{
      # if (!exists("MAP.pars")){
      #   cat(paste0("p=",length(par$gam),", n=",par$n,"\n", str.grid[[i,3]][[j]],str.grid[[i,4]][[j]]))
      #   save(.Random.seed, file = "error_random_seed.RData")
      # }
      old.gt = if (gprior == "t") {
        ifelse(length(MAP.pars$g.theta)==0,NA,exp(MAP.pars$g.theta))
      }else if (gprior == "robust"){
        ifelse(length(MAP.pars$g.theta)==0,NA,exp(MAP.pars$g.theta)+(1+par$n)/(1+sum(MAP.pars$gam==1 | MAP.pars$gam==3| MAP.pars$gam==4))-1)}
      old.ge = if (gprior == "t") {
        ifelse(length(MAP.pars$g.eta)==0,NA,exp(MAP.pars$g.eta))
      }else if (gprior == "robust"){
        ifelse(length(MAP.pars$g.eta)==0,NA,exp(MAP.pars$g.eta)+(1+par$n)/(1+sum(MAP.pars$gam==2 | MAP.pars$gam==3))-1)}
    }
    model.dict[nrow(model.dict)+1,] <<- c(old.model.string, old.gt, old.ge, post.old)
  }else{
    idx = which(model.dict[,1] == old.model.string)
    post.old = as.numeric(model.dict[idx,4])
  }
  
  # propose new gamma
  step.new = proposal.dist(par$gam)
  gam.new = step.new$gam.proposal
  log.q = step.new$lograt
  par.new$gam = gam.new
  
  new.model.string = paste0(gam.new, collapse="")
  if (!(new.model.string %in% model.dict[,1])){
    post.new = target.dist(par.new,data=data)
    if (post.new == -Inf){
      new.gt = NA
      new.ge = NA
    }else{
      # if (!exists("MAP.pars")){
      #   cat(paste0("p=",length(par$gam),", n=",par$n,"\n", str.grid[[i,3]][[j]],str.grid[[i,4]][[j]]))
      #   save(.Random.seed, file = "error_random_seed.RData")
      # }
      new.gt = if (gprior == "t") {
        ifelse(length(MAP.pars$g.theta)==0,NA,exp(MAP.pars$g.theta))
      }else if (gprior == "robust"){
        ifelse(length(MAP.pars$g.theta)==0,NA,exp(MAP.pars$g.theta)+(1+par$n)/(1+sum(MAP.pars$gam==1 | MAP.pars$gam==3| MAP.pars$gam==4))-1)}
      new.ge = if (gprior == "t") {
        ifelse(length(MAP.pars$g.eta)==0,NA,exp(MAP.pars$g.eta))
      }else if (gprior == "robust"){
        ifelse(length(MAP.pars$g.eta)==0,NA,exp(MAP.pars$g.eta)+(1+par$n)/(1+sum(MAP.pars$gam==2 | MAP.pars$gam==3))-1)}
    }
    model.dict[nrow(model.dict)+1,] <<- c(new.model.string, new.gt, new.ge, post.new)
  }else{
    idx = which(model.dict[,1] == new.model.string)
    post.new = as.numeric(model.dict[idx,4])
  }
  
  # accept probability
  log.acc.prob = post.new - post.old + log.q
  log.acc.prob = min(0, log.acc.prob)
  acc.prob = exp(log.acc.prob)
  if (is.na(acc.prob)){
    acc.prob = 0
  }
  
  par.out = par
  if (runif(1) > acc.prob){  # reject
    par.out$gam = par.old$gam
    move.list[nrow(move.list)+1,] <<- c(current_move, 0)
    postprob = post.old
  }else{  # accept
    par.out$gam = par.new$gam
    move.list[nrow(move.list)+1,] <<- c(current_move, 1)
    postprob = post.new
  }
  return(list(par = par.out, acc = acc.prob, postprob = postprob))
}

MH.step.fixed <- function(target.dist, proposal.dist, par, data,...){
  # initialize
  par.new = par.old = par
  postprob = 0
  
  # ensure the old model is in the dictionary
  old.model.string = paste0(par$gam, collapse="")
  if (!(old.model.string %in% model.dict[,1])){
    post.old = target.dist(par.old,data=data)
    model.dict[nrow(model.dict)+1,] <<- c(old.model.string,post.old)
  }else{
    idx = which(model.dict[,1] == old.model.string)
    post.old = as.numeric(model.dict[idx,2])
  }
  
  # propose new gamma
  step.new = proposal.dist(par$gam)
  gam.new = step.new$gam.proposal
  log.q = step.new$lograt
  par.new$gam = gam.new
  
  new.model.string = paste0(gam.new, collapse="")
  if (!(new.model.string %in% model.dict[,1])){
    post.new = target.dist(par.new,data=data)
    model.dict[nrow(model.dict)+1,] <<- c(new.model.string,post.new)
  }else{
    idx = which(model.dict[,1] == new.model.string)
    post.new = as.numeric(model.dict[idx,2])
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