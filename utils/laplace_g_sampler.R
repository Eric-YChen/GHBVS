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
