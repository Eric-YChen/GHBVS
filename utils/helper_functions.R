################################################################################
# Helper functions
################################################################################

lBeta = function(...){  # log of Multivariate Beta function
  a = c(...)
  return(sum(lgamma(a)) - lgamma(sum(a)))
}

# Helper function to parse parameter vector and update pars object
# Reduces code duplication across EB_gprior_fixed.R, EB_gprior_robust.R, gprior_setup_fixed.R
parse_par_vec <- function(par.vec, pars, data2) {
  pxt <- ncol(data2$xt)
  px <- ncol(data2$x)
  theta.keep <- as.logical(pars$gam == 1 | pars$gam == 3 | pars$gam == 4)
  eta.keep <- as.logical(pars$gam == 2 | pars$gam == 3)
  ptheta <- sum(theta.keep)
  peta <- sum(eta.keep)

  pars$nu <- par.vec[1]
  pars$theta0 <- par.vec[2]

  if (ptheta == 0) {
    pars$theta <- rep(0, pxt)
    etas <- eta.keep
    etas[etas == 1] <- par.vec[(3 + ptheta):(2 + ptheta + peta)]
    pars$eta <- etas
  } else if (peta == 0) {
    pars$eta <- rep(0, px)
    thetas <- theta.keep
    thetas[thetas == 1] <- par.vec[3:(2 + ptheta)]
    pars$theta <- thetas
  } else {
    thetas <- theta.keep
    thetas[thetas == 1] <- par.vec[3:(2 + ptheta)]
    etas <- eta.keep
    etas[etas == 1] <- par.vec[(3 + ptheta):(2 + ptheta + peta)]
    pars$theta <- thetas
    pars$eta <- etas
  }

  list(
    pars = pars,
    theta.keep = theta.keep,
    eta.keep = eta.keep,
    ptheta = ptheta,
    peta = peta
  )
}

# Posterior Inclusion Probabilities
PIP = function(chains){
  chain.long = do.call(rbind,chains)
  p = ncol(chain.long)
  PIP = matrix(0, nrow = 2, ncol = p)
  for (i in 1:p){
    PIP[1,i] = sum(chain.long[,i] ==1 | chain.long[,i] ==3  | chain.long[,i] ==4)/nrow(chain.long)
    PIP[2,i] = sum(chain.long[,i] ==2 | chain.long[,i] ==3  | chain.long[,i] ==4)/nrow(chain.long)
  }
  PIP = as.data.frame(PIP)
  rownames(PIP) = c("Time","PH")
  colnames(PIP) = paste0("x",1:p)
  return(PIP)
}

# Model frequencies
model.freq = function(chain){
  Models.Freq = apply(do.call(rbind, chain), 1, function(x){as.factor(paste0(x, collapse = ""))})
  return(sort(table(Models.Freq)/length(Models.Freq),decreasing = T))
}

# Model postprob
model.postprob = function(chain, postprob){
  chain = do.call(rbind, chain)
  postprob = c(postprob)
  all.models = unique(data.frame(Vars=apply(chain,1,function(x){paste0(x,collapse="")}),prob=postprob))
  all.models$prob = exp(all.models$prob - mclust::logsumexp(all.models$prob))
  Models.Prob = all.models[order(all.models$prob,decreasing = T),]
  Models.Prob$haz <- sapply(Models.Prob$Vars,gam2haz)
  return(Models.Prob)
}


model.postprob2 = function(fullchain){
  chain = fullchain$gam
  postprob = fullchain$post.prob
  all.models = unique(data.frame(Vars=chain,prob=postprob))
  all.models$prob = exp(all.models$prob - mclust::logsumexp(all.models$prob))
  Models.Prob = all.models[order(all.models$prob,decreasing = T),]
  return(Models.Prob)
}

# Hazard structure frequencies
haz.freq = function(haz.chain){
  Models.Freq = c(haz.chain)
  sort(table(Models.Freq)/length(Models.Freq),decreasing = T)
}

# subsetting parameters according to gamma
par.gamma <- function(par){
  # 0: (0,0), 1: (1,0), 2: (0,1), 3: (1,1)
  # (theta,eta) / (alpha,beta) / (time,PH) / (xt,x)
  gam = par$gam
  theta.keep = as.logical(gam == 1 | gam == 3 | gam == 4)
  eta.keep = as.logical(gam == 2 | gam == 3)
  par$theta = par$theta[theta.keep]
  par$eta = par$eta[eta.keep]
  return(par)
}

# subsetting data according to gamma
data.gamma <- function(par, data){
  # 0: (0,0), 1: (1,0), 2: (0,1), 3: (1,1)
  # (theta,eta) / (alpha,beta) / (time,PH) / (xt,x)
  
  # prevent accidental reusing of function
  if (ncol(data$xt) < length(par$gam) | ncol(data$x) < length(par$gam)){
    return(data)
  }
  
  gam = par$gam
  xt.keep = as.logical(gam == 1 | gam == 3 | gam == 4)
  x.keep = as.logical(gam == 2 | gam == 3)
  if (sum(xt.keep) == 0){
    data$xt = numeric(0)
    data$xt2 = numeric(0)
  }else{
    data$xt = data$xt[,xt.keep,drop=FALSE]
    data$xt2 = t(data$xt)%*%data$xt
    data$xt2.inv = solve(data$xt2)
  }
  if (sum(x.keep) == 0){
    data$x = numeric(0)
    data$x2 = numeric(0)
  }else{
    data$x = data$x[,x.keep,drop=FALSE]
    data$x2 = t(data$x)%*%data$x
    data$x2.inv = solve(data$x2)
  }
  return(data)
}

# vector of gamma (0,1,2,3) to matrix of boolean
gam2bools <- function(gam){
  # 0: (0,0), 1: (1,0), 2: (0,1), 3: (1,1)
  # (theta,eta) / (alpha,beta) / (time,PH) / (xt,x)
  gam.bools = matrix(0, nrow=2, ncol=length(gam))
  gam.bools[1, gam == 1 | gam == 3 | gam == 4] = 1
  gam.bools[2, gam == 2 | gam == 3 | gam == 4] = 1
  return(gam.bools)
}

bools2gam <- function(gam.bools){
  # 0: (0,0), 1: (1,0), 2: (0,1), 3: (1,1)
  # (theta,eta) / (alpha,beta) / (time,PH) / (xt,x)
  gam = rep(0, ncol(gam.bools))
  gam[gam.bools[1,] == 1 & gam.bools[2,] == 0] = 1
  gam[gam.bools[1,] == 0 & gam.bools[2,] == 1] = 2
  gam[gam.bools[1,] == 1 & gam.bools[2,] == 1] = 3
  return(gam)
}

# transform parameters to original scale (sigma, mu, alpha, beta)
par.transf = function(par){
  x = par
  nu = x$nu
  sigma = exp(-nu)
  theta0 = x$theta0*sigma
  theta = -x$theta*sigma
  eta = -x$eta
  return(list(sigma=sigma,mu=theta0,alpha=theta,beta=eta))
}

# enumerate all possible models
modelprob.grid <- function(p,gtheta=1,geta=1){
  grid = expand.grid(rep(list(0:3), p))
  grid = grid[-(grid[,1]==0 & grid[,2]==0),]
  record = matrix(0,nrow=nrow(grid),ncol=2)
  itervec<<-0
  #colnames(record) = c("theta1","theta2","eta1","eta2","g","postprob")
  for (i in 1:nrow(grid)){
    modelpar = list(inits = list(nu=rnorm(1), theta0=rnorm(1),
                                 theta=rnorm(p,0,0.5), eta=rnorm(p,0,0.5),
                                 g.theta = gtheta, g.eta = geta,
                                 gam = grid[i,],
                                 nu.alpha = 0.01, nu.beta = 0.01),
                    data = list(y=y,x=x,xt=xt,status=status),
                    target.dist = targetdist)
    par1 = modelpar$inits
    theta.keep = as.logical(par1$gam == 1 | par1$gam == 3)
    eta.keep = as.logical(par1$gam == 2 | par1$gam == 3)
    MAP.pars <<- numeric(0)
    
    record[i,1] = paste0(grid[i,],collapse = "")
    record[i,2] = as.numeric(targetdist(par = modelpar$inits,data = modelpar$data))
    #rm(MAP.pars)
  }
  # rank by posterior probability
  #record[,3] = itervec
  #record$rank = rank(-record[,2])
  return(record)
}

robust.g <- function(n, p, n_samples=1) {
  # Apply inverse CDF:
  # g = (g_min + b)/((1-u)^2) - b
  g_samples <- ((1 + n)/(p + 1)) / ((1-pmin(runif(n_samples), 1 - 1e-9))^2) - 1
  return(g_samples)
}


metrics <- function(gam, tru.gam){  # consider time & hazard part separately!!
    if (is.character(gam)){
      gam = as.numeric(strsplit(gam,split="")[[1]])
    }
    if (is.character(tru.gam)){
      tru.gam = as.numeric(strsplit(tru.gam,split="")[[1]])
    }
    TP = sum(gam !=0 & tru.gam !=0)
    TN = sum(gam == 0 & tru.gam == 0)
    FP = sum(gam != 0 & tru.gam == 0)
    FN = sum(gam == 0 & tru.gam != 0)
    Spec = TN/(FP + TN)
    Sens = TP/(TP + FN)
    FDR = FP/(TP + FP)
    n = length(gam)
    phi = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
    return(data.frame(Sens = Sens, Spec = Spec, FDR = FDR, phi = phi))
}

gam2haz = function(gam){
  if (is.character(gam)){
    gam = as.numeric(strsplit(gam,split="")[[1]])
  }
  if (all(gam==0)){
    return("Null")
  }else if (all(gam %in% c(0,1))){
    return("AH")
  }else if (all(gam %in% c(0,2))){
    return("PH")
  }else if (all(gam!=4)){
    return("GH")
  }else if (all(gam %in% c(0,4))){
    return("AFT")
  }else{
    return("Invalid model")
  }
}

simdat <- function(n,p,mu = rep(0,p),sigma = diag(p), theta = "default",
                   beta = c(1, -1, 0.25, -0.25, rep(0,p-4)), 
                   beta_t = c(1, 0.25, -1, -0.25, rep(0,p-4)), 
                   beta_h = c(0.25, -1, -0.25, 1, rep(0,p-4)),
                   hstr = "GH", baseline = "LN", cens = "0.25default",
                   n.simset=1, seeds = c(1:n.simset), out.format="list"){
  beta.active <- beta[1:4]
  beta <- rep(0,p)
  beta[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- beta.active
  
  beta_t.active <- beta_t[1:4]
  beta_t <- rep(0,p)
  beta_t[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- beta_t.active
  
  beta_h.active <- beta_h[1:4]
  beta_h <- rep(0,p)
  beta_h[as.integer(quantile(c(1:p),c(0.2,0.4,0.6,0.8)))] <- beta_h.active
  
  # Check inputs
  if (typeof(theta) == "character"){
    if (theta == "default"){
      if (baseline == "PGW"){
        theta = c(1, 1, 2)
      }else if (baseline == "LN"){
        theta = c(1.55,0.7)
      }
    }
  }
  if (cens == "0.25default"){  # around 25% censoring
    if (hstr == "GH" & baseline == "PGW"){
      cens = 11
    }else if (hstr == "GH" & baseline == "LN"){
      cens = 17
    }else if (hstr == "PH" & baseline == "PGW"){
      cens = 8
    }else if (hstr == "PH" & baseline == "LN"){
      cens = 10
    }else if (hstr == "AH" & baseline == "PGW"){
      cens = 6
    }else if (hstr == "AH" & baseline == "LN"){
      cens = 10
    }else if (hstr == "AFT" & baseline == "PGW"){
      cens = 6
    }else if (hstr == "AFT" & baseline == "LN"){
      cens = 12
    }
  }else if (cens == "0.5default"){  # around 50% censoring
    if (hstr == "GH" & baseline == "PGW"){
      cens = 1.8
    }else if (hstr == "GH" & baseline == "LN"){
      cens = 6.5
    }else if (hstr == "PH" & baseline == "PGW"){
      cens = 1.8
    }else if (hstr == "PH" & baseline == "LN"){
      cens = 4.5
    }else if (hstr == "AH" & baseline == "PGW"){
      cens = 2.2
    }else if (hstr == "AH" & baseline == "LN"){
      cens = 6.4
    }else if (hstr == "AFT" & baseline == "PGW"){
      cens = 1.8
    }else if (hstr == "AFT" & baseline == "LN"){
      cens = 4.8
    }
  }
  
  out = list()
  cens.rat = rep(0, n.simset)
  for (i in 1:n.simset){
    set.seed(seeds[i])
    x <- rmvnorm(n = n, mean = mu, sigma = sigma)
    X = x
    Xt = x
    
    if (hstr == "GH"){
      dat <- simGH(seed = seeds[i], n = n, des_h = X, des_t = Xt,
                   theta = theta, beta_h = beta_h, beta_t = beta_t,
                   hstr = hstr, baseline = baseline)
    }else{
      dat <- simGH(seed = seeds[i], n = n, des = X, theta = theta,
                   beta = beta, hstr = hstr, baseline = baseline)
    }

    if (cens <= 1){
      cens = quantile(dat, max(min(cens + runif(1, -0.05, 0.05), 1), 0.01))
    }
    times <- ifelse(dat < cens, dat, cens)
    status = ifelse(dat < cens, 1, 0)
    cens.rat[i] = mean(status)
    # Create data frame
    if (out.format == "list"){
      out[[i]] <- list(x=X,xt=Xt,y=times,status=status)
    }else{
      out[[i]] <- data.frame(x,y=times,status=status)
    }
  }
  return(list(data = out, cens.rat = cens.rat))
}

init.gam <- function(data, crit = "BIC", history.keep = F, simplify = F){
  cv.coxlasso.fit = cv.glmnet(x = data$x, y = Surv(data$y, data$status), nfolds=10, 
                              family = "cox", maxit=10000, alpha = 1)
  coxlasso.fit = glmnet(x = data$x, y = Surv(data$y, data$status), nfolds=10, 
                        family = "cox", maxit=10000, alpha = 1)
  init.phvars <- coef(coxlasso.fit, s = cv.coxlasso.fit$lambda.1se)!=0
  if (sum(init.phvars)==0){
    init.phvars <- coef(coxlasso.fit, s = cv.coxlasso.fit$lambda.min)!=0
    if (sum(init.phvars)==0){
      init.phvars <- c(1,rep(0,ncol(data$x)-1))
      data$x = data$x[,1, drop = F]
      x.coefs = 1
    }else{
      x.keep = which(init.phvars)
      data$x = data$x[,x.keep, drop = F]
      x.coefs = coef(coxlasso.fit, s = cv.coxlasso.fit$lambda.min)[x.keep]
    }
  }else{
    x.keep = which(init.phvars)
    data$x = data$x[,x.keep, drop = F]
    x.coefs = coef(coxlasso.fit, s = cv.coxlasso.fit$lambda.1se)[x.keep]
  }
  
  if (!simplify){
    init.vals = step.GH(data, x.coefs = x.coefs, crit = crit)
    history = init.vals$history
    init.vals = bools2gam(rbind(init.vals$xt.incl, as.vector(init.phvars)*1))
    if (history.keep){
      return(list(
        xt.incl = init.vals,
        history = history
      ))
    }else{
      return(init.vals)
    }
  }else{
    init.vals = as.vector(3*init.phvars + rep(0,ncol(data$x)))
    return(init.vals)
  }
  
}

step.GH <- function(data, x.coefs = NULL, crit = "BIC"){
  if (crit == "BIC"){
    crit = function(objective, n, p){
      return(-2 * objective + p * log(n))
    }
  }else if (crit == "AIC"){
    crit = function(objective, n, p){
      return(-2 * objective + 2 * p)
    }
  }
  if (is.null(x.coefs)){
    x.coefs = rep(0, ncol(data$x))
  }
  # Initialize variables
  xt.incl <- rep(0, ncol(data$xt))  # Start with no variables included
  if (ncol(data$x) == 0){
    current.crit <- Inf
  }else{
    null.model <- GHMLE(init = c(rep(0, 2), x.coefs), times = data$y,
                        status = data$status, hstr = "PH", dist = "LogNormal",
                        des = data$x, method = "nlminb", maxit = 10000)
    current.crit <- crit(-null.model$OPT$objective, nrow(data$x), length(null.model$OPT$par))
  }
  improved <- TRUE
  
  # Initialize tracking variables
  history <- data.frame(
    iteration = 0,
    added_var = 0,
    criterion = current.crit
  )
  iter_count <- 0
  
  # Iterate until no improvement is found
  while(improved) {
    iter_count <- iter_count + 1
    
    # Perform one step of forward selection
    step.result <- onestep.vs(xt.incl, data, x.coefs, crit)
    
    # Check if we got an improvement
    if ((step.result$crit < current.crit) & (!is.nan((step.result$crit-current.crit)/current.crit)) & (abs((step.result$crit-current.crit)/current.crit) > 1e-3)){
      added_var <- which(step.result$xt.incl != xt.incl)[1]
      
      # Update tracking
      history <- rbind(history, data.frame(
        iteration = iter_count,
        added_var = added_var,
        criterion = step.result$crit
      ))
      
      # Update current state
      current.crit <- step.result$crit
      xt.incl <- step.result$xt.incl
    } else {
      improved <- FALSE
    }
  }
  
  # Return the final model and history
  return(list(
    xt.incl = xt.incl,
    criterion = current.crit,
    history = history
  ))
}

onestep.vs <- function(xt.incl, data, x.coefs, crit){
  n = nrow(data$x)
  iter.vars = which(xt.incl == 0)
  best.crit = Inf
  best.xt.incl = xt.incl
  for (i in iter.vars){
    xt.incl.new = xt.incl
    xt.incl.new[i] = 1
    xt = data$xt[,as.logical(xt.incl.new),drop = F]
    if (ncol(data$x) != 0){
      model <- GHMLE(init = c(rep(0, 2),x.coefs, rep(0,ncol(xt))), times = data$y, 
                     status = data$status, hstr = "GH", dist = "LogNormal", 
                     des = data$x, des_t = xt, method = "nlminb", maxit = 5000)
    }else{
      model <- GHMLE(init = c(rep(0, 2),rep(0,ncol(xt))), times = data$y, 
                     status = data$status, hstr = "AH", dist = "LogNormal", 
                     des_t = xt, method = "nlminb", maxit = 5000)
    }
   
    crit.model <- crit(-model$OPT$objective, n, length(model$OPT$par))
    if (crit.model < best.crit){
      best.crit <- crit.model
      best.xt.incl <- xt.incl.new
    }
  }
  return(list(xt.incl = best.xt.incl, crit = best.crit))
}

g.robtrans <- function(g = NULL, g.trans = NULL, par, h.type = "t"){
  n = par$n
  if (h.type == "t" | h.type == "af"){
    p = sum(par$gam == 1 | par$gam == 3 | par$gam == 4)
  }else if (h.type == "p" | h.type == "h"  | h.type == "ph"){
    p = sum(par$gam == 2 | par$gam == 3)
  }
  if (is.null(g) & is.null(g.trans)){
    print("Provide one of g / g.trans.")
    return(NULL)
  }
  if (is.null(g.trans)){
    g.trans = log(g - (n+1)/(p+1) + 1)
    return(g.trans)
  }else if (is.null(g)){
    g = exp(g.trans) + (n+1)/(p+1) - 1
    return(g)
  }
}

compute_fisher_info <- function(par, data) {
  # Create zero parameter vector for Fisher info evaluation
  par$nu <- 0
  par$theta0 <- 0
  par$theta <- rep(0, length(par$theta))
  par$eta <- rep(0, length(par$eta))
  type <- gam2haz(par$gam)
  
  # Get Hessian at zero (this gives us the negative Fisher information)
  if (type == "GH" | type == "AH" | type == "PH"){
    #par$gam <- rep(3, length(par$gam))
    hess <- hessLLS(par, data, dist = "lnorm")
  }else if (type == "AFT"){
    #par$gam <- rep(4, length(par$gam))
    hess <- hessLLS(par, data, dist = "lnorm")
  }
  return(-hess)  # Return the negative Hessian as Fisher information
}
