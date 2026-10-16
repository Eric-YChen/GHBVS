lBeta = function(...){
  a = c(...)
  return(sum(lgamma(a)) - lgamma(sum(a)))
}

.logsumexp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

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

model.freq = function(chain){
  Models.Freq = apply(do.call(rbind, chain), 1, function(x){as.factor(paste0(x, collapse = ""))})
  return(sort(table(Models.Freq)/length(Models.Freq),decreasing = T))
}

model.postprob = function(chain, postprob){
  chain = do.call(rbind, chain)
  postprob = c(postprob)
  all.models = unique(data.frame(Vars=apply(chain,1,function(x){paste0(x,collapse="")}),prob=postprob))
  all.models$prob = exp(all.models$prob - .logsumexp(all.models$prob))
  Models.Prob = all.models[order(all.models$prob,decreasing = T),]
  Models.Prob$haz <- sapply(Models.Prob$Vars,gam2haz)
  return(Models.Prob)
}

# keys: model strings; logpost: log un-normalised posteriors; visits: post-burn-in visits
model.postprob.dict = function(keys, logpost, visits = rep(0L, length(keys))){
  ok = is.finite(logpost)
  out = data.frame(Vars = keys[ok], logpost = logpost[ok],
                   prob = exp(logpost[ok] - .logsumexp(logpost[ok])),
                   visits = visits[ok], stringsAsFactors = FALSE)
  out = out[order(out$prob, decreasing = TRUE), ]
  out$haz = vapply(out$Vars, gam2haz, character(1))
  rownames(out) = NULL
  return(out)
}

model.postprob2 = function(fullchain){
  chain = fullchain$gam
  postprob = fullchain$post.prob
  all.models = unique(data.frame(Vars=chain,prob=postprob))
  all.models$prob = exp(all.models$prob - .logsumexp(all.models$prob))
  Models.Prob = all.models[order(all.models$prob,decreasing = T),]
  return(Models.Prob)
}

# post.prob: output of model.postprob.dict; p: number of covariates
PIP.postprob = function(post.prob, p){
  G = do.call(rbind, lapply(strsplit(as.character(post.prob$Vars), ""), as.integer))
  if (is.null(dim(G))) G = matrix(G, ncol = p)
  w = post.prob$prob
  PIP = rbind(Time = colSums(w * (G == 1 | G == 3 | G == 4)),
              PH   = colSums(w * (G == 2 | G == 3 | G == 4)))
  PIP = as.data.frame(PIP)
  colnames(PIP) = paste0("x", 1:p)
  return(PIP)
}

# post.prob: output of model.postprob.dict
haz.postprob = function(post.prob){
  hz = tapply(post.prob$prob, factor(post.prob$haz, levels = c("Null","AH","PH","AFT","GH")), sum)
  hz[is.na(hz)] = 0
  return(hz)
}

haz.freq = function(haz.chain){
  Models.Freq = c(haz.chain)
  sort(table(Models.Freq)/length(Models.Freq),decreasing = T)
}

par.gamma <- function(par){
  gam = par$gam
  theta.keep = as.logical(gam == 1 | gam == 3 | gam == 4)
  eta.keep = as.logical(gam == 2 | gam == 3)
  par$theta = par$theta[theta.keep]
  par$eta = par$eta[eta.keep]
  return(par)
}

data.gamma <- function(par, data){
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

gam2bools <- function(gam){
  gam.bools = matrix(0, nrow=2, ncol=length(gam))
  gam.bools[1, gam == 1 | gam == 3 | gam == 4] = 1
  gam.bools[2, gam == 2 | gam == 3 | gam == 4] = 1
  return(gam.bools)
}

bools2gam <- function(gam.bools){
  gam = rep(0, ncol(gam.bools))
  gam[gam.bools[1,] == 1 & gam.bools[2,] == 0] = 1
  gam[gam.bools[1,] == 0 & gam.bools[2,] == 1] = 2
  gam[gam.bools[1,] == 1 & gam.bools[2,] == 1] = 3
  return(gam)
}

par.transf = function(par){
  x = par
  nu = x$nu
  sigma = exp(-nu)
  theta0 = x$theta0*sigma
  theta = -x$theta*sigma
  eta = -x$eta
  return(list(sigma=sigma,mu=theta0,alpha=theta,beta=eta))
}

robust.g <- function(n, p, n_samples=1) {
  g_samples <- ((1 + n)/(p + 1)) / ((1-pmin(runif(n_samples), 1 - 1e-9))^2) - 1
  return(g_samples)
}

metrics <- function(gam, tru.gam){
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

simdat <- function(n,p,mu = rep(0,p),sigma = diag(p), theta = "default", noise.sd = 0,
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

  if (typeof(theta) == "character"){
    if (theta == "default"){
      if (baseline == "PGW"){
        theta = c(1, 1, 2)
      }else if (baseline == "LN"){
        theta = c(1.55,0.7)
      }
    }
  }
  if (cens == "0.25default"){
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
  }else if (cens == "0.5default"){
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

    seed.t <- seeds[i] + 1000000L
    if (hstr == "GH"){
      dat <- simGH(seed = seed.t, n = n, des_h = X, des_t = Xt,
                   theta = theta, beta_h = beta_h, beta_t = beta_t,
                   hstr = hstr, baseline = baseline)
    }else{
      dat <- simGH(seed = seed.t, n = n, des = X, theta = theta,
                   beta = beta, hstr = hstr, baseline = baseline)
    }
    dat <- dat + rnorm(n, 0, noise.sd)
    dat <- pmax(dat, 1e-6)
    if (cens <= 1){
      cens = quantile(dat, max(min(cens + runif(1, -0.05, 0.05), 1), 0.01))
    }
    times <- ifelse(dat < cens, dat, cens)
    status = ifelse(dat < cens, 1, 0)
    cens.rat[i] = mean(status)
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
  n.ev = sum(as.numeric(data$status) == 1)
  max.vars = max(1, floor(n.ev / 5))
  ok.lambda = coxlasso.fit$lambda[coxlasso.fit$df <= max.vars]
  feasible = function(s){
    if (sum(coef(coxlasso.fit, s = s) != 0) <= max.vars) s else min(ok.lambda[ok.lambda >= s])
  }
  s.1se = feasible(cv.coxlasso.fit$lambda.1se); s.min = feasible(cv.coxlasso.fit$lambda.min)
  init.phvars <- coef(coxlasso.fit, s = s.1se)!=0
  if (sum(init.phvars)==0){
    init.phvars <- coef(coxlasso.fit, s = s.min)!=0
    if (sum(init.phvars)==0){
      init.phvars <- c(1,rep(0,ncol(data$x)-1))
      data$x = data$x[,1, drop = F]
      x.coefs = 1
    }else{
      x.keep = which(init.phvars)
      data$x = data$x[,x.keep, drop = F]
      x.coefs = coef(coxlasso.fit, s = s.min)[x.keep]
    }
  }else{
    x.keep = which(init.phvars)
    data$x = data$x[,x.keep, drop = F]
    x.coefs = coef(coxlasso.fit, s = s.1se)[x.keep]
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
  xt.incl <- rep(0, ncol(data$xt))
  if (ncol(data$x) == 0){
    current.crit <- Inf
  }else{
    null.model <- GHMLE(init = c(rep(0, 2), x.coefs), times = data$y,
                        status = data$status, hstr = "PH", dist = "LogNormal",
                        des = data$x, method = "nlminb", maxit = 10000)
    current.crit <- crit(-null.model$OPT$objective, nrow(data$x), length(null.model$OPT$par))
  }
  improved <- TRUE

  history <- data.frame(
    iteration = 0,
    added_var = 0,
    criterion = current.crit
  )
  iter_count <- 0

  while(improved) {
    iter_count <- iter_count + 1

    step.result <- onestep.vs(xt.incl, data, x.coefs, crit)

    if ((step.result$crit < current.crit) & (!is.nan((step.result$crit-current.crit)/current.crit)) & (abs((step.result$crit-current.crit)/current.crit) > 1e-3)){
      added_var <- which(step.result$xt.incl != xt.incl)[1]

      history <- rbind(history, data.frame(
        iteration = iter_count,
        added_var = added_var,
        criterion = step.result$crit
      ))

      current.crit <- step.result$crit
      xt.incl <- step.result$xt.incl
    } else {
      improved <- FALSE
    }
  }

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
  par$nu <- 0
  par$theta0 <- 0
  par$theta <- rep(0, length(par$theta))
  par$eta <- rep(0, length(par$eta))
  type <- gam2haz(par$gam)

  if (type == "GH" | type == "AH" | type == "PH"){
    hess <- hessLLS(par, data, dist = "lnorm")
  }else if (type == "AFT"){
    hess <- hessLLS(par, data, dist = "lnorm")
  }
  return(-hess)
}
