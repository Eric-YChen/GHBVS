# Set BVS_CORE to the path of this folder before sourcing (default: working directory).
if (!exists("BVS_CORE")) BVS_CORE <- "."
source(file.path(BVS_CORE, "utils", "routines2.R"))
source(file.path(BVS_CORE, "utils", "helper_functions.R"))
source(file.path(BVS_CORE, "HazReg", "routines.R"))

draw1 <- function(x) x[sample.int(length(x), 1)]

get_move_probs <- function(gam) {
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
  } else {
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

  if ("A/D_GH" %in% names(probs) && length(get_valid_idx(gam)) == 0) {
    probs <- probs[names(probs) != "A/D_GH"]
    probs <- probs / sum(probs)
  }

  return(probs)
}

get_move_prob <- function(gam, move_type) {
  probs <- get_move_probs(gam)
  if (move_type %in% names(probs)) return(probs[move_type])
  return(0)
}

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

  move.probs <- get_move_probs(gam.old)
  if (length(move.probs) == 1){
    move <- names(move.probs)
  }else{
    move <- sample(names(move.probs), 1, prob = move.probs)
  }

  current_move <<- move

  if (move == "A/D_GH") {
    valid_idx <- get_valid_idx(gam.old)
    n_valid <- length(valid_idx)
    idx <- draw1(valid_idx)

    if (gam.candid[idx] > 0) {
      gam.candid[idx] <- 0
      current_move <<- "D_GH"
      new_haz <- gam2haz(gam.candid)

      if (new_haz == "Null") {
        q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid
        q_rev <- 1.0 / p / 4
      } else {
        n_valid_new <- length(get_valid_idx(gam.candid))
        q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid
        q_rev <- ifelse(n_valid_new == 0, 0,
                        get_move_prob(gam.candid, "A/D_GH") / n_valid_new / 3)
      }
      log.q <- log(q_rev) - log(q_fwd)
    } else {
      new_val <- sample.int(3, 1)
      gam.candid[idx] <- new_val
      current_move <<- "A_GH"

      n_valid_new <- length(get_valid_idx(gam.candid))
      q_fwd <- get_move_prob(gam.old, "A/D_GH") / n_valid / 3
      q_rev <- ifelse(n_valid_new == 0, 0,
                      get_move_prob(gam.candid, "A/D_GH") / n_valid_new)
      log.q <- log(q_rev) - log(q_fwd)
    }
  } else if (move == "A/D") {
    idx <- sample.int(p, 1)

    if (gam.candid[idx] > 0) {
      gam.candid[idx] <- 0
      current_move <<- "D"
      new_haz <- gam2haz(gam.candid)

      if (new_haz == "Null") {
        q_fwd <- get_move_prob(gam.old, "A/D") / p
        q_rev <- 1.0 / p / 4
      } else {
        q_fwd <- get_move_prob(gam.old, "A/D") / p
        q_rev <- get_move_prob(gam.candid, "A/D") / p
      }
      log.q <- log(q_rev) - log(q_fwd)
    } else {
      gam.candid[idx] <- c(1,2,4)[haz == c("AH","PH","AFT")]
      current_move <<- "A"

      q_fwd <- get_move_prob(gam.old, "A/D") / p
      q_rev <- get_move_prob(gam.candid, "A/D") / p
      log.q <- log(q_rev) - log(q_fwd)
    }
  } else if (move == "C") {
    nonzero_idx <- which(gam.candid != 0)
    n_nonzero <- length(nonzero_idx)
    idx <- draw1(nonzero_idx)
    gam.candid[idx] <- draw1(setdiff(1:3, gam.candid[idx]))

    q_fwd <- get_move_prob(gam.old, "C") / n_nonzero / 2
    q_rev <- get_move_prob(gam.candid, "C") / n_nonzero / 2
    log.q <- log(q_rev) - log(q_fwd)
  } else if (move == "S") {
    idx1 <- sample.int(p, 1)
    gam.1 <- gam.candid[idx1]
    other_idx <- which(gam.candid != gam.1)
    n_other <- length(other_idx)
    idx2 <- draw1(other_idx)
    gam.2 <- gam.candid[idx2]

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

    if (gam.1 + gam.2 == 3) {
      q_fwd <- q_fwd / 2
    }
    if (gam.candid[idx1] + gam.candid[idx2] == 3) {
      q_rev <- q_rev / 2
    }

    log.q <- log(q_rev) - log(q_fwd)
  } else if (move == "C*") {
    if (haz == "AFT") {
      new_val <- sample.int(3, 1)
      gam.candid[gam.candid == 4] <- new_val
      new_haz <- gam2haz(gam.candid)

      q_fwd <- get_move_prob(gam.old, "C*") / 3
      if (new_haz == "GH") {
        q_rev <- get_move_prob(gam.candid, "C*")
        current_move <<- "AFT_GH"
      } else {
        q_rev <- get_move_prob(gam.candid, "C*") / 2
        current_move <<- ifelse(new_val == 1, "AFT_AH", "AFT_PH")
      }
      log.q <- log(q_rev) - log(q_fwd)
    } else if (haz == "GH") {
      gam.candid[gam.candid == 3] <- 4
      current_move <<- "GH_AFT"

      q_fwd <- get_move_prob(gam.old, "C*")
      q_rev <- get_move_prob(gam.candid, "C*") / 3
      log.q <- log(q_rev) - log(q_fwd)
    } else {
      nonzero_idx <- which(gam.candid != 0)
      current_val <- gam.candid[nonzero_idx[1]]
      new_val <- draw1(setdiff(c(1,2,4), current_val))
      gam.candid[nonzero_idx] <- new_val
      new_haz <- gam2haz(gam.candid)

      q_fwd <- get_move_prob(gam.old, "C*") / 2
      if (new_haz == "AFT") {
        q_rev <- get_move_prob(gam.candid, "C*") / 3
        current_move <<- ifelse(haz == "PH", "PH_AFT", "AH_AFT")
      } else {
        q_rev <- get_move_prob(gam.candid, "C*") / 2
        current_move <<- ifelse(haz == "PH", "PH_AH", "AH_PH")
      }
      log.q <- log(q_rev) - log(q_fwd)
    }
  } else if (move == "A_null") {
    idx <- sample.int(p, 1)
    new_val <- sample.int(4, 1)
    gam.candid[idx] <- new_val
    new_haz <- gam2haz(gam.candid)

    q_fwd <- 1.0 / p / 4
    if (new_haz == "GH") {
      n_valid_new <- length(get_valid_idx(gam.candid))
      q_rev <- ifelse(n_valid_new == 0, 0,
                      get_move_prob(gam.candid, "A/D_GH") / n_valid_new)
    } else {
      q_rev <- get_move_prob(gam.candid, "A/D") / p
    }
    log.q <- log(q_rev) - log(q_fwd)
  }

  return(list(gam.proposal = gam.candid, lograt = log.q))
}

dict.get <- function(key) get0(key, envir = model.env, inherits = FALSE)
dict.set <- function(key, val) assign(key, val, envir = model.env)

MH.step.fixed <- function(target.dist, proposal.dist, par, data, ...){
  par.new = par.old = par

  old.key = paste0(par$gam, collapse="")
  post.old = dict.get(old.key)
  if (is.null(post.old)){
    post.old = target.dist(par.old, data=data)
    dict.set(old.key, post.old)
  }

  step.new = proposal.dist(par$gam)
  par.new$gam = step.new$gam.proposal
  log.q = step.new$lograt

  new.key = paste0(par.new$gam, collapse="")
  post.new = dict.get(new.key)
  if (is.null(post.new)){
    post.new = target.dist(par.new, data=data)
    dict.set(new.key, post.new)
  }

  log.acc.prob = min(0, post.new - post.old + log.q)
  acc.prob = exp(log.acc.prob)
  if (is.nan(acc.prob)) acc.prob = 0

  accepted = runif(1) <= acc.prob
  par.out = par
  par.out$gam = if (accepted) par.new$gam else par.old$gam
  return(list(par = par.out, acc = acc.prob, accepted = accepted,
              postprob = if (accepted) post.new else post.old,
              key = if (accepted) new.key else old.key,
              move = current_move))
}

# model.par: list(inits, data, gam.param); mcmc.par: list(n.iter, n.burnin, n.thin, n.chain, printevery)
# gprior.mat: "LCM", "LCM-hyper" or "Product"; gamprior: "default", "custom" (model.par$gam.param) or "unif"
# null.model: admit the null model; target: optional function(par, data, ...) returning the log posterior
BVS.MCMC <- function(model.par, mcmc.par,
                     gprior.mat = "LCM", gamprior = "default",
                     null.model = TRUE,
                     verbose = T, target = NULL){
  inits <- model.par$inits
  data <- model.par$data

  for (nm in c("x", "xt")){
    if (is.null(dim(data[[nm]]))){
      warning("data$", nm, " is a vector; coercing to a one-column matrix. ",
              "Pass x/xt as matrices, e.g. matrix(v, ncol = 1).", call. = FALSE)
      data[[nm]] <- as.matrix(data[[nm]])
    }else if (!is.matrix(data[[nm]])){
      warning("data$", nm, " is a ", class(data[[nm]])[1],
              "; coercing with as.matrix().", call. = FALSE)
      data[[nm]] <- as.matrix(data[[nm]])
    }
    if (!is.numeric(data[[nm]])){
      stop("data$", nm, " must be numeric; got ", typeof(data[[nm]]),
           ". Convert factor/character columns to numeric first.", call. = FALSE)
    }
  }

  admit.null <<- isTRUE(null.model)

  inits$g <- nrow(data$x)*inits$g
  data$xt2 = t(data$xt)%*%data$xt; data$x2 = t(data$x)%*%data$x
  data$xt2.inv = solve(data$xt2); data$x2.inv = solve(data$x2)
  if (is.null(inits$n)){inits$n <- nrow(data$x)}
  inits$n.events <- sum(as.numeric(data$status) == 1)

  n.iter <- mcmc.par$n.iter
  n.burn <- if (!is.null(mcmc.par[["n.burnin"]])) mcmc.par[["n.burnin"]] else mcmc.par[["n.burn"]]
  n.thin <- mcmc.par$n.thin
  if (is.null(n.burn)){stop("mcmc.par$n.burnin must be given")}
  if (n.thin < 1){stop("n.thin must be integer >= 1")}
  if (n.iter < n.burn){stop("n.iter must be greater than n.burnin")}
  if (n.iter%%n.thin != 0){stop("n.iter must be divisible by n.thin")}
  printevery <- mcmc.par$printevery; n.chain <- mcmc.par$n.chain
  n.save = n.iter%/%n.thin - n.burn%/%n.thin

  gamprior <<- gamprior
  gam.param <<- model.par$gam.param
  if (gamprior == "custom" && is.null(gam.param)) stop("gamprior = 'custom' requires model.par$gam.param")

  model.env <<- new.env(hash = TRUE)
  if (gprior.mat == "LCM"){
    source(file.path(BVS_CORE, "utils", "EB_gprior_fixed.R"))
  }else if (gprior.mat == "LCM-hyper"){
    source(file.path(BVS_CORE, "utils", "EB_gprior_robust.R"))
  }else if (gprior.mat == "Product"){
    source(file.path(BVS_CORE, "utils", "gprior_setup_fixed.R"))
  }else{
    stop("gprior.mat must be one of 'LCM', 'LCM-hyper', 'Product'")
  }
  target.dist <- if (is.null(target)) targetdist else target
  gam.proposal <- ADS.proposal.aft
  MH.step <- MH.step.fixed

  if (gam2haz(inits$gam) == "Invalid model" || !is.finite(prior.gam(inits, log = TRUE)) ||
      (all(inits$gam == 0) && !admit.null)){
    stop("the initial model has zero prior probability")
  }

  chain.list = list()
  haz.list = matrix(NA, nrow=n.save, ncol=n.chain)
  modelsize.df = matrix(NA_real_, nrow=n.iter, ncol=n.chain)
  postprob.list = matrix(NA_real_, nrow=n.save, ncol=n.chain)
  acc.list = matrix(NA_real_, nrow=n.iter, ncol=n.chain)
  moves <- character(n.iter*n.chain); accepts <- logical(n.iter*n.chain)
  visit.env <- new.env(hash = TRUE)

  for(chain in 1:n.chain){
    par = inits
    current.chain = matrix(NA_integer_, nrow=n.save, ncol=length(par$gam))
    current.haz = rep(NA, n.save)
    iter.save = 0
    for(iter in 1:n.iter){
      update = MH.step(target.dist, proposal.dist=gam.proposal, par, data)
      par.new = update$par
      k = (chain - 1)*n.iter + iter
      moves[k] = update$move; accepts[k] = update$accepted

      acc.list[iter,chain] = update$acc
      modelsize.df[iter,chain] = sum(par.new$gam!=0)

      if (iter > n.burn){
        visit.env[[update$key]] <- get0(update$key, envir = visit.env, inherits = FALSE, ifnotfound = 0L) + 1L
      }

      if(iter > n.burn & iter %% n.thin == 0){
        iter.save = iter.save + 1
        current.chain[iter.save,] = par.new$gam
        current.haz[iter.save] = gam2haz(par.new$gam)
        postprob.list[iter.save,chain] = update$postprob
      }
      if(verbose & iter %% printevery == 0){
        print(paste0("Chain:", chain, ", Iter:", iter,
        ", Acc:", round(mean(acc.list[(iter-printevery+1):iter,chain]),3),
        ", modl.size: ", round(mean(modelsize.df[(iter-printevery+1):iter,chain]),3) ))
      }
      par = par.new
    }
    chain.list[[chain]] = as.data.frame(current.chain)
    haz.list[,chain] = current.haz
  }

  keys <- ls(model.env)
  lp.all <- vapply(keys, function(k) get(k, envir = model.env), numeric(1))
  model.dict <- data.frame(gam = keys, postprob = lp.all, row.names = NULL, stringsAsFactors = FALSE)
  post.prob <- model.postprob.dict(keys, logpost = lp.all,
                                   visits = vapply(keys, function(k) get0(k, envir = visit.env, inherits = FALSE,
                                                                         ifnotfound = 0L), integer(1)))
  p <- length(inits$gam)

  return(list(chain = chain.list,
              acc = as.data.frame(acc.list), modelsize = as.data.frame(modelsize.df),
              PIP = PIP(chain.list),
              PIP.post = PIP.postprob(post.prob, p),
              hazPIP = haz.freq(haz.list),
              haz.post = haz.postprob(post.prob),
              model.freq = model.freq(chain.list),
              model.prob = model.dict,
              post.prob = post.prob,
              move.list = data.frame(move = moves, accept = accepts),
              opt.log = opt.log))
}

log.gh.model.weight <- function(m, k, q){
  kk <- 0:m
  Nk <- choose(m, kk) * 2^(m - kk) - 2 * (kk == 0)
  ck <- ifelse(kk == 0, 1 - 2^(1 - m), 1)
  ok <- Nk > 0
  lw <- rep(-Inf, m + 1)
  lw[ok] <- -(dbinom(kk[ok], m, q, log = TRUE) + log(ck[ok]))
  lw[kk == k] - .logsumexp(log(Nk[ok]) + lw[ok])
}

# par4: parameter list with gam (and n.events); type: "default", "custom" or "unif"
# parameters (custom): list(layer1.betabin.par = c(a, b) or "uniform",
#   layer2.mulnom.par = c(h_AH, h_PH, h_GH, h_AFT), layer3.beta.par = q)
prior.gam = function(par4, log=T, type = gamprior,
                     parameters = gam.param){
  out <- function(lp) if (log) lp else exp(lp)
  gam = par4$gam; p = length(gam)
  m = sum(gam != 0); p0 = p - m
  curr.model = gam2haz(gam)
  if (curr.model == "Invalid model") return(out(-Inf))

  d.gam = sum(gam %in% c(1, 2, 4)) + 2 * sum(gam == 3)
  if (!is.null(par4$n.events) && d.gam + 2 > par4$n.events) return(out(-Inf))

  if (type == "unif") return(out(0))

  ab = c(1, 1); h = c(AH = 1, PH = 1, GH = 1, AFT = 1); q = 1/3
  if (type == "custom"){
    if (!is.null(parameters$layer1.betabin.par)) ab = parameters$layer1.betabin.par
    if (!is.null(parameters$layer2.mulnom.par))
      h = setNames(parameters$layer2.mulnom.par, c("AH", "PH", "GH", "AFT"))
    if (!is.null(parameters$layer3.beta.par)) q = parameters$layer3.beta.par
  }else if (type != "default"){
    stop("unknown model-space prior type: ", type)
  }

  lp = if (identical(ab, "uniform")) -p * log(2) else lbeta(ab[1] + m, ab[2] + p0) - lbeta(ab[1], ab[2])
  if (m == 0) return(out(lp))

  lp = lp + log(h[[curr.model]]) - log(sum(h))

  if (curr.model == "GH") lp = lp + log.gh.model.weight(m, sum(gam == 3), q)

  out(lp)
}
