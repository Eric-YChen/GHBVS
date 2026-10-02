# ---------------------------------------------------------------------------- #
#  MLE cache (per-worker, shared across BIC chains within one simulation)      #
# ---------------------------------------------------------------------------- #
mle_cache <- new.env(hash = TRUE, parent = emptyenv())

reset_mle_cache <- function() {
  rm(list = ls(envir = mle_cache), envir = mle_cache)
  invisible(NULL)
}

# Look up or compute the MLE log-likelihood and effective dimension for a
# given gamma. Returns list(log_lik_map, d_gam, ok).
get_or_compute_mle <- function(par, data, max.iter = 100) {
  key <- paste0(par$gam, collapse = "")
  if (exists(key, envir = mle_cache, inherits = FALSE)) {
    return(get(key, envir = mle_cache))
  }
  par_mle <- tryCatch(mle.gh(par, data, max.iter = max.iter),
                      error = function(e) NULL)
  if (is.null(par_mle)) {
    res <- list(log_lik_map = NA_real_, d_gam = NA_integer_, ok = FALSE)
  } else {
    ll <- tryCatch(lglik(par_mle, data), error = function(e) NA_real_)
    theta_keep <- as.logical(par$gam == 1 | par$gam == 3 | par$gam == 4)
    eta_keep   <- as.logical(par$gam == 2 | par$gam == 3)
    d_gam <- sum(theta_keep) + sum(eta_keep)
    res <- list(log_lik_map = ll, d_gam = d_gam, ok = is.finite(ll))
  }
  assign(key, res, envir = mle_cache)
  res
}

# ---------------------------------------------------------------------------- #
#  BIC log-marginal-likelihood approximation                                   #
# ---------------------------------------------------------------------------- #
marg.loglik.bic <- function(par, data,
                            n_type = c("n", "ev"),
                            max.iter = 100) {
  n_type <- match.arg(n_type)
  res <- get_or_compute_mle(par, data, max.iter)
  if (!isTRUE(res$ok)) return(-Inf)
  k_gam <- res$d_gam + 2L                       # +2 for common pars (nu, theta0)
  n_eff <- if (n_type == "n") nrow(data$x)
           else max(sum(data$status == 1L), 1L)
  res$log_lik_map - 0.5 * k_gam * log(n_eff)
}

# ---------------------------------------------------------------------------- #
#  Target-distribution factory                                                 #
# ---------------------------------------------------------------------------- #
#   prior_type:
#     "flat" - uniform model prior (still enforces the n.threshold cap on
#              models that violate p > n_threshold)
#     "bb"   - generalised Beta-Binomial prior used by the proposed methods
#              (prior.gam, type = "default")
make.targetdist.bic <- function(prior_type = c("flat", "bb"),
                                n_type     = c("n",    "ev")) {
  prior_type <- match.arg(prior_type)
  n_type     <- match.arg(n_type)
  gp_type    <- if (prior_type == "flat") "unif" else "default"
  function(par, data, max.iter = 100, ...) {
    if (sum(par$gam) == 0) return(-Inf)
    log_bic <- marg.loglik.bic(par, data,
                               n_type   = n_type,
                               max.iter = max.iter)
    if (!is.finite(log_bic)) return(-Inf)
    prior.gam(par, log = TRUE, type = gp_type) + log_bic
  }
}

# ---------------------------------------------------------------------------- #
#  BVS.MCMC.bic                                                                 #
#                                                                              #
#  Slim clone of BVS.MCMC (ADS_sampler.R) that accepts a user-supplied target  #
#  function. Re-uses the identical ADS proposal, MH step, and summarisation,   #
#  so the only difference vs the proposed methods is the embedded marginal-   #
#  likelihood approximation.                                                   #
# ---------------------------------------------------------------------------- #
BVS.MCMC.bic <- function(model.par, mcmc.par, target.fn,
                         gamprior_type = "default",
                         verbose = FALSE) {

  ## ----- Initialisation ----- ##
  inits <- model.par$inits
  data  <- model.par$data
  data$xt2     <- t(data$xt) %*% data$xt
  data$x2      <- t(data$x)  %*% data$x
  data$xt2.inv <- solve(data$xt2)
  data$x2.inv  <- solve(data$x2)
  if (is.null(inits$n))           inits$n <- nrow(data$x)
  if (is.null(inits$n.threshold)) inits$n.threshold <-
      inits$n * (1 - 0.5 * (1 - mean(data$status)))

  n.iter <- mcmc.par$n.iter
  n.burn <- mcmc.par$n.burnin
  n.thin <- mcmc.par$n.thin
  if (n.thin < 1)            stop("n.thin must be integer >= 1")
  if (n.iter < n.burn)       stop("n.iter must exceed n.burnin")
  if (n.iter %% n.thin != 0) stop("n.iter must be divisible by n.thin")
  printevery <- mcmc.par$printevery
  n.chain    <- mcmc.par$n.chain
  n.save     <- floor((n.iter - n.burn) / n.thin)

  ## Globals required by MH.step.fixed / ADS.proposal.aft / prior.gam
  move.list  <<- data.frame(matrix(NA, nrow = 0, ncol = 2))
  colnames(move.list) <<- c("move", "accept")
  model.dict <<- data.frame(matrix(NA, nrow = 0, ncol = 2))
  colnames(model.dict) <<- c("gam", "postprob")
  gamprior   <<- gamprior_type
  gam.param  <<- model.par$gam.param

  gam.proposal <- ADS.proposal.aft
  MH.step      <- MH.step.fixed

  chain.list    <- list()
  haz.list      <- matrix(NA, nrow = n.save, ncol = n.chain)
  modelsize.df  <- data.frame(matrix(NA, nrow = n.iter, ncol = n.chain))
  postprob.list <- matrix(NA, nrow = n.save, ncol = n.chain)
  acc.list      <- data.frame(matrix(NA, nrow = n.iter, ncol = n.chain))

  ## ----- MCMC ----- ##
  for (chain in 1:n.chain) {
    par <- inits
    current.chain <- data.frame(matrix(NA, nrow = n.save,
                                       ncol = length(par$gam)))
    current.haz   <- rep(NA, n.save)
    iter.save <- 0
    for (iter in 1:n.iter) {
      update   <- MH.step(target.fn, proposal.dist = gam.proposal, par, data)
      par.new  <- update$par
      acc.list[iter, chain]      <- update$acc
      modelsize.df[iter, chain]  <- sum(par.new$gam != 0)
      postprob <- update$postprob

      if (iter > n.burn && iter %% n.thin == 0) {
        iter.save <- iter.save + 1
        current.chain[iter.save, ]       <- par.new$gam
        current.haz[iter.save]           <- gam2haz(par.new$gam)
        postprob.list[iter.save, chain]  <- postprob
      }
      if (verbose && iter %% printevery == 0) {
        cat(sprintf("Chain %d, Iter %d, Acc %.3f, size %.2f\n",
                    chain, iter,
                    mean(acc.list[(iter - printevery + 1):iter, chain]),
                    mean(modelsize.df[(iter - printevery + 1):iter, chain])))
      }
      par <- par.new
    }
    chain.list[[chain]] <- current.chain
    haz.list[, chain]   <- current.haz
  }

  list(chain        = chain.list,
       acc          = acc.list,
       modelsize    = modelsize.df,
       PIP          = PIP(chain.list),
       hazPIP       = haz.freq(haz.list),
       model.freq   = model.freq(chain.list),
       model.prob   = model.dict,                  ## (gam string, log-target)
       post.prob    = model.postprob(chain.list, postprob.list),
       postprob.raw = postprob.list)               ## per-iteration log-target
}

# ---------------------------------------------------------------------------- #
#  build_fullchain                                                              #
#                                                                              #
#  Produces a tidy per-iteration data.frame with columns                       #
#     gam       (model string, e.g. "0203010030")                              #
#     acc.gam   (acceptance indicator at the saved iteration)                  #
#     post.prob (log-target at the saved model)                                #
#                                                                              #
#  Works for BVS.MCMC output (via model.dict lookup) and BVS.MCMC.bic output   #
#  (which exposes postprob.raw directly).                                      #
# ---------------------------------------------------------------------------- #
build_fullchain <- function(chains, mcmcpar) {
  gam_strs <- apply(chains$chain[[1]], 1,
                    function(x) paste0(x, collapse = ""))
  acc_idx  <- seq(mcmcpar$n.burnin + mcmcpar$n.thin,
                  mcmcpar$n.iter,
                  mcmcpar$n.thin)
  acc_vec  <- as.numeric(chains$acc[acc_idx, 1])

  if (!is.null(chains$postprob.raw)) {
    pp_vec <- as.numeric(chains$postprob.raw[, 1])
  } else {
    ## BVS.MCMC path: look up per-iteration log-target via model.dict
    dict <- chains$model.prob
    pp_map <- setNames(as.numeric(dict$postprob),
                       as.character(dict$gam))
    pp_vec <- pp_map[gam_strs]
  }

  data.frame(gam       = gam_strs,
             acc.gam   = acc_vec,
             post.prob = pp_vec,
             stringsAsFactors = FALSE)
}
