# Usage (from real_application/): Rscript nki_reference.R [mombf|BVSNLP]   (no argument: both)
suppressMessages({library(survival); library(parallel)})
source("prepare_data.R")
args <- commandArgs(trailingOnly = TRUE)
methods <- if (length(args) >= 1) args[1] else c("mombf", "BVSNLP")
stopifnot(all(methods %in% c("mombf", "BVSNLP")))
quick <- nzchar(Sys.getenv("BVS_QUICK"))
res.dir <- if (quick) "results_quick" else "results"
dir.create(res.dir, showWarnings = FALSE)
n.runs <- if (quick) 2 else 10
lse <- function(x){ m <- max(x); m + log(sum(exp(x - m))) }

dat <- prep.nki()
X <- dat$x

if ("mombf" %in% methods){
  suppressMessages(library(mombf))
  ys <- Surv(log(dat$y), event = dat$status)
  mdat <- data.frame(ys = ys, X, check.names = TRUE)
  xnames <- colnames(mdat)[-1]
  fy <- reformulate(xnames, response = "ys")
  reduce.ms <- function(m){
    S <- m$postSample[, match(xnames, colnames(m$postSample)), drop = FALSE]
    d <- data.frame(model = apply(S, 1, function(r) paste(colnames(X)[r == 1], collapse = ",")),
                    logpost = m$postProb, stringsAsFactors = FALSE)
    u <- d[!duplicated(d$model), ]
    u$visits <- as.vector(table(factor(d$model, levels = u$model)))
    list(models = u, margpp = m$margpp, postMode = m$postMode, postModeProb = m$postModeProb)
  }
  ms.list <- list(); time.ms <- numeric(0)
  for (i in seq_len(n.runs)){
    t0 <- Sys.time()
    m <- modelSelection(fy, data = mdat, priorCoef = zellnerprior(tau = nrow(X)),
                        priorDelta = modelbbprior(1, 1), enumerate = FALSE, method = "Laplace",
                        niter = if (quick) 2000 else 110000, burnin = if (quick) 1000 else 10000,
                        thinning = 1, center = FALSE, scale = FALSE, verbose = FALSE)
    ms.list[[i]] <- reduce.ms(m); rm(m)
    time.ms[i] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("mombf chain %d: %.1f min, %d distinct models\n", i, time.ms[i], nrow(ms.list[[i]]$models)))
    save(ms.list, time.ms, file = file.path(res.dir, "nki_ref_mombf.RData"))
  }
  allm <- do.call(rbind, lapply(ms.list, `[[`, "models"))
  lpv <- tapply(allm$logpost, allm$model, function(v) v[1])
  ms.models <- data.frame(model = names(lpv), logpost = as.numeric(lpv),
                          visits = as.numeric(tapply(allm$visits, allm$model, sum)[names(lpv)]),
                          stringsAsFactors = FALSE)
  ms.models$prob <- exp(ms.models$logpost - lse(ms.models$logpost))
  ms.models <- ms.models[order(-ms.models$prob), ]; rownames(ms.models) <- NULL
  save(ms.list, ms.models, time.ms, file = file.path(res.dir, "nki_ref_mombf.RData"))
  cat("AFT (block) Zellner, top models:\n"); print(head(ms.models[, c("model", "prob")], 10))
}

if ("BVSNLP" %in% methods){
  df <- data.frame(X, check.names = FALSE)
  t0 <- Sys.time()
  cl <- makeCluster(min(n.runs, detectCores()))
  bvs.list <- tryCatch({
    clusterEvalQ(cl, suppressMessages(library(BVSNLP)))
    clusterExport(cl, c("df", "dat", "quick"), envir = environment())
    parLapply(cl, seq_len(n.runs), function(i)
      bvs(X = df, resp = cbind(dat$y, dat$status), family = "survival", nlptype = "piMOM",
          mod_prior = "unif", niter = if (quick) 50 else 11000, inseed = i, ncpu = 1))
  }, finally = stopCluster(cl))
  time.bvs <- difftime(Sys.time(), t0, units = "mins")
  bvs.models <- do.call(rbind, lapply(bvs.list, function(b)
    data.frame(model = sapply(b$max_models, function(idx) paste(colnames(X)[sort(idx)], collapse = ",")),
               logpost = b$max_prob_vec, stringsAsFactors = FALSE)))
  bvs.models <- bvs.models[!duplicated(bvs.models$model), ]
  bvs.models$prob <- exp(bvs.models$logpost - lse(bvs.models$logpost))
  bvs.models <- bvs.models[order(-bvs.models$prob), ]; rownames(bvs.models) <- NULL
  bvs.nvis <- sapply(bvs.list, `[[`, "num_vis_models")
  bvs.list <- lapply(bvs.list, function(b) b[setdiff(names(b), "des_mat")])
  save(bvs.list, bvs.models, bvs.nvis, time.bvs, file = file.path(res.dir, "nki_ref_bvsnlp.RData"))
  cat("Cox piMOM, models evaluated per run:", bvs.nvis, "\n")
  cat("Cox piMOM, top models:\n"); print(head(bvs.models[, c("model", "prob")], 10))
}
