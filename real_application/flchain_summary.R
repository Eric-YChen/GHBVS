# Usage (from real_application/): Rscript flchain_summary.R   (after flchain.R)
BVS_CORE <- normalizePath("../core")
suppressMessages({library(mvtnorm); library(survival); library(glmnet)})
source(file.path(BVS_CORE, "ADS_sampler.R"))
source(file.path(BVS_CORE, "utils", "map_estimates.R"))
source("prepare_data.R")

res.dir <- if (nzchar(Sys.getenv("BVS_QUICK"))) "results_quick" else "results"
priors <- c("LCM", "Product", "LCM-hyper")
res <- list()
for (d in c("orig", "trans")) for (p in priors){
  e <- new.env()
  load(file.path(res.dir, paste0("flchain_", d, "_", p, ".RData")), envir = e)
  ch <- e$chains
  vis <- setNames(ch$post.prob$visits, ch$post.prob$Vars)[ch$model.prob$gam]
  ch$post.prob <- model.postprob.dict(ch$model.prob$gam, ch$model.prob$postprob,
                                      visits = ifelse(is.na(vis), 0L, as.integer(vis)))
  ch$PIP.post <- PIP.postprob(ch$post.prob, ncol(e$dat$x)); ch$haz.post <- haz.postprob(ch$post.prob)
  res[[d]][[p]] <- list(chains = ch, dat = e$dat, runtime = e$runtime)
}

role <- c("0" = "", "1" = "^time", "2" = "^haz", "3" = "", "4" = "")
gam2label <- function(gs, vars){
  g <- strsplit(gs, "")[[1]]
  paste0("{", paste0(vars[g != "0"], role[g[g != "0"]], collapse = ", "), "}^", gam2haz(gs))
}
gam2size <- function(g){
  g <- as.integer(strsplit(g, "")[[1]]); sum(g %in% c(1, 2, 4)) + 2 * sum(g == 3)
}
hpd.set <- function(pp, level = 0.9){
  pp <- pp[order(pp$prob, decreasing = TRUE), ]
  pp[seq_len(which(cumsum(pp$prob) >= level)[1]), ]
}

sink(file.path(res.dir, "flchain_summary.txt"), split = TRUE)
for (d in c("orig", "trans")){
  vars <- colnames(res[[d]][[1]]$dat$x)
  cat("\n================ design:", d, "================\n")
  sets <- lapply(res[[d]], function(r) hpd.set(r$chains$post.prob))
  models <- unique(unlist(lapply(sets, function(s) s$Vars)))
  tab <- data.frame(model = sapply(models, gam2label, vars = vars),
                    effects = sapply(models, gam2size), row.names = NULL)
  for (p in priors){
    pp <- res[[d]][[p]]$chains$post.prob
    tab[[p]] <- round(pp$prob[match(models, pp$Vars)], 3)
    tab[[paste0(p, ".in90")]] <- models %in% sets[[p]]$Vars
  }
  cat("\n-- 90% HPP sets (renormalised over all evaluated models)\n"); print(tab)
  for (p in priors){
    ch <- res[[d]][[p]]$chains
    cat("\n--", p, "| runtime", round(as.numeric(res[[d]][[p]]$runtime), 1), "min",
        "| acceptance", round(mean(ch$acc[, 1]), 3), "| evaluated models", nrow(ch$post.prob), "\n")
    cat("hazard structure:\n"); print(round(ch$haz.post, 4))
    pip <- ch$PIP.post; colnames(pip) <- vars; cat("PIP:\n"); print(round(pip, 3))
    if (length(ch$opt.log)) { cat("optimisation problems:\n"); print(ch$opt.log) }
  }
  cat("\n-- MAP estimates (alpha, beta) of the models in the 90% sets\n")
  dat <- res[[d]][[1]]$dat
  for (gs in models){
    gam <- as.integer(strsplit(gs, "")[[1]]); used <- which(gam != 0)
    fmt <- function(M) paste0(vars[used], " (", sprintf("%.3f", M$alpha[used]), ", ",
                              sprintf("%.3f", M$beta[used]), ")", collapse = "; ")
    cat("\n", gam2label(gs, vars), "\n")
    cat("  LCM (g_E = 1): ", fmt(map.lcm(gam, dat)), "\n")
    cat("  Product       : ", fmt(map.product(gam, dat)), "\n")
    cat("  LCM (robust g): ", fmt(map.lcm(gam, dat, hyper = TRUE)), "\n")
  }
}

cat("\n================ data ================\n")
fl <- read.csv("data/flchain.csv", stringsAsFactors = FALSE)
dc <- fl[, c("age", "sex", "kappa", "lambda", "creatinine", "mgus", "futime", "death")]
cat("rows:", nrow(fl), "| complete cases:", sum(complete.cases(dc)), "\n")
dc <- dc[complete.cases(dc), ]; dc <- dc[dc$futime > 0, ]
cat("n =", nrow(dc), "| censored =", round(1 - mean(dc$death), 3),
    "| longest follow-up (years) =", round(max(dc$futime) / 365.24, 1), "\n")
cat("cor(kappa, lambda) =", round(cor(dc$kappa, dc$lambda), 3), "\n")
cat("cor(kappa + lambda, kappa / lambda) =", round(cor(dc$kappa + dc$lambda, dc$kappa / dc$lambda), 3), "\n")
sink()
