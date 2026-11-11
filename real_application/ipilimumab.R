# Usage (from real_application/): Rscript ipilimumab.R
BVS_CORE <- normalizePath("../core")
suppressMessages({library(mvtnorm); library(survival); library(glmnet)})
source(file.path(BVS_CORE, "ADS_sampler.R"))
source(file.path(BVS_CORE, "utils", "EB_gprior_fixed.R"))
source(file.path(BVS_CORE, "utils", "map_estimates.R"))
source("prepare_data.R")
quick <- nzchar(Sys.getenv("BVS_QUICK"))
res.dir <- if (quick) "results_quick" else "results"
dir.create(res.dir, showWarnings = FALSE)

specs <- list(GH = c(1, 1, 1, 1), PH = c(0, 1, 0, 0), AFT = c(0, 0, 0, 1))
gp <- function(h) list(layer1.betabin.par = c(1, 1), layer2.mulnom.par = h, layer3.beta.par = 1/3)
mcmcpar <- if (quick) list(n.iter = 2000, n.burnin = 1000, n.thin = 2, n.chain = 1, printevery = 2000) else
  list(n.iter = 20000, n.burnin = 10000, n.thin = 2, n.chain = 1, printevery = 20000)

role <- c("1" = "^time", "2" = "^haz", "3" = "", "4" = "")
label <- function(s, vars){
  g <- strsplit(s, "")[[1]]; i <- which(g != "0")
  if (!length(i)) return("null")
  paste0("{", paste0(vars[i], role[g[i]], collapse = ", "), "}^", gam2haz(s))
}
base.par <- function(dat, gam){
  n <- nrow(dat$x); p <- ncol(dat$x)
  list(nu = 0, theta0 = 0, theta = rep(0, p), eta = rep(0, p), g = n, n = n,
       n.events = sum(dat$status), a_g = NULL, gam = gam, nu.alpha = 0.01, nu.beta = 0.01)
}

analyse <- function(dat){
  p <- ncol(dat$x); vars <- colnames(dat$x)
  grid <- as.matrix(expand.grid(rep(list(0:4), p)))
  grid <- grid[apply(grid, 1, gam2haz) != "Invalid model", , drop = FALSE]
  models <- apply(grid, 1, paste0, collapse = "")
  admit.null <<- TRUE
  logml <- apply(grid, 1, function(g) marg.loglik(base.par(dat, g), dat))
  gamprior <<- "custom"
  exact <- sapply(specs, function(h){
    gam.param <<- gp(h)
    lpr <- apply(grid, 1, function(g) prior.gam(base.par(dat, g), log = TRUE))
    lp <- lpr + logml; w <- exp(lp - max(lp)); setNames(w / sum(w), models)
  })
  mcmc <- lapply(specs, function(h){
    mp <- list(inits = c(base.par(dat, rep(0, p))[c("nu", "theta0", "theta", "eta")],
                         list(g = 1, a_g = NULL, gam = rep(0, p), nu.alpha = 0.01, nu.beta = 0.01)),
               gam.param = gp(h), data = dat)
    BVS.MCMC(mp, mcmcpar, gprior.mat = "LCM", gamprior = "custom", null.model = TRUE, verbose = FALSE)
  })
  freq <- sapply(mcmc, function(fit){
    f <- setNames(rep(0, length(models)), models); f[names(fit$model.freq)] <- as.numeric(fit$model.freq); f
  })
  renorm <- sapply(mcmc, function(fit){
    f <- setNames(rep(0, length(models)), models); f[fit$post.prob$Vars] <- fit$post.prob$prob; f
  })
  logBF <- setNames(logml - logml[models == paste(rep(0, p), collapse = "")], models)
  hz <- apply(exact, 2, function(w) tapply(w, factor(sapply(models, gam2haz),
                                                     levels = c("Null", "AH", "PH", "AFT", "GH")), sum))
  pip <- matrix(sapply(colnames(exact), function(s) sapply(seq_len(p), function(j)
    sum(exact[grid[, j] != 0, s]))), nrow = p, dimnames = list(vars, colnames(exact)))
  list(grid = grid, models = models, vars = vars, logml = logml, exact = exact,
       freq = freq, renorm = renorm, mcmc = mcmc, logBF = logBF, hz = hz, pip = pip)
}

set.seed(1)
S1 <- prep.ipilimumab(noise = FALSE)
S2 <- prep.ipilimumab(noise = TRUE, noise.seed = 1)
A1 <- analyse(S1$dat)
A2 <- analyse(S2$dat)

km    <- survfit(Surv(time, status) ~ trt, data = S1$df)
cox1  <- coxph(Surv(time, status) ~ trt, data = S1$df)
aft1  <- survreg(Surv(time, status) ~ trt, data = S1$df, dist = "lognormal")
cox2  <- coxph(Surv(time, status) ~ trt + spur, data = S2$df)
aft2  <- survreg(Surv(time, status) ~ trt + spur, data = S2$df, dist = "lognormal")

map1 <- lapply(c(null = 0, AH = 1, PH = 2, GH = 3, AFT = 4), function(k) map.lcm(k, S1$dat))
S.bma <- function(t, x){
  w <- A1$exact[, "GH"]
  Reduce(`+`, lapply(names(map1), function(m){
    M <- map1[[m]]; w[[as.character(M$gam)]] * S.gh.lnorm(t, x, M$mu, M$sigma, M$alpha, M$beta)
  }))
}

save(S1, S2, A1, A2, km, cox1, aft1, cox2, aft2, map1, file = file.path(res.dir, "ipilimumab.RData"))

fmt.models <- function(A, s, level = 0.9){
  w <- sort(A$exact[, s], decreasing = TRUE); k <- which(cumsum(w) >= level)[1]
  paste0(sapply(names(w)[1:k], label, vars = A$vars), " (", sprintf("%.3f", w[1:k]), ")", collapse = "; ")
}
sink(file.path(res.dir, "ipilimumab_summary.txt"), split = TRUE)
cat("n =", nrow(S1$df), "| events =", sum(S1$df$status), "| censored =",
    round(100 * mean(S1$df$status == 0), 1), "% | longest follow-up =", max(S1$df$time), "months\n")
for (k in 1:2){
  A <- list(A1, A2)[[k]]
  cat("\n================ Study", k, if (k == 1) "(treatment only)" else "(treatment + noise)", "================\n")
  if (k == 1){ cat("\n-- Cox PH\n"); print(summary(cox1)$coefficients)
               cat("\n-- log-normal AFT\n"); print(summary(aft1)$table) }
  else       { cat("\n-- Cox PH\n"); print(summary(cox2)$coefficients)
               cat("\n-- log-normal AFT\n"); print(summary(aft2)$table) }
  for (s in names(specs)){
    cat("\n--", s, "analysis\n")
    cat("90% HPM set:", fmt.models(A, s), "\n")
    cat("PIP:", paste(A$vars, sprintf("%.3f", A$pip[, s]), collapse = ", "), "\n")
    cat("hazard structure (Null, AH, PH, AFT, GH):", sprintf("%.3f", A$hz[, s]), "\n")
    cat("max |MCMC frequency - exact|:", sprintf("%.4f", max(abs(A$freq[, s] - A$exact[, s]))),
        "| max |renormalised MCMC - exact|:", sprintf("%.4f", max(abs(A$renorm[, s] - A$exact[, s]))), "\n")
  }
  cat("\n-- log Bayes factors against the null model\n")
  print(round(sort(setNames(A$logBF, sapply(A$models, label, vars = A$vars)), decreasing = TRUE), 2))
}
cat("\n-- MAP (mu, sigma, alpha, beta) under the LCM prior, treatment-only models\n")
print(round(t(sapply(map1, function(M) c(mu = M$mu, sigma = M$sigma, alpha = M$alpha[[1]], beta = M$beta[[1]]))), 3))
sink()

if (!quick){
  dir.create("../tables", showWarnings = FALSE)
  tex.model <- function(s, vars, structure){
    g <- strsplit(s, "")[[1]]; i <- which(g != "0")
    if (!length(i)) return("Null")
    v <- ifelse(vars[i] == "trt", "\\textbf{trt}", vars[i])
    paste0("\\{", paste(v, collapse = ", "), "\\}", if (structure == "GH") paste0("$^{", gam2haz(s), "}$") else "")
  }
  sets <- lapply(names(specs), function(s){
    w <- sort(A2$exact[, s], decreasing = TRUE); k <- which(cumsum(w) >= 0.9)[1]
    sprintf("%d: %s (%.3f)", 1:k, sapply(names(w)[1:k], tex.model, vars = A2$vars, structure = s), w[1:k])
  })
  K <- max(lengths(sets)); sets <- lapply(sets, function(x) c(x, rep("", K - length(x))))
  hpm <- sapply(seq_len(K), function(k) paste(c(if (k == 1) sprintf("\\multirow{%d}{*}{90\\%% HPM set}", K) else "",
                                                 sapply(sets, `[`, k)), collapse = " & "))
  pip <- paste(c("$(\\text{PIP}_\\text{trt},\\text{PIP}_\\text{spur})$",
                 sapply(names(specs), function(s) sprintf("$(%.3f,%.3f)$", A2$pip["trt", s], A2$pip["spur", s]))), collapse = " & ")
  writeLines(c("\\begin{table}[h]", "\\centering",
    paste0("\\caption{$90\\%$ highest probability models sets and posterior inclusion probabilities of the variables under each ",
           "hazard structure assumption using the proposed methodology under LCM prior. }"),
    "\\label{tab:motiv_hpm}", "\\begin{tabular}{@{}cccc@{}}", "\\toprule",
    "Summary & GH & PH & AFT \\\\ \\midrule", paste(hpm, "\\\\"), "\\midrule", paste(pip, "\\\\ \\bottomrule"),
    "\\end{tabular}", "\\end{table}"), "../tables/tab_motiv_hpm.tex")
  cc <- summary(cox2)$coefficients; at <- summary(aft2)$table
  writeLines(c("\\begin{table}[h]", "\\centering",
    "\\caption{Fitted coefficients and p-values under semi-parametric PH and log-normal AFT model.}",
    "\\label{tab:motiv_standard}", "\\begin{tabular}{@{}ccccc@{}}", "\\toprule",
    "Model & coef$_\\text{trt}$ & p-value & coef$_\\text{spur}$ & p-value \\\\ \\midrule",
    sprintf("Cox PH & %.3f & %.3f & %.3f & %.3f \\\\", cc["trt", "coef"], cc["trt", "Pr(>|z|)"], cc["spur", "coef"], cc["spur", "Pr(>|z|)"]),
    sprintf("AFT & %.3f & %.3f & %.3f & %.3f \\\\ \\bottomrule", at["trt", "Value"], at["trt", "p"], at["spur", "Value"], at["spur", "p"]),
    "\\end{tabular}", "\\end{table}"), "../tables/tab_motiv_standard.tex")
  cat("written ../tables/tab_motiv_hpm.tex and ../tables/tab_motiv_standard.tex\n")
}

col.arm <- c(placebo = "#2a78d6", ipilimumab = "#eb6834")
ink <- "#0b0b0b"; ink2 <- "#52514e"; rule <- "#d9d8d4"
tr  <- function(col, a) adjustcolor(col, alpha.f = a)
tmax <- 36; tg <- seq(0.01, tmax, length.out = 600)
cox.sf <- survfit(cox1, newdata = data.frame(trt = 0:1))
aft.lp <- predict(aft1, newdata = data.frame(trt = 0:1), type = "lp")
p.eff <- function(s) 1 - A1$exact["0", s]
fmt.p <- function(x) if (x > 0.999) "P(effect) > 0.999" else sprintf("P(effect) = %.2f", x)
fits <- list(
  list(title = "(a) Proportional hazards",
       S0 = function(t) summary(cox.sf, times = t, extend = TRUE)$surv[, 1],
       S1 = function(t) summary(cox.sf, times = t, extend = TRUE)$surv[, 2],
       note = c(sprintf("Cox: HR %.2f (p = %.3f)", exp(coef(cox1)), summary(cox1)$coefficients[, "Pr(>|z|)"]),
                paste("PH only:", fmt.p(p.eff("PH"))))),
  list(title = "(b) Accelerated failure time",
       S0 = function(t) plnorm(t, aft.lp[1], aft1$scale, lower.tail = FALSE),
       S1 = function(t) plnorm(t, aft.lp[2], aft1$scale, lower.tail = FALSE),
       note = c(sprintf("Log-normal: TR %.2f (p = %.2f)", exp(coef(aft1)[2]), summary(aft1)$table["trt", "p"]),
                paste("AFT only:", fmt.p(p.eff("AFT"))))),
  list(title = "(c) Joint selection (GH family)",
       S0 = function(t) S.bma(t, 0), S1 = function(t) S.bma(t, 1),
       note = c(fmt.p(p.eff("GH")),
                sprintf("AH %.2f, GH %.2f", A1$hz["AH", "GH"], A1$hz["GH", "GH"]))))
draw <- function(){
  op <- par(mfrow = c(1, 3), mar = c(3.2, 3.2, 2.2, 0.6), oma = c(1.4, 0, 0, 0),
            mgp = c(2, 0.55, 0), tcl = -0.25, cex = 0.72, col.axis = ink2, col.lab = ink, fg = ink2, family = "sans")
  on.exit(par(op))
  for (k in seq_along(fits)){
    f <- fits[[k]]
    plot(NA, xlim = c(0, tmax), ylim = c(0, 1), xaxs = "i", yaxs = "i",
         xlab = "Months since randomisation", ylab = if (k == 1) "Overall survival" else "", axes = FALSE)
    abline(h = seq(0.2, 0.8, 0.2), col = rule, lwd = 0.6)
    axis(1, at = seq(0, 36, 6), lwd = 0.6); axis(2, at = seq(0, 1, 0.2), las = 1, lwd = 0.6)
    box(bty = "l", lwd = 0.6)
    mtext(f$title, side = 3, line = 0.6, adj = 0, cex = 0.78, font = 2, col = ink)
    lines(km[1], col = tr(col.arm[1], 0.55), lwd = 1, conf.int = FALSE, mark.time = FALSE)
    lines(km[2], col = tr(col.arm[2], 0.55), lwd = 1, conf.int = FALSE, mark.time = FALSE)
    lines(tg, f$S0(tg), col = col.arm[1], lwd = 2)
    lines(tg, f$S1(tg), col = col.arm[2], lwd = 2)
    text(tmax - 0.8, c(0.93, 0.84), f$note, adj = c(1, 1), cex = 0.88, col = ink)
    if (k == 1){
      text(21, 0.40, "Ipilimumab", col = ink, adj = c(0, 0), cex = 0.88)
      text(13, 0.17, "Placebo", col = ink, adj = c(0, 1), cex = 0.88)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.85, text.col = ink,
         legend = c("Placebo", "Ipilimumab", "Kaplan-Meier (thin)", "Model fit (thick)"),
         col = c(col.arm, tr(ink2, 0.55), ink2), lwd = c(2, 2, 1, 2), lty = 1, seg.len = 1.6,
         x.intersp = 0.6, text.width = NA)
}
pdf(file.path(res.dir, "fig_ipilimumab.pdf"), width = 7.2, height = 2.9)
draw(); invisible(dev.off())
