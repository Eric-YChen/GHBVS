# Usage (from real_application/): Rscript nki_summary.R   (after nki_proposed.R and nki_reference.R)
BVS_CORE <- normalizePath("../core")
source(file.path(BVS_CORE, "utils", "helper_functions.R"))
res.dir <- if (nzchar(Sys.getenv("BVS_QUICK"))) "results_quick" else "results"
lse <- function(x){ m <- max(x); m + log(sum(exp(x - m))) }
role <- c("1" = "^time", "2" = "^haz", "3" = "", "4" = "")

summarise.proposed <- function(prior){
  e <- new.env(); load(file.path(res.dir, paste0("nki_", prior, ".RData")), envir = e)
  vars <- colnames(e$dat$x)
  ok <- sapply(e$chains.list, function(x) is.null(x$error))
  ch <- e$chains.list[ok]
  D <- lapply(ch, function(x) x$model.prob[is.finite(x$model.prob$postprob), ])
  alld <- do.call(rbind, D); alld <- alld[!duplicated(alld$gam), ]
  res <- data.frame(model = alld$gam, lp = alld$postprob, stringsAsFactors = FALSE)
  res$pp <- exp(res$lp - lse(res$lp))
  res$haz <- vapply(res$model, gam2haz, character(1))
  res <- res[order(-res$pp), ]; rownames(res) <- NULL
  allv <- do.call(rbind, lapply(ch, function(x) x$post.prob[x$post.prob$visits > 0, c("Vars", "logpost", "visits")]))
  vis <- tapply(allv$visits, allv$Vars, sum)
  res$mf <- as.numeric(vis[res$model]) / sum(allv$visits); res$mf[is.na(res$mf)] <- 0
  lpv <- tapply(allv$logpost, allv$Vars, function(v) v[1])
  res$pv <- exp(as.numeric(lpv[res$model]) - lse(lpv)); res$pv[is.na(res$pv)] <- 0
  top <- seq_len(min(50, nrow(res)))
  res$vars <- NA_character_; res$label <- NA_character_
  res$vars[top] <- sapply(res$model[top], function(s){ g <- strsplit(s, "")[[1]]; paste(vars[g != "0"], collapse = ",") })
  res$label[top] <- sapply(res$model[top], function(s){ g <- strsplit(s, "")[[1]]; i <- which(g != "0")
    paste0("{", paste0(vars[i], role[g[i]], collapse = ", "), "}^", gam2haz(s)) })
  per.chain <- sapply(D, function(d){
    x <- setNames(exp(d$postprob - lse(d$postprob)), d$gam)[res$model[1:10]]; x[is.na(x)] <- 0; x })
  pip <- vapply(seq_along(vars), function(j) sum(res$pp[substr(res$model, j, j) != "0"]), numeric(1))
  names(pip) <- vars
  list(res = res, per.chain = per.chain, vars = vars, n.ok = sum(ok), n = length(ok),
       n.visited = length(vis), runtime = sapply(ch, function(x) as.numeric(x$runtime)),
       acc = sapply(ch, function(x) mean(x$acc[, 1])), pip = pip,
       opt.log = unique(unlist(lapply(ch, function(x) names(x$opt.log)))))
}
hz.probs <- function(res){
  h <- tapply(res$pp, factor(res$haz, levels = c("AH", "PH", "AFT", "GH")), sum); h[is.na(h)] <- 0; h
}
hpd.size <- function(p) which(cumsum(sort(p, decreasing = TRUE)) >= 0.5)[1]

P <- list(LCM = summarise.proposed("LCM"), Product = summarise.proposed("Product"))
e <- new.env()
load(file.path(res.dir, "nki_ref_mombf.RData"), envir = e)
load(file.path(res.dir, "nki_ref_bvsnlp.RData"), envir = e)
ref <- list(AFT.Zellner = e$ms.models, Cox.piMOM = e$bvs.models)
ms.per.chain <- sapply(e$ms.list, function(m){
  d <- m$models; v <- setNames(exp(d$logpost - lse(d$logpost)), d$model)[ref$AFT.Zellner$model[1:10]]
  v[is.na(v)] <- 0; v })

sink(file.path(res.dir, "nki_summary.txt"), split = TRUE)
for (s in names(P)){
  S <- P[[s]]
  cat("\n================", s, "prior ================\n")
  cat("chains:", S$n.ok, "of", S$n, "| distinct evaluated models:", nrow(S$res), "| distinct visited after burn-in:", S$n.visited,
      "| runtime per chain (min):", round(S$runtime), "| acceptance:", round(S$acc, 3), "\n")
  top <- head(S$res, 15); top[c("pp", "pv", "mf")] <- lapply(top[c("pp", "pv", "mf")], round, 3)
  print(top[, c("label", "pp", "pv", "mf")], row.names = FALSE)
  cat("per-chain range of the top-10 probabilities:", round(max(apply(S$per.chain, 1, function(v) diff(range(v)))), 4), "\n")
  cat("hazard structure (AH, PH, AFT, GH):", round(hz.probs(S$res), 3), "\n")
  cat("50% HPP set size:", hpd.size(S$res$pp), "\n")
  cat("PIPs > 0.1:\n"); print(round(sort(S$pip[S$pip > 0.1], decreasing = TRUE), 3))
  if (length(S$opt.log)) cat("models with optimisation problems:", S$opt.log, "\n")
}
for (s in names(ref)){
  cat("\n================", s, "================\n")
  cat("distinct models:", nrow(ref[[s]]), "\n")
  print(head(ref[[s]][, c("model", "prob")], 10), row.names = FALSE)
  cat("50% HPP set size:", hpd.size(ref[[s]]$prob), "\n")
  if (s == "AFT.Zellner")
    cat("per-chain range of the top-10 probabilities:", round(max(apply(ms.per.chain, 1, function(v) diff(range(v)))), 4),
        "| runtime per chain (min):", round(e$time.ms), "\n")
  if (s == "Cox.piMOM")
    cat("models evaluated per run:", e$bvs.nvis, "| wall time (min):", round(as.numeric(e$time.bvs), 1), "\n")
}

cat("\n================ Table 2: rank (probability) of the top LCM models ================\n")
rank.of <- function(probs, key, keys){
  o <- order(-probs); r <- match(key, keys[o])
  if (is.na(r)) "-" else sprintf("%d (%.3f)", r, probs[o][r])
}
tab <- do.call(rbind, lapply(1:6, function(i){
  m <- P$LCM$res[i, ]
  data.frame(model = m$label,
             LCM = rank.of(P$LCM$res$pp, m$model, P$LCM$res$model),
             Product = rank.of(P$Product$res$pp, m$model, P$Product$res$model),
             AFT.Zellner = if (m$haz == "AFT") rank.of(ref$AFT.Zellner$prob, m$vars, ref$AFT.Zellner$model) else "-",
             stringsAsFactors = FALSE)
}))
tab$Cox.piMOM <- sprintf("%d (%.3f): {%s}", 1:6, head(ref$Cox.piMOM$prob, 6), head(ref$Cox.piMOM$model, 6))
print(tab, row.names = FALSE)
cat("\nProduct prior, top 6 models:\n")
print(head(P$Product$res[, c("label", "pp")], 6), row.names = FALSE)
sink()
