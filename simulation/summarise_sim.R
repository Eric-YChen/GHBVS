# Usage (from simulation/): Rscript summarise_sim.R
res.dir <- "results"; out.dir <- file.path(res.dir, "summary")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
ok <- function(r) is.list(r) && is.null(r[["error"]])
gam_vec <- function(s) as.integer(strsplit(s, "", fixed = TRUE)[[1]])
haz_of <- function(s) {
  g <- gam_vec(s)
  if (all(g == 0)) "Null" else if (all(g %in% 0:1)) "AH" else if (all(g %in% c(0, 2))) "PH" else
    if (all(g != 4)) "GH" else if (all(g %in% c(0, 4))) "AFT" else "Invalid"
}
sens_spec <- function(model, truth) {
  g <- gam_vec(model) != 0; t <- gam_vec(truth) != 0
  c(sens = sum(g & t) / sum(t), spec = sum(!g & !t) / sum(!t))
}
mse <- function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
load_res <- function(f) { e <- new.env(); load(f, envir = e); list(r = Filter(ok, e$results_parallel), d = e$design) }

# fit: summarise_fit() output; truth: true model string
fit_metrics <- function(fit, truth) {
  ss <- sens_spec(fit$hpm, truth); act <- gam_vec(truth) != 0
  hz <- fit$haz.post[c("AH", "PH", "AFT", "GH")]
  c(ss, haz = as.numeric(names(which.max(hz)) == haz_of(truth)),
    haz.freq = as.numeric(names(which.max(fit$haz.freq[c("AH", "PH", "AFT", "GH")])) == haz_of(truth)),
    pip.active = mean(fit$PIP.post[act]), pip.inactive = mean(fit$PIP.post[!act]))
}
scenario_table <- function(R, design, get_fit, runtime = NULL, init_time = NULL) {
  g <- design$str.grid
  rows <- lapply(seq_len(nrow(g)), function(i) {
    ri <- Filter(function(r) r$i == i, R); if (!length(ri)) return(NULL)
    M <- t(vapply(ri, function(r) fit_metrics(get_fit(r), g$true.mod[i]), numeric(6)))
    out <- data.frame(g[i, c("n", "p", "hstr", "baseline", "cens")], n.rep = nrow(M),
                      sens = mean(M[, "sens"]), sens.se = sd(M[, "sens"]) / sqrt(nrow(M)),
                      spec = mean(M[, "spec"]), spec.se = sd(M[, "spec"]) / sqrt(nrow(M)),
                      haz = mean(M[, "haz"]), haz.se = sd(M[, "haz"]) / sqrt(nrow(M)),
                      haz.freq = mean(M[, "haz.freq"]), pip.active = mean(M[, "pip.active"]),
                      pip.inactive = mean(M[, "pip.inactive"]),
                      events = mean(vapply(ri, function(r) r$cens_rat, numeric(1))))
    if (!is.null(runtime)) out$minutes <- mean(vapply(ri, runtime, numeric(1))) / 60
    if (!is.null(init_time)) out$init.minutes <- mean(vapply(ri, init_time, numeric(1))) / 60
    out
  })
  do.call(rbind, rows)
}
replicate_table <- function(R, design) {
  g <- design$str.grid
  sel <- vapply(R, function(r) { hz <- r$haz.post[c("AH", "PH", "AFT", "GH")]; names(hz)[which.max(hz)] }, character(1))
  pip <- t(vapply(R, function(r) { act <- which(gam_vec(g$true.mod[r$i]) != 0)
    c(r$PIP.post[act], mean(r$PIP.post[-act])) }, numeric(5)))
  i <- vapply(R, function(r) r$i, numeric(1))
  data.frame(g[i, c("n", "p", "hstr", "baseline", "cens")], j = vapply(R, function(r) r$j, numeric(1)), sel_hstr = sel,
             pip_1 = pip[, 1], pip_2 = pip[, 2], pip_3 = pip[, 3], pip_4 = pip[, 4], pip_inactive = pip[, 5], row.names = NULL)
}
main_outputs <- function(x, stem) {
  tab <- scenario_table(x$r, x$d, function(r) r, function(r) r$runtime[["sampler"]], function(r) r$runtime[["init"]])
  write.csv(tab, file.path(out.dir, paste0(stem, ".csv")), row.names = FALSE)
  write.csv(replicate_table(x$r, x$d), file.path(out.dir, paste0("replicates_", stem, ".csv")), row.names = FALSE)
}

convert_old <- function(r) {
  pp <- r$post.prob
  incl <- do.call(rbind, lapply(strsplit(pp$Vars, "", fixed = TRUE), function(s) s != "0"))
  w <- pp$prob; v <- pp$visits / sum(pp$visits)
  hz <- factor(pp$haz, levels = c("Null", "AH", "PH", "AFT", "GH"))
  by_haz <- function(wt) { h <- tapply(wt, hz, sum); h[is.na(h)] <- 0; h }
  list(i = r$i, j = r$j, cens_rat = r$cens_rat, hpm = pp$Vars[1], hpm.visited = pp$Vars[1],
       PIP.post = colSums(w * incl), PIP.freq = colSums(v * incl), haz.post = by_haz(w), haz.freq = by_haz(v),
       PIP.haz = t(sapply(c("AH", "PH", "AFT", "GH"), function(l) { s <- pp$haz == l
         if (!any(s)) rep(NA_real_, ncol(incl)) else colSums(w[s] * incl[s, , drop = FALSE]) / sum(w[s]) })),
       runtime = r$runtime)
}
load_main_old <- function(prior, p = 10) {
  fn <- file.path(res.dir, sprintf("sim_main_%s_p%d_default.RData", prior, p))
  fo <- file.path(res.dir, sprintf("mpi_simulation_raw_%s_p%d.RData", prior, p))
  if (file.exists(fn)) return(load_res(fn))
  if (!file.exists(fo)) return(NULL)
  e <- new.env(); load(fo, envir = e)
  g <- e$str.grid; for (k in names(g)) if (is.factor(g[[k]])) g[[k]] <- as.character(g[[k]])
  g$n <- as.integer(g$n); g$p <- as.integer(g$p)
  list(r = lapply(Filter(function(r) is.list(r) && !is.null(r$fullchain), e$results_parallel), convert_old),
       d = list(str.grid = g))
}
for (p in c(10, 50)) for (prior in c("LCM", "Product")) {
  if (!file.exists(file.path(res.dir, sprintf("sim_main_%s_p%d_default.RData", prior, p))) &&
      !is.null(x <- load_main_old(prior, p))) {
    main_outputs(x, sprintf("main_%s_p%d_default", prior, p))
    cat("summarised the finished p =", p, prior, "runs (earlier format)\n")
  }
}

files <- list.files(res.dir, pattern = "^sim_.*\\.RData$", full.names = TRUE)
for (f in files) {
  b <- sub("\\.RData$", "", basename(f)); x <- load_res(f)
  if (!length(x$r)) next
  if (startsWith(b, "sim_main_")) {
    main_outputs(x, sub("^sim_", "", b))
  } else if (startsWith(b, "sim_bic_")) {
    tab <- do.call(rbind, lapply(names(x$d$methods), function(m) {
      R <- Filter(function(r) !is.null(r$fits[[m]]), x$r)
      cbind(method = m, scenario_table(R, x$d, function(r) r$fits[[m]], function(r) r$runtime[[m]]))
    }))
    write.csv(tab, file.path(out.dir, paste0(sub("^sim_", "", b), ".csv")), row.names = FALSE)
  } else if (startsWith(b, "sim_freq_")) {
    g <- x$d$str.grid
    tab <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      ri <- Filter(function(r) r$i == i, x$r); if (!length(ri)) return(NULL)
      M <- t(vapply(ri, function(r) r$metrics[c("Sens", "Spec", "haz_correct")], numeric(3)))
      data.frame(g[i, c("n", "p", "hstr", "baseline", "cens")], n.rep = nrow(M),
                 sens = mean(M[, 1]), sens.se = sd(M[, 1]) / sqrt(nrow(M)),
                 spec = mean(M[, 2]), spec.se = sd(M[, 2]) / sqrt(nrow(M)),
                 haz = mean(M[, 3]), haz.se = sd(M[, 3]) / sqrt(nrow(M)))
    }))
    write.csv(tab, file.path(out.dir, paste0(sub("^sim_", "", b), ".csv")), row.names = FALSE)
  } else if (startsWith(b, "sim_conv_")) {
    g <- x$d$str.grid
    tab <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
      ri <- Filter(function(r) r$i == i, x$r); if (!length(ri)) return(NULL)
      v <- function(s, k) vapply(ri, function(r) as.numeric(r$gr[[s]][[k]]), numeric(1))
      cbind(g[i, c("n", "p", "hstr", "baseline", "cens")], n.rep = length(ri),
            do.call(cbind, lapply(names(x$d$gr_settings), function(s) {
              m <- rbind(mse(v(s, "mean")), mse(v(s, "max")), mse(v(s, "mean.practical")))
              setNames(as.data.frame(t(as.vector(t(m)))),
                       paste0(s, ".", c("mean", "mean.se", "max", "max.se", "mean.practical", "mean.practical.se")))
            })))
    }))
    write.csv(tab, file.path(out.dir, paste0(sub("^sim_", "", b), ".csv")), row.names = FALSE)
  }
  cat("summarised", basename(f), "\n")
}

fr <- file.path(res.dir, "sim_main_LCM_p10_restricted.RData")
if (file.exists(fr) && !is.null(M <- load_main_old("LCM", 10))) {
  R <- load_res(fr)
  key <- function(g, i) paste(g$n[i], g$hstr[i], g$baseline[i], g$cens[i])
  mkey <- vapply(M$r, function(r) paste(key(M$d$str.grid, r$i), r$j), character(1))
  pairs <- do.call(rbind, lapply(R$r, function(r) {
    g <- R$d$str.grid; m <- M$r[[match(paste(key(g, r$i), r$j), mkey)]]
    if (is.null(m)) return(NULL)
    data.frame(g[r$i, c("n", "hstr", "baseline", "cens")], j = r$j, variable = seq_along(r$PIP.post),
               active = gam_vec(g$true.mod[r$i]) != 0,
               pip.full.given.true = m$PIP.haz[g$hstr[r$i], ], pip.restricted = r$PIP.post,
               haz.full.true = m$haz.post[[g$hstr[r$i]]], row.names = NULL)
  }))
  write.csv(pairs, file.path(out.dir, "restricted_pairs.csv"), row.names = FALSE)
  cat("summarised restricted vs full model space\n")
}
