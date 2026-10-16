# Usage (from simulation/): Rscript make_figures.R   (after summarise_sim.R)
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(grid)
  library(dplyr); library(tidyr); library(ggh4x) })
source("utils/figure_style.R")
S <- file.path("results", "summary"); OUT <- file.path("results", "figures")
fig_dir <- function(d) { p <- file.path(OUT, d); dir.create(p, showWarnings = FALSE, recursive = TRUE); p }
cens_pct <- function(x) ifelse(x == "0.25default", 25L, 50L)

cells <- function(files, method = NULL) {
  x <- do.call(rbind, lapply(file.path(S, files), read.csv, stringsAsFactors = FALSE))
  if (!is.null(method)) x <- x[x$method == method, ]
  data.frame(n = x$n, p = x$p, hstr = x$hstr, baseline = x$baseline, c = cens_pct(x$cens), n.rep = x$n.rep,
             model_sensitivity = x$sens, model_specificity = x$spec, hstr_sensitivity = x$haz,
             model_sensitivity_se = x$sens.se, model_specificity_se = x$spec.se, hstr_sensitivity_se = x$haz.se,
             chain_time = if (is.null(x$minutes)) NA else x$minutes,
             init_time = if (is.null(x$init.minutes)) NA else x$init.minutes, stringsAsFactors = FALSE)
}
three <- function(d, label, dir, hstr_file = "plot_hstr_sensitivity.pdf") {
  create_heatmap(d, "model_sensitivity", sprintf("Mean Model Sensitivity (%s)", label), dir, "plot_model_sensitivity.pdf")
  create_heatmap(d, "model_specificity", sprintf("Mean Model Specificity (%s)", label), dir, "plot_model_specificity.pdf",
                 col.values = scales::rescale(c(0.9, 0.997, 1), to = c(0, 1)))
  create_heatmap(d, "hstr_sensitivity", sprintf("Hazard Structure Sensitivity (%s)", label), dir, hstr_file)
}
time_plot <- function(d, metric, title, dir, file)
  create_heatmap(d, metric, title, dir, file, reverse_scale = TRUE,
                 custom_limits = c(0, max(d[[metric]], na.rm = TRUE)), col.values = c(0, 0.5, 1))

for (pr in c("LCM", "Product")) {
  dir <- fig_dir(if (pr == "LCM") "EB-sims" else "FB-sims")
  d <- cells(sprintf("main_%s_p%d_default.csv", pr, c(10, 50, 100)))
  cat(sprintf("=== main study, %s: %d scenarios ===\n", pr, nrow(d)))
  three(d, pr, dir)
  create_heatmap_se(d, "model_sensitivity_se", sprintf("Model Sensitivity SE (%s)", pr), dir, "plot_model_sensitivity_SE.pdf", se_rule = "simple")
  create_heatmap_se(d, "model_specificity_se", sprintf("Model Specificity SE (%s)", pr), dir, "plot_model_specificity_SE.pdf", se_rule = "simple")
  time_plot(d, "chain_time", "Mean Chain Running Time (minutes)", dir, "plot_chain_time.pdf")
  if (pr == "LCM") time_plot(d, "init_time", "Mean Initialization Time (minutes)", dir, "plot_init_time.pdf")
}

gr <- do.call(rbind, lapply(c(10, 50, 100), function(pp) {
  f <- file.path(S, sprintf("conv_p%d.csv", pp)); if (!file.exists(f)) return(NULL); read.csv(f, stringsAsFactors = FALSE) }))
gr$c <- cens_pct(gr$cens)
dir <- fig_dir("EB-sims")
for (w in list(c("s50k", "", "plot_gelman_rubin"), c("s20k", ", first 20,000 iterations", "plot_gelman_rubin_20k"))) {
  g <- gr; g$psrf <- g[[paste0(w[1], ".mean")]]; g$psrf_se <- g[[paste0(w[1], ".mean.se")]]
  create_heatmap(g, "psrf", sprintf("Mean Gelman-Rubin Statistic%s (LCM)", w[2]), dir, paste0(w[3], ".pdf"),
                 reverse_scale = TRUE, col.values = c(0, 0.5, 1), custom_limits = c(1, max(g$psrf)), digits = 4, text_size = 5.4)
  create_heatmap_se(g, "psrf_se", sprintf("Gelman-Rubin Statistic SE%s (LCM)", w[2]), dir, paste0(w[3], "_SE.pdf"),
                    digits = 5, text_size = 5.0)
}

three(cells(sprintf("freq_bic_p%d.csv", c(10, 50, 100))), "Frequentist BIC", fig_dir("freq-sims"), "plot_hstr_accuracy.pdf")
have <- function(f) f[file.exists(file.path(S, f))]
three(cells(have(c("bic_p10.csv", "bic_p10_n500.csv")), "bic_flat"), "BIC-flat", fig_dir("BIC-unif_sims"))
bg <- cells(have(c("bic_p10.csv", "bic_p10_n500.csv", "bic_p50.csv", "bic_p50_n500.csv")), "bic_gamma"); dir <- fig_dir("BIC_sims")
three(bg, "BIC-gam", dir)
create_heatmap_se(bg, "model_sensitivity_se", "Model Sensitivity SE (BIC-gam)", dir, "plot_model_sensitivity_SE.pdf")
create_heatmap_se(bg, "model_specificity_se", "Model Specificity SE (BIC-gam)", dir, "plot_model_specificity_SE.pdf")
create_heatmap_se(bg, "hstr_sensitivity_se", "Hazard Structure Sensitivity SE (BIC-gam)", dir, "plot_hstr_sensitivity_SE.pdf")
uf <- sprintf("main_LCM_p%d_unif.csv", c(10, 50)); uf <- uf[file.exists(file.path(S, uf))]
three(cells(uf), "LCM - L1_Unif", fig_dir("EB-unif-sims"))

HSTR_LEV <- c("AH", "PH", "AFT", "GH")
facet_cols <- function(z) z %>% mutate(
  row_label = factor(paste0("p=", p), levels = paste0("p=", c(10, 50, 100))),
  row_group = factor(paste0("n=", n), levels = paste0("n=", c(250, 500, 1000))),
  col_label = factor(paste0(baseline, " ", c, "%"), levels = c("LN 25%", "LN 50%", "PGW 25%", "PGW 50%")),
  col_group = factor(hstr, levels = HSTR_LEV))
strip <- strip_nested(
  background_x = list(element_rect(fill = c("#FFE4E4", "#E4E4FF", "#E4FFE4", "#FFFFE4"), color = "gray40"),
                      element_rect(fill = "gray95", color = "gray70")),
  background_y = list(element_rect(fill = "gray95", color = "gray70"),
                      element_rect(fill = c("#FFF0E0", "#E0F0FF", "#F0FFE0"), color = "gray40")),
  text_x = list(element_text(size = 12, face = "bold"), element_text(size = 10, face = "bold")),
  text_y = list(element_text(size = 10, face = "bold"), element_text(size = 12, face = "bold")))
panel_theme <- function(x_text) theme_bw(base_size = 9) + theme(
  strip.background = element_blank(), axis.text.x = x_text, axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 11, face = "bold", margin = margin(0, 8, 0, 0)),
  panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(linewidth = 0.2, color = "gray90"),
  panel.spacing.x = unit(0.15, "lines"), panel.spacing.y = unit(0.15, "lines"),
  panel.border = element_rect(linewidth = 0.3, color = "gray50"), legend.position = "bottom",
  legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9),
  legend.key.size = unit(0.5, "cm"), plot.margin = margin(10, 10, 10, 10))
save10 <- function(p, path) { ggsave(path, p, width = 10, height = 10, units = "in", dpi = 300); cat("  saved:", path, "\n") }
var_colors <- c("X[1]" = "#E41A1C", "X[2]" = "#377EB8", "X[3]" = "#4DAF4A", "X[4]" = "#984EA3", "bar(X)[inactive]" = "#999999")
var_labels <- c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]), expression(bar(X)[inactive]))

for (pr in c("LCM", "Product")) {
  dir <- fig_dir(if (pr == "LCM") "EB-sims" else "FB-sims")
  d <- do.call(rbind, lapply(sprintf("replicates_main_%s_p%d_default.csv", pr, c(10, 50, 100)),
                             function(f) read.csv(file.path(S, f), stringsAsFactors = FALSE)))
  d$c <- cens_pct(d$cens)
  props <- d %>% count(p, n, hstr, baseline, c, sel_hstr, name = "k") %>% group_by(p, n, hstr, baseline, c) %>%
    mutate(Proportion = k / sum(k)) %>% ungroup() %>%
    complete(nesting(p, n, hstr, baseline, c), sel_hstr = HSTR_LEV, fill = list(k = 0, Proportion = 0)) %>%
    mutate(Hstr_Selected = factor(sel_hstr, levels = HSTR_LEV), x_dummy = "") %>% facet_cols()
  g <- ggplot(props, aes(x = x_dummy, y = Proportion, fill = Hstr_Selected)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.85) +
    geom_text(aes(label = ifelse(Proportion >= 0.10, sprintf("%.0f", Proportion * 100), "")),
              position = position_stack(vjust = 0.5), size = 2, fontface = "bold", color = "white") +
    facet_nested(rows = vars(row_label, row_group), cols = vars(col_group, col_label),
                 nest_line = element_line(colour = "gray50", linewidth = 0.5), strip = strip) +
    scale_fill_manual(values = c(AH = "#E41A1C", PH = "#377EB8", AFT = "#4DAF4A", GH = "#984EA3"), drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(x = NULL, y = "Proportion", fill = "Selected Hazard Structure") +
    panel_theme(element_blank()) + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_line(linewidth = 0.3)) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
  save10(g, file.path(dir, "hstr_proportion_stacked.pdf"))
  long <- d %>% select(p, n, hstr, baseline, c, pip_1, pip_2, pip_3, pip_4, pip_inactive) %>%
    pivot_longer(cols = starts_with("pip"), names_to = "Variable", values_to = "Inclusion") %>%
    mutate(Inclusion = pmin(pmax(Inclusion, 0), 1)) %>%
    mutate(Variable = factor(Variable, levels = c("pip_1", "pip_2", "pip_3", "pip_4", "pip_inactive"),
                             labels = names(var_colors))) %>% facet_cols()
  g <- ggplot(long, aes(x = Variable, y = Inclusion, fill = Variable)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8, linewidth = 0.3) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.5, linewidth = 0.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.3) +
    facet_nested(rows = vars(row_label, row_group), cols = vars(col_group, col_label),
                 nest_line = element_line(colour = "gray50", linewidth = 0.5), strip = strip) +
    scale_fill_manual(values = var_colors, labels = var_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_x_discrete(labels = var_labels) + labs(x = NULL, y = "Inclusion Proportion", fill = "Variable") +
    panel_theme(element_text(size = 6, angle = 45, hjust = 1, vjust = 1)) +
    theme(axis.ticks = element_line(linewidth = 0.3), legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-5, 0, 0, 0)) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
  save10(g, file.path(dir, "violin_inclusion_proportions.pdf"))
}

rp <- read.csv(file.path(S, "restricted_pairs.csv"), stringsAsFactors = FALSE)
rp <- rp[order(rp$n, rp$hstr, rp$baseline, rp$cens, rp$j, rp$variable), ]
types <- c("Strong Positive Effect", "Strong Negative Effect", "Small Positive Effect", "Small Negative Effect")
rp$type <- NA_character_
for (k in split(seq_len(nrow(rp)), paste(rp$n, rp$hstr, rp$baseline, rp$cens, rp$j))) rp$type[k[rp$active[k]]] <- types
for (s in list(c("LN", "0.25default", "EB-LNc75-sim2"), c("PGW", "0.5default", "EB-PGWc50-sim2"),
               c("PGW", "0.25default", "EB-PGWc25-sim2"))) for (h in c("PH", "AFT")) {
  z0 <- rp[rp$baseline == s[1] & rp$cens == s[2] & rp$hstr == h & !is.na(rp$pip.full.given.true), ]
  if (!nrow(z0)) next
  f <- file.path(fig_dir(s[3]), sprintf("comparison_%s.pdf", h))
  pdf(f, width = 15, height = 3); par(mfrow = c(1, 5), mar = c(4.5, 4.5, 2.5, 1))
  for (tp in c(types, "spurious")) {
    if (tp == "spurious") {
      z <- aggregate(cbind(pip.full.given.true, pip.restricted) ~ n + j, data = z0[!z0$active, ], FUN = mean)
      ttl <- sprintf("%s: Mean PIP for Spurious Variables", h); xl <- "Mean PIP (Conditioned)"; yl <- "Mean PIP (Reference)"
    } else { z <- z0[z0$active & z0$type == tp, ]; ttl <- sprintf("%s: %s", h, tp); xl <- "PIP (Conditioned)"; yl <- "PIP (Reference)" }
    plot(z$pip.full.given.true, z$pip.restricted, col = ifelse(z$n == 250, "blue", "red"), pch = 16, cex = 0.6,
         xlim = c(0, 1), ylim = c(0, 1), xlab = xl, ylab = yl, main = ttl)
    abline(0, 1, lty = 2, col = "grey")
    legend("bottomright", c("n=250", "n=500"), col = c("blue", "red"), pch = 16, cex = 0.9)
  }
  dev.off(); cat("  saved:", f, sprintf("(%d replicates without a visited model of the true structure left out)\n",
                                        length(unique(paste(rp$n, rp$j)[rp$baseline == s[1] & rp$cens == s[2] & rp$hstr == h & is.na(rp$pip.full.given.true)]))))
}
cat("\nDone: figures in", OUT, "\n")
