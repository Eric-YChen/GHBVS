suppressPackageStartupMessages({ library(ggplot2); library(grid) })

CELL_W <- 7.6745 / 10
CELL_H <- 6.1646 /  9

.is_null_unit <- function(u)
  vapply(seq_along(u), function(i) grepl("null", as.character(u[i])), logical(1))

.fit_title <- function(title, max_in) {
  tw <- function(s) convertWidth(
    grobWidth(textGrob(s, gp = gpar(fontsize = 24, fontface = "bold"))),
    "in", valueOnly = TRUE)
  if (tw(title) <= max_in) return(title)
  words <- strsplit(title, " ", fixed = TRUE)[[1]]
  if (length(words) < 2L) return(title)
  cand <- lapply(seq_len(length(words) - 1L), function(k) {
    l1 <- paste(words[seq_len(k)], collapse = " ")
    l2 <- paste(words[-seq_len(k)], collapse = " ")
    list(txt = paste0(l1, "\n", l2), wide = max(tw(l1), tw(l2)),
         comma = grepl(",$", l1))
  })
  fits <- Filter(function(z) z$wide <= max_in, cand)
  pool <- if (length(fits)) fits else cand
  at_comma <- Filter(function(z) z$comma, pool)
  if (length(at_comma)) pool <- at_comma
  pool[[which.min(vapply(pool, function(z) z$wide, numeric(1)))]]$txt
}

# p: ggplot object; path: output file; nrow, ncol: grid size in heatmap cells
save_sized <- function(p, path, nrow, ncol) {
  title <- p$labels$title
  pdf(NULL, width = 7, height = 7)
  ws <- ggplotGrob(p)$widths
  ow <- sum(convertWidth(ws[!.is_null_unit(ws)], "in", valueOnly = TRUE))
  w  <- ow + CELL_W * ncol
  wrapped <- FALSE
  if (is.character(title) && length(title) == 1L && nzchar(title)) {
    ttl <- .fit_title(title, w - 0.10)
    if (!identical(ttl, title)) { p <- p + labs(title = ttl); wrapped <- TRUE }
  }
  hs <- ggplotGrob(p)$heights
  oh <- sum(convertHeight(hs[!.is_null_unit(hs)], "in", valueOnly = TRUE))
  dev.off()
  h <- oh + CELL_H * nrow
  pdf(path, width = w, height = h); print(p); dev.off()
  cat(sprintf("  saved: %-52s %5.2f x %5.2f in  (%d x %d)%s\n",
              path, w, h, ncol, nrow, if (wrapped) "  [title wrapped]" else ""))
  invisible(p)
}

numformat <- function(val, digits = 3)
  sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", digits, "f"), val))

.axis_setup <- function(d) {
  hstr_order <- c("AH", "PH", "AFT", "GH")
  d$x_label <- paste0(d$hstr, "-", d$baseline, "\n(", d$c, "%)")
  d$y_label <- paste0("n=", d$n, "\np=", d$p)
  o <- order(match(d$hstr, hstr_order), d$baseline, d$c)
  x_levels <- unique(d$x_label[o])
  y_levels <- unique(d$y_label[order(d$p * 10000 + d$n)])
  d$x_label <- factor(d$x_label, levels = x_levels)
  d$y_label <- factor(d$y_label, levels = rev(y_levels))
  xg <- sub("-[A-Z]+\n.*$", "", x_levels)
  yg <- sub("^.*\np=", "", rev(y_levels))
  list(d  = d,
       xb = head(cumsum(rle(xg)$lengths), -1) + 0.5,
       yb = head(cumsum(rle(yg)$lengths), -1) + 0.5,
       nx = length(x_levels), ny = length(y_levels))
}

.base_theme <- function()
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "none", panel.grid = element_blank(),
        legend.title = element_blank())

create_heatmap <- function(data, metric_name, title, out_dir, filename,
                           col.values = c(0, 0.75, 1), reverse_scale = FALSE,
                           custom_limits = NULL, digits = 3, text_size = 6) {
  s  <- .axis_setup(data)
  hd <- data.frame(x_label = s$d$x_label, y_label = s$d$y_label,
                   value = s$d[[metric_name]])
  p <- ggplot(hd, aes(x = x_label, y = y_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    { if (length(s$yb)) geom_hline(yintercept = s$yb, color = "black", linewidth = 1.5) } +
    { if (length(s$xb)) geom_vline(xintercept = s$xb, color = "black", linewidth = 1.5) } +
    geom_text(aes(label = sapply(value, function(x) {
      if (is.na(x)) "NA"
      else if (is.null(custom_limits) && round(x, digits) == 1) "1"
      else numformat(x, digits) })), size = text_size, color = "black") +
    scale_fill_gradientn(
      colors = if (reverse_scale) c("green", "yellow", "red")
               else c("red", "yellow", "green"),
      limits = if (is.null(custom_limits)) c(0, 1) else custom_limits,
      values = col.values, name = metric_name) +
    labs(title = title, x = "", y = "") + .base_theme()
  save_sized(p, file.path(out_dir, filename), s$ny, s$nx)
}

# se_rule: "simple" (colour limit rounded up to 0.01) or "auto" (rounded up at the leading digit)
create_heatmap_se <- function(data, metric_name, title, out_dir, filename,
                              digits = 3, text_size = 6, se_rule = "auto") {
  s  <- .axis_setup(data)
  hd <- data.frame(x_label = s$d$x_label, y_label = s$d$y_label,
                   value = s$d[[metric_name]])
  mx <- max(hd$value, na.rm = TRUE)
  max_scale <- if (se_rule == "simple") {
    ceiling(mx * 100) / 100
  } else {
    step <- 10^floor(log10(mx))
    max(ceiling(mx / step) * step, ceiling(mx * 100) / 100 * (step >= 0.01))
  }
  p <- ggplot(hd, aes(x = x_label, y = y_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    { if (length(s$yb)) geom_hline(yintercept = s$yb, color = "black", linewidth = 1.5) } +
    { if (length(s$xb)) geom_vline(xintercept = s$xb, color = "black", linewidth = 1.5) } +
    geom_text(aes(label = sapply(value, numformat, digits = digits)),
              size = text_size, color = "black") +
    scale_fill_gradient2(low = "green", mid = "yellow", high = "red",
                         midpoint = max_scale / 2, limits = c(0, max_scale),
                         name = "SE") +
    labs(title = title, x = "", y = "") + .base_theme()
  save_sized(p, file.path(out_dir, filename), s$ny, s$nx)
}
