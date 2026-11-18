# ============================================================
# 04_r-distribution.R
#
# PURPOSE
#   Histograms of the individual CRRA risk parameter r_i (posterior mean of
#   the Bayesian MPL model) by treatment. Participants with a consistent
#   Holt-Laury choice pattern are shown pooled over sex and split by sex. A
#   further figure compares consistent and inconsistent participants.
#
# INPUT
#   model "mpl": scored-<tr>.csv      posterior mean of r_i and consistency flag
#   path_src/participants.csv         sex
#
# OUTPUT (data/output/figures/)
#   r-distribution-consistent.png        consistent participants, pooled over sex
#   r-distribution-consistent-sex.png    consistent participants, by sex
#   r-distribution-by-consistency.png    consistent and inconsistent participants,
#                                        pooled over sex
#
# NOTES
#   - One panel per treatment of cfg$run$treatment, coloured by the
#     treatment style (_reporting-style.R); where a figure compares two
#     groups, the first (women, consistent) is drawn in the colour of the
#     treatment and the second (men, inconsistent) in its light shade
#   - Bars show the share of the participants of the group in bins of width
#     0.25, so that groups of different size can be compared; the bins are
#     the same in all figures. The mean of each group is marked by a vertical
#     line and printed next to it
#   - Participants whose sex is neither F nor M are shown in the pooled
#     figures but not in the figure by sex; the group sizes are printed in
#     the panels
#   - The vertical axes of the two rows of a panel share the scale, except in
#     the figure by consistency, where the rows have their own scale because
#     the inconsistent estimates are concentrated in a few bins
#   - The axis range covers all values of r (at least [-0.5, 2], in steps of 0.5)
#   - If the bins with participants are separated by a stretch of empty bins
#     of at least 0.75 (the extreme low estimates), the horizontal axis is
#     broken there: the empty stretch is compressed to a short gap marked by
#     "//", the tick labels remain in units of r, and distances are
#     preserved on both sides of the break
# ============================================================

r_distribution <- function(cfg) {
  
  style  <- treatment_style(cfg)
  tr_vec <- style$treat
  
  # ---- Posterior means of r and consistency flag ----
  dt <- rbindlist(lapply(tr_vec, function(tr) {
    sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
    sc[, .(pid = as.character(pid), treat = tr, r_mean, inconsistent = as.integer(inconsistent))]
  }))
  
  participants <- fread(file.path(path_src, "participants.csv"), encoding = "UTF-8")
  participants[, pid := as.character(pid)]
  dt <- merge(dt, participants[, .(pid, sex = toupper(trimws(sex)))], by = "pid", all.x = TRUE)
  dt[, sex_label := fcase(sex == "F", "Women", sex == "M", "Men", default = NA_character_)]
  
  # Axis range and bins common to all panels and figures
  r_lo <- min(-0.5, floor(min(dt$r_mean) * 2) / 2)
  r_hi <- max(2,    ceiling(max(dt$r_mean) * 2) / 2)
  bw   <- 0.25
  
  # ---- Axis break over the largest empty stretch of bins ----
  # Occupied bins are located from the posterior means of all participants, so
  # that the break is the same in every figure
  occ  <- sort(unique(floor((dt$r_mean - r_lo) / bw) * bw + r_lo))
  gaps <- diff(occ) - bw
  i_g  <- which.max(gaps)
  has_break <- length(gaps) > 0 && gaps[i_g] >= 0.75
  g_len <- 0.2                                   # width of the compressed stretch
  if (has_break) {
    cut_lo <- occ[i_g] + bw + bw / 2             # empty stretch that is removed
    cut_hi <- occ[i_g + 1] - bw / 2
    shift  <- (cut_hi - cut_lo) - g_len
  }
  # position on the broken axis (identity above the break, shifted below it)
  f <- function(x) {
    if (!has_break) return(x)
    ifelse(x <= cut_lo, x + shift,
           ifelse(x >= cut_hi, x, cut_hi - g_len + (x - cut_lo) * g_len / (cut_hi - cut_lo)))
  }
  brk <- seq(r_lo, r_hi, by = 0.5)
  brk <- brk[brk >= min(dt$r_mean) - bw]         # no tick below the data
  if (has_break) brk <- brk[brk <= cut_lo | brk >= cut_hi]
  
  # ---- Panel: histograms of one treatment, one row per group ----
  # d has the columns r_mean and grp; grp_order lists the groups from top to
  # bottom and fills their colours. With one group the row label is dropped
  # and the group size is printed in the title.
  hist_panel <- function(d, st, grp_order, fills, free_y = FALSE) {
    d <- copy(d)[grp %in% grp_order]
    n_grp <- d[, .N, by = grp]
    lab   <- paste0(grp_order, "\nn = ", n_grp[match(grp_order, grp), N])
    d[, grp_lab := factor(lab[match(grp, grp_order)], levels = lab)]
    
    # shares of the group in each bin (bins start at r_lo)
    d[, bin := floor((r_mean - r_lo) / bw) * bw + r_lo + bw / 2]
    h <- d[, .N, by = .(grp_lab, bin)]
    h[, share := N / sum(N), by = grp_lab]
    h[, bin_x := f(bin)]
    
    gm   <- d[, .(m = mean(r_mean)), by = grp_lab]
    gm[, m_x := f(m)]
    refs <- data.table(grp_lab = factor(lab[1], levels = lab))   # labels in the top row
    
    one   <- length(grp_order) == 1L
    # baseline of the x axis, drawn by hand so that it can be interrupted at the break
    x_lim <- f(c(r_lo, r_hi))
    if (has_break) {
      x_mid <- cut_hi - g_len / 2
      base  <- data.table(x0 = c(x_lim[1], x_mid + 0.08), x1 = c(x_mid - 0.08, x_lim[2]))
    } else {
      base  <- data.table(x0 = x_lim[1], x1 = x_lim[2])
    }
    
    title <- if (one) paste0(st$type, " treatment (n = ", nrow(d), ")") else paste0(st$type, " treatment")
    
    g <- ggplot(h, aes(x = bin_x, y = share, fill = grp_lab)) +
      geom_col(width = bw, colour = "white", linewidth = 0.25, alpha = 0.9) +
      geom_vline(xintercept = f(c(0, 1)), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
      geom_vline(data = gm, aes(xintercept = m_x), colour = "grey20", linewidth = 0.7) +
      geom_text(data = gm, aes(x = m_x, y = Inf, label = sprintf("mean %.2f", m)), inherit.aes = FALSE,
                size = 3, colour = "grey20", hjust = -0.08, vjust = 3.0) +
      geom_text(data = refs, aes(x = f(0), y = Inf), label = "Risk neutral", inherit.aes = FALSE,
                size = 2.6, colour = "grey45", hjust = 1.05, vjust = 1.6) +
      geom_text(data = refs, aes(x = f(1), y = Inf), label = "Log utility", inherit.aes = FALSE,
                size = 2.6, colour = "grey45", hjust = -0.05, vjust = 1.6) +
      scale_fill_manual(values = setNames(fills, lab), guide = "none") +
      geom_segment(data = base, aes(x = x0, xend = x1, y = 0, yend = 0), inherit.aes = FALSE,
                   colour = "black", linewidth = 0.5) +
      scale_x_continuous(name = expression(italic(r)[i] ~ "(posterior mean)"),
                         limits = x_lim, breaks = f(brk), labels = sprintf("%.1f", brk)) +
      scale_y_continuous(name = "Share of participants", expand = expansion(mult = c(0, 0.30)),
                         labels = function(x) paste0(round(100 * x), "%")) +
      facet_grid(grp_lab ~ ., switch = "y", scales = if (free_y) "free_y" else "fixed") +
      coord_cartesian(clip = "off") +
      labs(title = title) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(size = 11, face = "bold"),
            strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = if (one) element_blank() else element_text(angle = 0, hjust = 1, size = 10),
            axis.title = element_text(size = 10), axis.line.x = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
    
    # break marker on the axis
    if (has_break) {
      g <- g + geom_text(data = data.table(x = x_mid, y = 0), aes(x = x, y = y), label = "//",
                         inherit.aes = FALSE, size = 4, fontface = "bold")
    }
    g
  }
  
  panel_pooled <- function(d, st) {
    d <- copy(d)[, grp := "All"]
    hist_panel(d, st, "All", st$colour)
  }
  
  panel_sex <- function(d, st) {
    d <- copy(d)[!is.na(sex_label)][, grp := sex_label]
    hist_panel(d, st, c("Women", "Men"), c(st$colour, st$colour_light))
  }
  
  panel_consistency <- function(d, st) {
    d <- copy(d)[, grp := fifelse(inconsistent == 1L, "Inconsistent", "Consistent")]
    hist_panel(d, st, c("Consistent", "Inconsistent"), c(st$colour, st$colour_light), free_y = TRUE)
  }
  
  build_and_save <- function(d, panel_fun, f_name, height) {
    panels <- lapply(seq_len(nrow(style)), function(k)
      panel_fun(d[treat == style$treat[k]], style[k]))
    fig <- Reduce(`+`, panels)
    f_png <- file.path(path_fig, f_name)
    ggsave(f_png, fig, width = 5 * nrow(style), height = height, dpi = 300)
    msg("Saved: ", f_png)
  }
  
  build_and_save(dt[inconsistent == 0L], panel_pooled,      "r-distribution-consistent.png",       2.8)
  build_and_save(dt[inconsistent == 0L], panel_sex,         "r-distribution-consistent-sex.png",   4.2)
  build_and_save(dt,                     panel_consistency, "r-distribution-by-consistency.png",   4.2)
  
  invisible(TRUE)
}
