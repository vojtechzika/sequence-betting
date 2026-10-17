# ============================================================
# 03_r-distribution.R
#
# PURPOSE
#   Raincloud plots of the individual CRRA risk parameter r_i (posterior
#   mean of the Bayesian MPL model) by treatment: pooled over participants
#   and split by sex, each for all participants and for participants with a
#   consistent Holt-Laury choice pattern.
#
# INPUT
#   model "mpl": scored-<tr>.csv      posterior mean of r_i and consistency flag
#   path_src/participants.csv         sex
#
# OUTPUT (data/output/figures/)
#   r-distribution.png                   all participants, pooled over sex
#   r-distribution-sex.png               all participants, by sex
#   r-distribution-consistent.png        consistent participants, pooled over sex
#   r-distribution-consistent-sex.png    consistent participants, by sex
#
# NOTES
#   - One panel per treatment of cfg$run$treatment, coloured by the
#     treatment style (_reporting-style.R); in the figures by sex, men and
#     women share the colour of the treatment and the mean of each group is
#     marked by a segment in its own shade
#   - Participants whose sex is neither F nor M are shown in the pooled
#     figures but not in the figures by sex; the group sizes are printed in
#     the panels
#   - The axis range covers all values of r (at least [-0.5, 2], in steps of 0.5)
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
  
  # Axis range common to all panels and figures
  r_lo <- min(-0.5, floor(min(dt$r_mean) * 2) / 2)
  r_hi <- max(2,    ceiling(max(dt$r_mean) * 2) / 2)
  
  # ---- Panel: by sex, one treatment ----
  y_men <- 1; y_women <- 4
  
  panel_sex <- function(d, st) {
    d <- copy(d)[!is.na(sex_label)]
    lab <- paste0(st$type, " treatment (n = ", nrow(d), ")")
    d[, grp   := lab]
    d[, y_pos := fifelse(sex_label == "Men", y_men, y_women)]
    d[, seg_col := fifelse(sex_label == "Women", st$colour, st$colour_light)]
    
    group_means <- d[, .(r_mean_grp = mean(r_mean)), by = .(sex_label, y_pos, seg_col)]
    n_by_sex    <- d[, .N, by = sex_label]
    get_n <- function(sx) n_by_sex[sex_label == sx, N]
    
    ggplot(d, aes(x = r_mean, y = y_pos, group = sex_label)) +
      geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
      stat_halfeye(aes(fill = grp), adjust = 0.8, width = 0.6, .width = 0,
                   point_colour = NA, alpha = 0.80, position = position_nudge(y = 0.12)) +
      geom_jitter(aes(fill = grp), shape = 21, size = 1.8, alpha = 0.75, stroke = 0.25,
                  colour = "white", height = 0.07, width = 0) +
      geom_segment(data = group_means,
                   aes(x = r_mean_grp, xend = r_mean_grp, y = y_pos + 0.02, yend = y_pos + 0.72,
                       colour = seg_col), linewidth = 0.8) +
      scale_colour_identity() +
      geom_text(data = group_means,
                aes(x = r_mean_grp, y = y_pos - 0.25, label = round(r_mean_grp, 2)),
                size = 3, colour = "grey30", hjust = 0.5) +
      scale_fill_manual(values = setNames(st$colour, lab), name = NULL) +
      guides(fill = guide_legend(override.aes = list(shape = 21, colour = "white", size = 3))) +
      annotate("text", x = 0, y = 0.3, label = "Risk neutral", size = 2.6,
               colour = "grey45", hjust = 1.05, vjust = 0) +
      annotate("text", x = 1, y = 0.3, label = "Log utility", size = 2.6,
               colour = "grey45", hjust = -0.05, vjust = 0) +
      scale_x_continuous(name = expression(italic(r)[i] ~ "(posterior mean)"),
                         limits = c(r_lo, r_hi), breaks = seq(r_lo, r_hi, by = 0.5)) +
      scale_y_continuous(name = NULL, breaks = c(y_men, y_women),
                         labels = c(paste0("Men\nn = ", get_n("Men")),
                                    paste0("Women\nn = ", get_n("Women")))) +
      coord_cartesian(ylim = c(0.2, 7)) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = 11), axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 10),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
  }
  
  # ---- Panel: pooled over sex, one treatment ----
  panel_pooled <- function(d, st) {
    lab <- paste0(st$type, " treatment (n = ", nrow(d), ")")
    d <- copy(d)[, grp := lab]
    
    ggplot(d, aes(x = r_mean)) +
      geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
      stat_halfeye(aes(y = 1.1, fill = grp), adjust = 0.8, width = 1.3, .width = 0,
                   point_colour = NA, alpha = 0.80) +
      geom_jitter(aes(y = 0.85, fill = grp), shape = 21, size = 1.8, alpha = 0.75,
                  stroke = 0.25, colour = "white", height = 0.15, width = 0) +
      annotate("text", x = 0, y = 0.62, label = "Risk neutral", size = 2.6,
               colour = "grey45", hjust = 1.05, vjust = 0) +
      annotate("text", x = 1, y = 0.62, label = "Log utility", size = 2.6,
               colour = "grey45", hjust = -0.05, vjust = 0) +
      scale_fill_manual(values = setNames(st$colour, lab), name = NULL) +
      guides(fill = guide_legend(override.aes = list(shape = 21, colour = "white", size = 3))) +
      scale_x_continuous(name = expression(italic(r)[i] ~ "(posterior mean)"),
                         limits = c(r_lo, r_hi), breaks = seq(r_lo, r_hi, by = 0.5)) +
      scale_y_continuous(name = NULL, breaks = NULL) +
      coord_cartesian(ylim = c(0.55, 2.05)) +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = 11), axis.title.x = element_text(size = 10),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
  }
  
  build_and_save <- function(d, panel_fun, f_name, height) {
    set.seed(42)    # jitter
    panels <- lapply(seq_len(nrow(style)), function(k)
      panel_fun(d[treat == style$treat[k]], style[k]))
    fig <- Reduce(`+`, panels) + plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    f_png <- file.path(path_fig, f_name)
    ggsave(f_png, fig, width = 5 * nrow(style), height = height, dpi = 300)
    msg("Saved: ", f_png)
  }
  
  build_and_save(dt,                     panel_pooled, "r-distribution.png",                  2.8)
  build_and_save(dt,                     panel_sex,    "r-distribution-sex.png",              4.2)
  build_and_save(dt[inconsistent == 0L], panel_pooled, "r-distribution-consistent.png",       2.8)
  build_and_save(dt[inconsistent == 0L], panel_sex,    "r-distribution-consistent-sex.png",   4.2)
  
  invisible(TRUE)
}
