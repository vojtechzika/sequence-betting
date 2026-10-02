# ============================================================
# 03_r_raincloud.R
#
# PURPOSE
#   Raincloud plot of individual CRRA risk parameters (r_i)
#   by sex and treatment. Produces full-sample and
#   consistent-only versions, each both by sex and pooled
#   (no sex facet), by treatment.
#
#   Per-panel treatment titles are rendered as a bottom,
#   fill-mapped legend (colored swatches) rather than a plot
#   title, for both the by-sex and the pooled figures. Legends
#   from the two side-by-side panels are collected into one
#   shared legend via patchwork::plot_layout(guides = "collect").
#   Within a by-sex panel, men and women share the panel's
#   treatment color; sex is distinguished by the y-axis
#   position/label and by the color of each group's posterior-
#   mean segment (still an identity scale, independent of the
#   legend's fill scale).
#
# INPUT
#   path_mod/mpl_r_draws_m25.rds  -- posterior draws FN
#   path_mod/mpl_r_draws_m19.rds  -- posterior draws FP
#   path_src/participants.csv     -- for sex variable
#   path_out/mpl_scored_m25.csv   -- for inconsistency flag
#   path_out/mpl_scored_m19.csv   -- for inconsistency flag
#
# OUTPUT
#   path_fig/r-distribution-sex.png             -- full sample, by sex
#   path_fig/r-distribution-pooled.png          -- full sample, pooled (no sex facet)
#   path_fig/r-distribution-consistent-sex.png  -- consistent only, by sex
#   path_fig/r-distribution-consistent.png      -- consistent only, pooled (no sex facet)
# ============================================================

r_raincloud <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  f_draws_m25    <- file.path(path_mod, "mpl_r_draws_m25.rds")
  f_draws_m19    <- file.path(path_mod, "mpl_r_draws_m19.rds")
  f_participants <- file.path(path_src, "participants.csv")
  
  stopifnot(file.exists(f_draws_m25), file.exists(f_draws_m19),
            file.exists(f_participants))
  
  # ---- Load draws ----
  summarise_draws <- function(path, treat_label) {
    obj     <- readRDS(path)
    pids    <- as.character(obj$pid)
    r_draws <- obj$r_draws
    stopifnot(is.matrix(r_draws), ncol(r_draws) == length(pids))
    data.table(pid = pids, treat = treat_label,
               r_mean = apply(r_draws, 2, mean))
  }
  
  dt <- rbind(
    summarise_draws(f_draws_m25, "FN (m = 2.5)"),
    summarise_draws(f_draws_m19, "FP (m = 1.9)")
  )
  
  # ---- Attach sex ----
  participants <- fread(f_participants, encoding = "UTF-8")
  participants[, pid := as.character(pid)]
  dt <- merge(dt, participants[, .(pid, sex = toupper(trimws(sex)))],
              by = "pid", all.x = TRUE)
  dt[, sex_label := fcase(sex == "F", "Women", sex == "M", "Men",
                          default = NA_character_)]
  dt <- dt[!is.na(sex_label)]
  dt[, treat := factor(treat, levels = c("FN (m = 2.5)", "FP (m = 1.9)"))]
  
  # ---- Attach inconsistency flag ----
  f_scored_m25 <- file.path(path_out, "mpl_scored_m25.csv")
  f_scored_m19 <- file.path(path_out, "mpl_scored_m19.csv")
  stopifnot(file.exists(f_scored_m25), file.exists(f_scored_m19))
  
  scored <- rbind(
    fread(f_scored_m25)[, .(pid = as.character(pid), inconsistent)],
    fread(f_scored_m19)[, .(pid = as.character(pid), inconsistent)]
  )
  dt <- merge(dt, scored, by = "pid", all.x = TRUE)
  
  # ---- Layout ----
  y_men   <- 1
  y_women <- 4
  
  dt[, y_pos := fcase(sex_label == "Men",   y_men,
                      sex_label == "Women", y_women)]
  
  # ---- Colours ----
  col_fn_w <- "#2166AC"; col_fn_m <- "#74ADD1"
  col_fp_w <- "#C0392B"; col_fp_m <- "#E8927C"
  
  dt[, seg_col := fcase(
    treat == "FN (m = 2.5)" & sex_label == "Women", col_fn_w,
    treat == "FN (m = 2.5)" & sex_label == "Men",   col_fn_m,
    treat == "FP (m = 1.9)" & sex_label == "Women", col_fp_w,
    treat == "FP (m = 1.9)" & sex_label == "Men",   col_fp_m
  )]
  
  # ---- Panel builder (by sex, within one treatment) ----
  make_panel <- function(data_sub, title, rep_col) {
    
    data_sub <- copy(data_sub)
    lab <- paste0(title, " (n = ", nrow(data_sub), ")")
    data_sub[, grp := lab]
    
    n_counts    <- data_sub[, .N, by = sex_label]
    group_means <- data_sub[, .(r_mean_grp = mean(r_mean, na.rm = TRUE)),
                            by = .(sex_label, y_pos, seg_col)]
    
    get_n <- function(sx) n_counts[sex_label == sx, N]
    
    y_breaks <- c(y_men, y_women)
    y_labels <- c(paste0("Men\nn = ",   get_n("Men")),
                  paste0("Women\nn = ", get_n("Women")))
    
    ggplot(data_sub, aes(x = r_mean, y = y_pos, group = sex_label)) +
      
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      
      stat_halfeye(
        aes(fill = grp),
        adjust = 0.8, width = 0.6, .width = 0,
        point_colour = NA, alpha = 0.80,
        position = position_nudge(y = 0.12)
      ) +
      
      geom_jitter(
        aes(fill = grp),
        shape = 21, size = 1.8, alpha = 0.75,
        stroke = 0.25, colour = "white",
        height = 0.07, width = 0, seed = 42
      ) +
      
      geom_segment(
        data = group_means,
        aes(x = r_mean_grp, xend = r_mean_grp,
            y = y_pos + 0.02, yend = y_pos + 0.72,
            colour = seg_col),
        linewidth = 0.8
      ) +
      scale_colour_identity() +
      
      geom_text(
        data = group_means,
        aes(x = r_mean_grp, y = y_pos - 0.25,
            label = round(r_mean_grp, 2)),
        size = 3, colour = "grey30", hjust = 0.5
      ) +
      
      scale_fill_manual(values = setNames(rep_col, lab), name = NULL) +
      guides(fill = guide_legend(
        override.aes = list(shape = 21, colour = "white", size = 3)
      )) +
      
      annotate("text", x = 0, y = 0.3, label = "Risk neutral",
               size = 2.6, colour = "grey45", hjust = 1.05, vjust = 0) +
      annotate("text", x = 1, y = 0.3, label = "Log utility",
               size = 2.6, colour = "grey45", hjust = -0.05, vjust = 0) +
      
      scale_x_continuous(
        name   = expression(italic(r)[i] ~ "(posterior mean)"),
        limits = c(-0.5, 2),
        breaks = seq(-0.5, 2, by = 0.5)
      ) +
      scale_y_continuous(name = NULL, breaks = y_breaks, labels = y_labels) +
      coord_cartesian(ylim = c(0.2, 7)) +
      
      theme_classic(base_size = 11) +
      theme(
        legend.position     = "bottom",
        legend.title        = element_blank(),
        legend.text         = element_text(size = 11),
        axis.text.y         = element_text(size = 10),
        axis.title.x        = element_text(size = 10),
        panel.grid.major.x  = element_line(colour = "grey92", linewidth = 0.3)
      )
  }
  
  # ---- Panel builder (pooled, no sex facet, within one treatment) ----
  make_panel_pooled <- function(data_sub, title, fill_col) {
    
    data_sub <- copy(data_sub)
    lab <- paste0(title, " (n = ", nrow(data_sub), ")")
    data_sub[, grp := lab]
    
    ggplot(data_sub, aes(x = r_mean)) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      stat_halfeye(
        aes(y = 1.1, fill = grp),
        adjust = 0.8, width = 1.3, .width = 0,
        point_colour = NA, alpha = 0.80
      ) +
      geom_jitter(
        aes(y = 0.85, fill = grp),
        shape = 21, size = 1.8, alpha = 0.75,
        stroke = 0.25, colour = "white",
        height = 0.15, width = 0, seed = 42
      ) +
      annotate("text", x = 0, y = 0.62, label = "Risk neutral",
               size = 2.6, colour = "grey45", hjust = 1.05, vjust = 0) +
      annotate("text", x = 1, y = 0.62, label = "Log utility",
               size = 2.6, colour = "grey45", hjust = -0.05, vjust = 0) +
      scale_fill_manual(values = setNames(fill_col, lab), name = NULL) +
      guides(fill = guide_legend(
        override.aes = list(shape = 21, colour = "white", size = 3)
      )) +
      scale_x_continuous(
        name   = expression(italic(r)[i] ~ "(posterior mean)"),
        limits = c(-0.5, 2),
        breaks = seq(-0.5, 2, by = 0.5)
      ) +
      scale_y_continuous(name = NULL, breaks = NULL) +
      coord_cartesian(ylim = c(0.55, 2.05)) +
      theme_classic(base_size = 11) +
      theme(
        legend.position     = "bottom",
        legend.title        = element_blank(),
        legend.text         = element_text(size = 11),
        axis.title.x        = element_text(size = 10),
        panel.grid.major.x  = element_line(colour = "grey92", linewidth = 0.3)
      )
  }
  
  # ---- Build and save (by sex) ----
  make_and_save <- function(data_sub, f_name) {
    
    set.seed(42)
    
    fig <- (make_panel(data_sub[treat == "FN (m = 2.5)"], "FN treatment", col_fn_w) +
              make_panel(data_sub[treat == "FP (m = 1.9)"], "FP treatment", col_fp_w)) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    f_png <- file.path(path_fig, f_name)
    ggsave(f_png, fig, width = 10, height = 4.2, dpi = 300)
    msg("Saved: ", f_png)
    invisible(fig)
  }
  
  # ---- Build and save (pooled, no sex facet, consistent sample only) ----
  make_and_save_pooled <- function(data_sub, f_name) {
    
    fig <- (make_panel_pooled(data_sub[treat == "FN (m = 2.5)"],
                              "FN treatment", col_fn_w) +
              make_panel_pooled(data_sub[treat == "FP (m = 1.9)"],
                                "FP treatment", col_fp_w)) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    f_png <- file.path(path_fig, f_name)
    ggsave(f_png, fig, width = 10, height = 2.8, dpi = 300)
    msg("Saved: ", f_png)
    invisible(fig)
  }
  
  make_and_save(dt, "r-distribution-sex.png")
  make_and_save_pooled(dt, "r-distribution-pooled.png")
  if (isTRUE(cfg$run$consistent_only)) {
    make_and_save(dt[inconsistent == 0L], "r-distribution-consistent-sex.png")
    make_and_save_pooled(dt[inconsistent == 0L], "r-distribution-consistent.png")
  }
  
  invisible(TRUE)
}