# ============================================================
# 03_r_raincloud.R
#
# PURPOSE
#   Raincloud plot of individual CRRA risk parameters (r_i)
#   by sex and treatment. Produces full-sample and
#   consistent-only versions.
#
# INPUT
#   path_mod/mpl_r_draws_m25.rds  -- posterior draws FN
#   path_mod/mpl_r_draws_m19.rds  -- posterior draws FP
#   path_src/participants.csv     -- for sex variable
#   path_out/mpl_scored_m25.csv   -- for inconsistency flag
#   path_out/mpl_scored_m19.csv   -- for inconsistency flag
#
# OUTPUT
#   path_fig/r_distribution.png            -- full sample
#   path_fig/r_distribution_consistent.png -- consistent only
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
  
  dt[, fill_col := fcase(
    treat == "FN (m = 2.5)" & sex_label == "Women", col_fn_w,
    treat == "FN (m = 2.5)" & sex_label == "Men",   col_fn_m,
    treat == "FP (m = 1.9)" & sex_label == "Women", col_fp_w,
    treat == "FP (m = 1.9)" & sex_label == "Men",   col_fp_m
  )]
  
  # ---- Panel builder ----
  make_panel <- function(data_sub, title) {
    
    n_counts    <- data_sub[, .N, by = sex_label]
    group_means <- data_sub[, .(r_mean_grp = mean(r_mean, na.rm = TRUE)),
                            by = .(sex_label, y_pos, fill_col)]
    
    get_n <- function(sx) n_counts[sex_label == sx, N]
    
    y_breaks <- c(y_men, y_women)
    y_labels <- c(paste0("Men\nn = ",   get_n("Men")),
                  paste0("Women\nn = ", get_n("Women")))
    
    ggplot(data_sub, aes(x = r_mean, y = y_pos)) +
      
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      
      stat_halfeye(
        aes(fill = fill_col),
        adjust = 0.8, width = 0.6, .width = 0,
        point_colour = NA, alpha = 0.80,
        position = position_nudge(y = 0.12)
      ) +
      
      geom_jitter(
        aes(fill = fill_col),
        shape = 21, size = 1.8, alpha = 0.75,
        stroke = 0.25, colour = "white",
        height = 0.07, width = 0, seed = 42
      ) +
      
      geom_segment(
        data = group_means,
        aes(x = r_mean_grp, xend = r_mean_grp,
            y = y_pos + 0.02, yend = y_pos + 0.72,
            colour = fill_col),
        linewidth = 0.8
      ) +
      scale_colour_identity() +
      
      geom_text(
        data = group_means,
        aes(x = r_mean_grp, y = y_pos - 0.25,
            label = round(r_mean_grp, 2)),
        size = 3, colour = "grey30", hjust = 0.5
      ) +
      
      scale_fill_identity() +
      
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
      coord_cartesian(ylim = c(0.4, 7)) +
      labs(title = title) +
      
      theme_classic(base_size = 11) +
      theme(
        plot.title         = element_text(face = "bold", size = 11),
        axis.text.y        = element_text(size = 10),
        axis.title.x       = element_text(size = 10),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        legend.position    = "none"
      )
  }
  
  # ---- Build and save ----
  make_and_save <- function(data_sub, suffix) {
    
    set.seed(42)
    fig <- make_panel(data_sub[treat == "FN (m = 2.5)"],
                      "FN treatment (m = 2.5)") +
      make_panel(data_sub[treat == "FP (m = 1.9)"],
                 "FP treatment (m = 1.9)") +
      plot_annotation(
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      )
    
    f_png <- file.path(path_fig, paste0("r_distribution", suffix, ".png"))
    ggsave(f_png, fig, width = 10, height = 4, dpi = 300)
    msg("Saved: ", f_png)
    invisible(fig)
  }
  
  make_and_save(dt, "")
  if (isTRUE(cfg$run$consistent_only)) {
    make_and_save(dt[inconsistent == 0L], "_consistent")
  }
  
  invisible(TRUE)
}