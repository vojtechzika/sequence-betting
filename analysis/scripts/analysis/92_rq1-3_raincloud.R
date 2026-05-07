# ============================================================
# 92_rq1_3_raincloud.R
# ============================================================

rq1_3_raincloud <- function(cfg) {
  
  library(data.table)
  library(ggplot2)
  library(ggdist)
  library(patchwork)
  library(scales)
  
  tr  <- "m25"
  tag <- "confirmatory"
  
  f_out <- file.path(path_fig, "rq1_3_participants_raincloud.png")
  if (should_skip(f_out, cfg, "output", "RQ1-3 combined raincloud")) return(invisible(NULL))
  
  # ---- Helpers -----------------------------------------------------------
  load_pid <- function(rq, tg) {
    f <- file.path(path_out, paste0(rq, "_", tr, "_", tg, "_participants.csv"))
    if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
    fread(f, encoding = "UTF-8")
  }
  
  load_grand <- function(rq, param_label, tg = tag) {
    f <- file.path(path_out, paste0(rq, "_", tr, "_", tg, "_model_summary.csv"))
    if (!file.exists(f)) return(list(median = NA_real_, lo = NA_real_, hi = NA_real_))
    mod <- fread(f, encoding = "UTF-8")
    row <- mod[parameter == param_label]
    if (nrow(row) == 0) return(list(median = NA_real_, lo = NA_real_, hi = NA_real_))
    list(median = row$median, lo = row$q025, hi = row$q975)
  }
  
  # ---- Load data ---------------------------------------------------------
  rq1_conf <- load_pid("rq1", "confirmatory")
  rq1_full <- load_pid("rq1", "full")
  rq2_conf <- load_pid("rq2", "confirmatory")
  rq2_full <- load_pid("rq2", "full")
  rq3_conf <- load_pid("rq3", "confirmatory")
  rq3_full <- load_pid("rq3", "full")
  
  gm <- list(
    rq1      = load_grand("rq1", "grand mean bet probability"),
    rq2      = load_grand("rq2", "grand mean mu_a (absolute scale, proportion of endowment)"),
    rq3      = load_grand("rq3", "grand mean mu_c (proportion of endowment)"),
    rq1_full = load_grand("rq1", "grand mean bet probability",                                "full"),
    rq2_full = load_grand("rq2", "grand mean mu_a (absolute scale, proportion of endowment)", "full"),
    rq3_full = load_grand("rq3", "grand mean mu_c (proportion of endowment)",                 "full")
  )
  
  # ---- Shared colours ----------------------------------------------------
  col_conf  <- "#74ADD1"
  col_point <- "#2166AC"
  col_full  <- "grey60"
  col_gm    <- "steelblue"
  col_eu    <- "grey30"
  
  # ---- Panel builder -----------------------------------------------------
  make_panel <- function(dt_conf, dt_full, x_var,
                         gm_conf, gm_full,
                         eu_ref, eu_hjust,
                         x_label, x_limits = NULL, x_breaks = NULL,
                         jitter_fill = col_point,
                         jitter_aes  = NULL,
                         label_full_hjust =  0.1, label_full_vjust = 3.5,
                         label_conf_hjust = -0.1, label_conf_vjust = 1.5) {
    
    p <- ggplot(dt_conf, aes(x = .data[[x_var]], y = 1))
    
    # Ghost
    if (!is.null(dt_full))
      p <- p + stat_halfeye(
        data = dt_full, aes(x = .data[[x_var]], y = 1),
        fill = col_full, adjust = 0.8, width = 0.5, .width = 0,
        point_colour = NA, alpha = 0.40,
        position = position_nudge(y = 0.12), inherit.aes = FALSE
      )
    
    # Reference lines
    p <- p +
      geom_vline(xintercept = gm_full$median, colour = col_full, linewidth = 0.5, linetype = "solid") +
      geom_vline(xintercept = gm_conf$median, colour = col_gm,   linewidth = 0.6, linetype = "solid")
    
    # Annotations
    p <- p +
      annotate("text", x = eu_ref,         y = Inf, label = "",
               hjust = eu_hjust, vjust = 1.5, size = 2.8, colour = col_eu) +
      annotate("text", x = gm_full$median, y = Inf,
               label = sprintf("%.1f%%", 100 * gm_full$median),
               hjust = label_full_hjust, vjust = label_full_vjust,
               size = 2.8, colour = col_full) +
      annotate("text", x = gm_conf$median, y = Inf,
               label = sprintf("%.1f%%", 100 * gm_conf$median),
               hjust = label_conf_hjust, vjust = label_conf_vjust,
               size = 2.8, colour = col_gm)
    
    # Confirmatory distribution
    p <- p +
      stat_halfeye(
        fill = col_conf, adjust = 0.8, width = 0.5, .width = 0,
        point_colour = NA, alpha = 0.80, position = position_nudge(y = 0.12)
      )
    
    # Jitter
    if (is.null(jitter_aes)) {
      p <- p + geom_jitter(
        fill = col_point, shape = 21, size = 1.8, alpha = 0.75,
        stroke = 0.25, colour = "white", height = 0.07, width = 0
      )
    } else {
      p <- p + geom_jitter(
        jitter_aes,
        shape = 21, size = 1.8, alpha = 0.75,
        stroke = 0.25, colour = "white", height = 0.07, width = 0
      )
    }
    
    # Scales
    p <- p + if (!is.null(x_breaks)) {
      scale_x_continuous(
        name   = x_label,
        limits = x_limits,
        breaks = x_breaks,
        labels = scales::percent_format(accuracy = 1),
        expand = c(0, 0)
      )
    } else {
      scale_x_continuous(
        name   = x_label,
        limits = x_limits,
        labels = scales::percent_format(accuracy = 1),
        expand = c(0, 0)
      )
    }
    
    p + scale_y_continuous(name = NULL, breaks = NULL) +
      theme_classic(base_size = 11) +
      theme(
        axis.title.x       = element_text(size = 10),
        axis.text.y        = element_blank(),
        axis.line.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        legend.position    = "none"
      )
  }
  
  # ---- Build panels ------------------------------------------------------
  # Each panel has its own label_full_hjust, label_full_vjust,
  #                          label_conf_hjust, label_conf_vjust
  p1 <- make_panel(
    dt_conf  = rq1_conf, dt_full = rq1_full, x_var = "mu_i_median",
    gm_conf  = gm$rq1,   gm_full = gm$rq1_full,
    eu_ref   = 1.0,      eu_hjust = -0.1,
    x_label  = expression(hat(mu)[i]^b ~ "(posterior median betting probability)"),
    x_limits = c(-0.05, 1.05),
    x_breaks = seq(0, 1, by = 0.1),
    label_full_hjust =  -0.1, label_full_vjust = 1,
    label_conf_hjust = 1.1, label_conf_vjust = 1
  )
  
  p2 <- make_panel(
    dt_conf  = rq2_conf, dt_full = rq2_full, x_var = "mu_a_median",
    gm_conf  = gm$rq2,   gm_full = gm$rq2_full,
    eu_ref   = 0,        eu_hjust = 1.15,
    x_label  = expression(hat(mu)[i]^a ~ "(posterior median proportional stake deviation)"),
    label_full_hjust =  -0.1, label_full_vjust = 1,
    label_conf_hjust = 1.1, label_conf_vjust = 1
  )
  
  p3 <- make_panel(
    dt_conf    = rq3_conf, dt_full = rq3_full, x_var = "mu_c_median",
    gm_conf    = gm$rq3,   gm_full = gm$rq3_full,
    eu_ref     = 0,        eu_hjust = -0.15,
    x_label    = expression(hat(mu)[i]^c ~ "(posterior median welfare loss, share of endowment)"),
    jitter_aes = aes(fill = mu_b_median),
    label_full_hjust =  -0.1, label_full_vjust = 1,
    label_conf_hjust = 1.1, label_conf_vjust = 1
  ) +
    scale_fill_gradient(
      low    = "#E31A1C",
      high   = "#2166AC",
      name   = "Mean betting probability",
      labels = scales::percent_format(accuracy = 1)
    )
  
  # ---- Combine -----------------------------------------------------------
  combined <- (p1 + theme(plot.margin = margin(8, 8, 8, 0))) +
    (p2 + theme(plot.margin = margin(8, 0, 8, 8))) +
    (p3 + theme(plot.margin = margin(8, 8, 8, 8))) +
    plot_layout(ncol = 2, nrow = 2, heights = c(1, 1),
                design = "AB\nCC") +
    plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
  
  ggsave(f_out, combined, width = 8, height = 5, dpi = 600)
  msg("Saved: ", f_out)
  
  invisible(TRUE)
}