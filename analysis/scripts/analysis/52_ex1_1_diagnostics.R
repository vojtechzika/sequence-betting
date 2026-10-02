# ============================================================
# 52_ex1_1_diagnostics.R
#
# Diagnostics for EX1.I.1 anchor similarity outputs.
#
# INPUTS
#   path_out/ex1_1_anchors_<tr>.csv   (from 51_)
#
# OUTPUTS
#   path_out/ex1_1_anchor_check_<tr>.csv
#   path_fig/ex1_1_<tr>_wH_main.png
#   path_fig/ex1_1_<tr>_wT_main.png
#   path_fig/ex1_1_<tr>_direction_main.png
#   path_fig/ex1_1_<tr>_scatter_wH_wT.png        (per-treatment)
#   path_fig/ex1_1_scatter_panel.pdf              (two-treatment panel, journal)
#   path_fig/ex1_1_<tr>_eps_sensitivity.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

ex1_1_diagnostics <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  lab_H  <- as.character(cfg$design$seq$anchor_labels$pure_heads)   # "HHHHHH"
  lab_T  <- as.character(cfg$design$seq$anchor_labels$pure_tails)   # "TTTTTT"
  
  # Human-readable treatment names used in figure titles / strips
  tr_labels <- c(m25 = "FN treatment (m = 2.5)", m19 = "FP treatment (m = 1.9)")
  
  # ---- Shared journal theme ---------------------------------------------------
  theme_journal <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
      theme(
        # Axes
        axis.line        = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks       = element_line(colour = "grey30", linewidth = 0.35),
        axis.text        = element_text(colour = "grey20", size = base_size - 1),
        axis.title       = element_text(colour = "grey10", size = base_size,
                                        margin = margin(t = 4, r = 4)),
        # Panel
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        # Legend
        legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.title     = element_text(size = base_size - 1, colour = "grey20"),
        legend.text      = element_text(size = base_size - 1, colour = "grey20"),
        legend.key.size  = unit(0.45, "cm"),
        legend.margin    = margin(t = 2),
        legend.box.margin = margin(0, 0, 0, 0),
        # Strip (facet label)
        strip.background = element_rect(fill = "grey96", colour = "grey70",
                                        linewidth = 0.4),
        strip.text       = element_text(size = base_size, colour = "grey10",
                                        face = "bold", margin = margin(4, 4, 4, 4)),
        # Title / caption
        plot.title       = element_text(size = base_size + 1, face = "bold",
                                        colour = "grey10", margin = margin(b = 4)),
        plot.subtitle    = element_text(size = base_size - 1, colour = "grey40",
                                        margin = margin(b = 6)),
        plot.caption     = element_text(size = base_size - 2, colour = "grey50",
                                        hjust = 0, margin = margin(t = 6)),
        plot.margin      = margin(8, 10, 6, 8)
      )
  }
  
  # ---- Colour / shape palette -------------------------------------------------
  # Three point categories: anchor (HHHHHH / TTTTTT), extreme (top-5), other
  pal_fill <- c(
    "Anchor"   = "#1a1a2e",   # near-black navy
    "Extreme"  = "#c0392b",   # dark red
    "Other"    = "#aab4c4"    # muted steel-blue
  )
  pal_colour <- c(
    "Anchor"   = "#1a1a2e",
    "Extreme"  = "#9b2819",
    "Other"    = "#8494a8"
  )
  shape_map <- c(Anchor = 23, Extreme = 21, Other = 21)  # diamond vs circle
  size_map  <- c(Anchor = 3.2, Extreme = 2.4, Other = 1.6)
  
  # ---- Helper: save PNG -------------------------------------------------------
  save_plot <- function(p, fn, w = 10, h = 7) {
    ggsave(file.path(path_fig, fn), plot = p, width = w, height = h,
           units = "in", dpi = 300)
    msg("Saved: ", file.path(path_fig, fn))
  }
  
  # ---- Helper: label set for one treatment ------------------------------------
  # Returns the sequences to label: anchors + top-5 most H-like + top-5 most T-like
  # (excluding anchors from the extreme sets to avoid double-labelling).
  label_set <- function(dt, n_extreme = 5) {
    anchors <- c(lab_H, lab_T)
    non_anc <- dt[!sequence %in% anchors]
    top_H   <- non_anc[order(-dir_main)][seq_len(min(n_extreme, .N)), sequence]
    top_T   <- non_anc[order( dir_main)][seq_len(min(n_extreme, .N)), sequence]
    unique(c(anchors, top_H, top_T))
  }
  
  # ---- Per-treatment loop (all original plots) --------------------------------
  dt_list <- list()   # collect for joint panel
  
  for (tr in tr_vec) {
    
    f_csv <- file.path(path_out, paste0("ex1_1_anchors_", tr, ".csv"))
    if (!file.exists(f_csv)) {
      warning("EX1.1 diagnostics: missing input, skipping: ", f_csv)
      next
    }
    
    dt <- fread(f_csv, encoding = "UTF-8")
    dt[, sequence  := as.character(sequence)]
    dt[, dir_main  := wH_main - wT_main]
    dt[, is_anchor := sequence %in% c(lab_H, lab_T)]
    dt[, anchor_label := fcase(
      sequence == lab_H, lab_H,
      sequence == lab_T, lab_T,
      default = "Other"
    )]
    dt[, tr := tr]
    dt_list[[tr]] <- dt
    
    # ---- 0) Anchor check table -----------------------------------------------
    f_anchor <- file.path(path_out, paste0("ex1_1_anchor_check_", tr, ".csv"))
    if (!should_skip(f_anchor, cfg, "output",
                     paste0("EX1.1 anchor check (", tr, ")"))) {
      fwrite(dt[sequence %in% c(lab_H, lab_T)], f_anchor)
      msg("Saved: ", f_anchor)
    }
    
    # ---- 1) wH_main barplot --------------------------------------------------
    p1 <- ggplot(dt, aes(x = reorder(sequence, wH_main), y = wH_main)) +
      geom_col(fill = "grey55") +
      coord_flip() +
      labs(
        x     = "Sequence",
        y     = expression(paste("Hot-hand similarity weight (", w[H], ")")),
        title = paste0("EX1.1 / ", tr, " \u2014 similarity to ", lab_H)
      ) +
      theme_journal()
    save_plot(p1, paste0("ex1_1_", tr, "_wH_main.png"))
    
    # ---- 2) wT_main barplot --------------------------------------------------
    p2 <- ggplot(dt, aes(x = reorder(sequence, wT_main), y = wT_main)) +
      geom_col(fill = "grey55") +
      coord_flip() +
      labs(
        x     = "Sequence",
        y     = expression(paste("Gambler's fallacy similarity weight (", w[T], ")")),
        title = paste0("EX1.1 / ", tr, " \u2014 similarity to ", lab_T)
      ) +
      theme_journal()
    save_plot(p2, paste0("ex1_1_", tr, "_wT_main.png"))
    
    # ---- 3) Directional index barplot ----------------------------------------
    p3 <- ggplot(dt, aes(x = reorder(sequence, dir_main), y = dir_main)) +
      geom_col(fill = "grey55") +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30",
                 linewidth = 0.4) +
      coord_flip() +
      labs(
        x     = "Sequence",
        y     = expression(paste("Directional index (", w[H], " \u2212 ", w[T], ")")),
        title = paste0("EX1.1 / ", tr, " \u2014 hot-hand vs. gambler direction")
      ) +
      theme_journal()
    save_plot(p3, paste0("ex1_1_", tr, "_direction_main.png"))
    
    # ---- 4) Scatter (per-treatment, compact) ---------------------------------
    lbl_seqs <- label_set(dt)
    dt[, pt_role := fcase(
      sequence %in% c(lab_H, lab_T), "Anchor",
      sequence %in% lbl_seqs,        "Extreme",
      default = "Other"
    )]
    dt[, pt_label := fifelse(sequence %in% lbl_seqs, sequence, NA_character_)]
    
    p4 <- ggplot(dt, aes(x = wH_main, y = wT_main)) +
      geom_point(
        aes(fill = pt_role, colour = pt_role, shape = pt_role, size = pt_role),
        stroke = 0.5, alpha = 0.85
      ) +
      geom_text_repel(
        aes(label = pt_label),
        size          = 2.8,
        family        = "serif",
        colour        = "grey15",
        segment.color = "grey50",
        segment.size  = 0.3,
        segment.alpha = 0.7,
        box.padding   = 0.35,
        point.padding = 0.25,
        min.segment.length = 0.2,
        max.overlaps  = 20,
        na.rm         = TRUE
      ) +
      scale_fill_manual(values = pal_fill,
                        breaks = c("Anchor", "Extreme", "Other"),
                        labels = c("Anchor sequence", "Top-5 extreme", "Other")) +
      scale_colour_manual(values = pal_colour,
                          breaks = c("Anchor", "Extreme", "Other"),
                          labels = c("Anchor sequence", "Top-5 extreme", "Other")) +
      scale_shape_manual(values = shape_map,
                         breaks = c("Anchor", "Extreme", "Other"),
                         labels = c("Anchor sequence", "Top-5 extreme", "Other")) +
      scale_size_manual(values = size_map,
                        breaks = c("Anchor", "Extreme", "Other"),
                        labels = c("Anchor sequence", "Top-5 extreme", "Other")) +
      labs(
        x      = expression(paste("Hot-hand similarity weight (", w[H], ")")),
        y      = expression(paste("Gambler's fallacy similarity weight (", w[T], ")")),
        title  = paste0("EX1.1 / ", tr, " \u2014 anchor similarity space"),
        fill   = NULL, colour = NULL, shape = NULL, size = NULL
      ) +
      guides(fill   = guide_legend(override.aes = list(size = 3)),
             colour = guide_legend(override.aes = list(size = 3)),
             shape  = guide_legend(override.aes = list(size = 3)),
             size   = "none") +
      theme_journal()
    save_plot(p4, paste0("ex1_1_", tr, "_scatter_wH_wT.png"), w = 8, h = 7)
    
    # ---- 5) eps sensitivity --------------------------------------------------
    need <- c("wH_eps_005", "wT_eps_005", "wH_eps_003", "wT_eps_003",
              "wH_eps_008", "wT_eps_008")
    if (all(need %in% names(dt))) {
      
      ord    <- dt[order(dir_main), sequence]
      dt_long <- rbindlist(list(
        dt[, .(sequence, eps = 0.03, dir = wH_eps_003 - wT_eps_003)],
        dt[, .(sequence, eps = 0.05, dir = wH_eps_005 - wT_eps_005)],
        dt[, .(sequence, eps = 0.08, dir = wH_eps_008 - wT_eps_008)]
      ))
      dt_long[, sequence := factor(sequence, levels = ord)]
      
      p5 <- ggplot(dt_long, aes(x = eps, y = dir, group = sequence)) +
        geom_line(alpha = 0.25, colour = "grey40", linewidth = 0.4) +
        labs(
          x     = expression(paste("Tolerance (", epsilon, ")")),
          y     = expression(paste("Directional index (", w[H], " \u2212 ", w[T], ")")),
          title = paste0("EX1.1 / ", tr, " \u2014 sensitivity to \u03b5")
        ) +
        theme_journal()
      save_plot(p5, paste0("ex1_1_", tr, "_eps_sensitivity.png"), w = 7, h = 5)
      
      dir_range <- diff(range(dt$dir_main))
      if (dir_range < 0.3) {
        warning(
          "EX1.1 diagnostics: weight range for tr='", tr, "' is narrow (",
          round(dir_range, 3), "). Similarity weights lack discriminating power; ",
          "GHI results for this treatment should be interpreted with caution."
        )
      }
      
    } else {
      warning("EX1.1 diagnostics: missing eps columns for tr='", tr,
              "', skipping sensitivity plot.")
    }
  }  # end per-treatment loop
  
  # ============================================================
  # JOINT TWO-PANEL SCATTER: m25 (left) | m19 (right)
  # Journal quality — shared axes, unified legend
  # ============================================================
  
  required_tr <- c("m25", "m19")
  if (!all(required_tr %in% names(dt_list))) {
    warning("EX1.1 diagnostics: both m25 and m19 needed for joint panel; skipping.")
    return(invisible(TRUE))
  }
  
  # ---- Compute point roles and labels for each treatment --------------------
  panel_list <- lapply(required_tr, function(tr) {
    
    dt <- copy(dt_list[[tr]])
    
    lbl_seqs <- label_set(dt)
    
    # Determine top-5 H-extreme and top-5 T-extreme separately (for direction nudge)
    anchors <- c(lab_H, lab_T)
    non_anc <- dt[!sequence %in% anchors]
    top_H   <- non_anc[order(-dir_main)][seq_len(min(5L, .N)), sequence]
    top_T   <- non_anc[order( dir_main)][seq_len(min(5L, .N)), sequence]
    
    dt[, pt_role := fcase(
      sequence %in% anchors, "Anchor",
      sequence %in% c(top_H, top_T), "Extreme",
      default = "Other"
    )]
    dt[, pt_label := fifelse(sequence %in% lbl_seqs, sequence, NA_character_)]
    
    # Nudge labels: anchor + H-extreme nudge upward; T-extreme nudge downward
    dt[, nudge_x := 0]
    dt[, nudge_y := fcase(
      sequence %in% c(top_H, lab_H),  0.012,
      sequence %in% c(top_T, lab_T), -0.012,
      default = 0
    )]
    
    dt[, tr_label := tr_labels[tr]]
    dt
  })
  
  # Shared axis limits (union of both treatments)
  all_dt <- rbindlist(panel_list)
  xr <- range(all_dt$wH_main, na.rm = TRUE)
  yr <- range(all_dt$wT_main, na.rm = TRUE)
  pad <- 0.02
  xlim_shared <- c(xr[1] - pad, xr[2] + pad)
  ylim_shared <- c(yr[1] - pad, yr[2] + pad)
  
  # ---- Build one panel per treatment ----------------------------------------
  make_panel <- function(dt, is_left = TRUE) {
    
    p <- ggplot(dt, aes(x = wH_main, y = wT_main)) +
      
      # Reference lines at the anchor positions (subtle guides)
      geom_hline(yintercept = dt[sequence == lab_T, wT_main][1],
                 linetype = "dotted", colour = "grey70", linewidth = 0.35) +
      geom_vline(xintercept = dt[sequence == lab_H, wH_main][1],
                 linetype = "dotted", colour = "grey70", linewidth = 0.35) +
      
      # Points — "Other" drawn first so anchor/extreme sit on top
      geom_point(
        data = dt[pt_role == "Other"],
        aes(fill = pt_role, colour = pt_role, shape = pt_role, size = pt_role),
        stroke = 0.3, alpha = 0.60
      ) +
      geom_point(
        data = dt[pt_role == "Extreme"],
        aes(fill = pt_role, colour = pt_role, shape = pt_role, size = pt_role),
        stroke = 0.55, alpha = 0.90
      ) +
      geom_point(
        data = dt[pt_role == "Anchor"],
        aes(fill = pt_role, colour = pt_role, shape = pt_role, size = pt_role),
        stroke = 0.8, alpha = 1.00
      ) +
      
      # Labels via ggrepel
      geom_text_repel(
        aes(label = pt_label),
        size               = 2.55,
        family             = "serif",
        fontface           = "plain",
        colour             = "grey10",
        segment.color      = "grey55",
        segment.size       = 0.28,
        segment.alpha      = 0.75,
        segment.linetype   = 1,
        box.padding        = 0.40,
        point.padding      = 0.30,
        min.segment.length = 0.15,
        max.overlaps       = 30,
        seed               = 42,
        na.rm              = TRUE
      ) +
      
      # Scales
      scale_fill_manual(
        values = pal_fill,
        breaks = c("Anchor", "Extreme", "Other"),
        labels = c("Anchor (\u0048\u0048\u0048\u0048\u0048\u0048 / \u0054\u0054\u0054\u0054\u0054\u0054)",
                   "Top-5 extreme (each side)",
                   "All other sequences")
      ) +
      scale_colour_manual(
        values = pal_colour,
        breaks = c("Anchor", "Extreme", "Other"),
        labels = c("Anchor (\u0048\u0048\u0048\u0048\u0048\u0048 / \u0054\u0054\u0054\u0054\u0054\u0054)",
                   "Top-5 extreme (each side)",
                   "All other sequences")
      ) +
      scale_shape_manual(
        values = shape_map,
        breaks = c("Anchor", "Extreme", "Other"),
        labels = c("Anchor (\u0048\u0048\u0048\u0048\u0048\u0048 / \u0054\u0054\u0054\u0054\u0054\u0054)",
                   "Top-5 extreme (each side)",
                   "All other sequences")
      ) +
      scale_size_manual(
        values = size_map,
        breaks = c("Anchor", "Extreme", "Other"),
        labels = c("Anchor (\u0048\u0048\u0048\u0048\u0048\u0048 / \u0054\u0054\u0054\u0054\u0054\u0054)",
                   "Top-5 extreme (each side)",
                   "All other sequences")
      ) +
      coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
      
      labs(
        title    = unique(dt$tr_label),
        x        = expression(italic(w)[italic(H)] ~ "(hot-hand similarity)"),
        y        = expression(italic(w)[italic(T)] ~ "(gambler\u2019s fallacy similarity)"),
        fill     = NULL, colour = NULL, shape = NULL, size = NULL
      ) +
      guides(
        fill   = guide_legend(override.aes = list(size = 3.2, stroke = 0.6),
                              ncol = 1),
        colour = guide_legend(override.aes = list(size = 3.2, stroke = 0.6),
                              ncol = 1),
        shape  = guide_legend(override.aes = list(size = 3.2, stroke = 0.6),
                              ncol = 1),
        size   = "none"
      ) +
      theme_journal(base_size = 9.5)
    
    # Right panel: suppress y-axis title and labels (shared axis)
    if (!is_left) {
      p <- p + theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank()
      )
    }
    p
  }
  
  p_left  <- make_panel(panel_list[[1]], is_left = TRUE)
  p_right <- make_panel(panel_list[[2]], is_left = FALSE)
  
  # ---- Combine with patchwork -----------------------------------------------
  p_panel <- (p_left | p_right) +
    plot_layout(guides = "collect") +
    plot_annotation(
      # caption = paste0(
      #   "Notes: Each point represents one of the 64 six-toss sequences. ",
      #   "Similarity weights w\u2091 and w\u2090 measure the posterior probability ",
      #   "that a sequence\u2019s Heads-choice distribution is close (within \u03b5\u2009=\u20090.05) ",
      #   "to the pure hot-hand anchor HHHHHH and pure gambler\u2019s fallacy anchor TTTTTT, ",
      #   "respectively (normalised so w\u2091\u2009+\u2009w\u2090\u2009+\u2009w\u2080\u2009=\u20091). ",
      #   "FN: False-Negative treatment (m\u2009=\u20092.5); FP: False-Positive treatment (m\u2009=\u20091.9). ",
      #   "Labelled sequences: anchors and the five most hot-handish and five most gamblerish ",
      #   "sequences by directional index w\u2091\u2009\u2212\u2009w\u2090."
      # ),
      theme = theme(
        plot.caption     = element_text(size = 7.5, colour = "grey45",
                                        hjust = 0, lineheight = 1.35,
                                        margin = margin(t = 8)),
        plot.background  = element_rect(fill = "white", colour = NA)
      )
    ) &
    theme(legend.position = "none")
  
  # ---- Save panel (PDF for journal + PNG for drafts) -------------------------
  f_panel_pdf <- file.path(path_fig, "ex1_1_scatter_panel.pdf")
  f_panel_png <- file.path(path_fig, "ex1_1_scatter_panel.png")
  
  ggsave(f_panel_pdf, plot = p_panel, width = 14, height = 7.5,
         units = "in", device = cairo_pdf)
  msg("Saved: ", f_panel_pdf)
  
  ggsave(f_panel_png, plot = p_panel, width = 14, height = 7.5,
         units = "in", dpi = 300)
  msg("Saved: ", f_panel_png)
  
  invisible(TRUE)
}