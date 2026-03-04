# ============================================================
# scripts/analysis/52_ex1_1_diagnostics.R
#
# Diagnostics for EX1.I.1 anchor similarity outputs
#   - Reads:   ex1_1_anchors_<tr>.csv  (from 51_ex1_1_anchors.R)
#   - Writes:  figures into output dir
#
# Produces:
#   1) wH_main barplot (should rank HHHHHH high)
#   2) wT_main barplot (should rank OOOOOO high)
#   3) directional index (wH_main - wT_main) barplot
#   4) scatter: (wH_main, wT_main)
#   5) eps sensitivity: dir at eps=0.03/0.05/0.08 (line plot)
#   6) anchor table printout to console + CSV
#
# Notes:
# - Assumes the pipeline provides: cfg, path_out_ds(), msg()
# - Avoids heavy validation; intended as a downstream diagnostic script.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ex1_1_diagnostics <- function(cfg) {
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$run$treatment))
  
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  # helper: save plot
  save_plot <- function(p, fn, w = 10, h = 7) {
    ggsave(filename = file.path(fig_dir, fn), plot = p, width = w, height = h, units = "in")
    msg("Saved: ", file.path(fig_dir, fn))
  }
  
  for (tr in tr_vec) {
    
    f_csv <- file.path(out_dir, paste0("ex1_1_anchors_", tr, ".csv"))
    if (!file.exists(f_csv)) {
      msg("EX1.1 diagnostics: missing input, skipping: ", f_csv)
      next
    }
    
    dt <- fread(f_csv, encoding = "UTF-8")
    
    # main columns expected from 51_ex1_1_anchors.R
    dt[, sequence := as.character(sequence)]
    dt[, dir_main := wH_main - wT_main]
    
    # anchor labels (from design)
    lab_H <- as.character(cfg$design$seq$anchor_labels$pure_heads)
    lab_T <- as.character(cfg$design$seq$anchor_labels$pure_tails)
    
    # ------------------------------------------------------------
    # 0) Anchor check table (console + csv)
    # ------------------------------------------------------------
    anchor_tbl <- dt[sequence %in% c(lab_H, lab_T)]
    f_anchor <- file.path(out_dir, paste0("ex1_1_anchor_check_", tr, ".csv"))
    fwrite(anchor_tbl, f_anchor)
    msg("Saved: ", f_anchor)
    
    msg("EX1.1 anchor check (", ds, "/", tr, "):")
    print(anchor_tbl)
    
    # ------------------------------------------------------------
    # 1) wH_main barplot
    # ------------------------------------------------------------
    p1 <- ggplot(dt, aes(x = reorder(sequence, wH_main), y = wH_main)) +
      geom_col() +
      coord_flip() +
      labs(
        x = "Sequence",
        y = "Hot-hand similarity weight (wH_main)",
        title = paste0("EX1.1 / ", ds, " / ", tr, " — similarity to hot anchor (", lab_H, ")")
      )
    save_plot(p1, paste0("ex1_1_", tr, "_wH_main.png"))
    
    # ------------------------------------------------------------
    # 2) wT_main barplot
    # ------------------------------------------------------------
    p2 <- ggplot(dt, aes(x = reorder(sequence, wT_main), y = wT_main)) +
      geom_col() +
      coord_flip() +
      labs(
        x = "Sequence",
        y = "Gambler similarity weight (wT_main)",
        title = paste0("EX1.1 / ", ds, " / ", tr, " — similarity to gambler anchor (", lab_T, ")")
      )
    save_plot(p2, paste0("ex1_1_", tr, "_wT_main.png"))
    
    # ------------------------------------------------------------
    # 3) Directional index barplot
    # ------------------------------------------------------------
    p3 <- ggplot(dt, aes(x = reorder(sequence, dir_main), y = dir_main)) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = "Sequence",
        y = "Directional index (wH_main − wT_main)",
        title = paste0("EX1.1 / ", ds, " / ", tr, " — hot vs gambler direction")
      )
    save_plot(p3, paste0("ex1_1_", tr, "_direction_main.png"))
    
    # ------------------------------------------------------------
    # 4) Scatter: (wH_main, wT_main)
    # ------------------------------------------------------------
    dt[, is_anchor := sequence %in% c(lab_H, lab_T)]
    p4 <- ggplot(dt, aes(x = wH_main, y = wT_main)) +
      geom_point(aes(shape = is_anchor), size = 2) +
      labs(
        x = "Hot-hand similarity (wH_main)",
        y = "Gambler similarity (wT_main)",
        title = paste0("EX1.1 / ", ds, " / ", tr, " — anchor similarity space"),
        shape = "Anchor"
      )
    save_plot(p4, paste0("ex1_1_", tr, "_scatter_wH_wT.png"))
    
    # ------------------------------------------------------------
    # 5) eps sensitivity: dir at eps=0.03/0.05/0.08
    # ------------------------------------------------------------
    # expected columns: wH_eps_003 / 005 / 008 etc (naming is gsub(".", "", sprintf("%.2f", eps)))
    need <- c("wH_eps_005", "wT_eps_005", "wH_eps_003", "wT_eps_003", "wH_eps_008", "wT_eps_008")
    if (all(need %in% names(dt))) {
      
      dt_long <- rbindlist(list(
        dt[, .(sequence, eps = 0.03, dir = wH_eps_003 - wT_eps_003)],
        dt[, .(sequence, eps = 0.05, dir = wH_eps_005 - wT_eps_005)],
        dt[, .(sequence, eps = 0.08, dir = wH_eps_008 - wT_eps_008)]
      ))
      
      # order sequences by main direction for readability
      ord <- dt[order(dir_main), sequence]
      dt_long[, sequence := factor(sequence, levels = ord)]
      
      p5 <- ggplot(dt_long, aes(x = eps, y = dir, group = sequence)) +
        geom_line(alpha = 0.25) +
        labs(
          x = "epsilon (tolerance)",
          y = "Directional index (wH − wT)",
          title = paste0("EX1.1 / ", ds, " / ", tr, " — sensitivity to epsilon")
        )
      save_plot(p5, paste0("ex1_1_", tr, "_eps_sensitivity.png"))
      
    } else {
      msg("EX1.1 diagnostics: skipping eps sensitivity (missing eps columns) for tr=", tr)
    }
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Inspect anchor probabilities from RQ4 output
# ------------------------------------------------------------



# Example:
# ex1_1_diagnostics(cfg)