# ============================================================
# 25_rq2_stake_vs_r.R
#
# PURPOSE
#   Checks whether the bimodal split in individual stake deviation
#   (mu_a_median from 23_rq2_tables.R: under_solid vs. over_solid
#   participants, "the hat") is substantially explained by
#   heterogeneity in elicited CRRA risk aversion (r_i), rather than
#   an independent behavioral split. Correlates r_i against mu_a_median
#   at the participant level (Pearson + Spearman), and reports group
#   means of r_i by calib_label.
#
#   Motivates the footnote in Sec:EconConsequences ("Individual
#   Effects") flagging this mechanical link, and the pooled r_i
#   descriptive figure (03_r_raincloud.R).
#
# INPUT
#   path_out/rq2_<tr>_<tag>_participants.csv  -- mu_a_median, calib_label
#   path_out/mpl_scored_<tr>.csv               -- r_mean, inconsistent
#
# OUTPUT
#   path_out/rq2_stakes_vs_r_<tr>_<tag>_summary.csv  -- one row: Pearson r,
#     Spearman rho (+ p-values), n
#   path_out/rq2_stakes_vs_r_<tr>_<tag>_by_label.csv -- mean r_i (+ SD, n)
#     by calib_label group
#
# NOTES
#   - Merged by pid, inner join: the participant set is whatever the
#     rq2 <tag> table already restricts to (e.g. confirmatory), so no
#     separate consistent-only file is read here.
#   - Both tags are attempted per treatment; confirmatory is skipped
#     for a treatment with no normative-better classification, same
#     gate as 23_rq2_tables.R.
# ============================================================

rq2_stake_vs_r <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  design <- cfg$design
  tags   <- c("full", "confirmatory")
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  
  for (tr in tr_vec) {
    for (tag in tags) {
      
      if (tag == "confirmatory" && !isTRUE(design$a_flags$betting_normative[[tr]])) next
      
      f_part <- file.path(path_out, paste0("rq2_", tr, "_", tag, "_participants.csv"))
      f_mpl  <- file.path(path_out, paste0("mpl_scored_", tr, ".csv"))
      
      if (!file.exists(f_part) || !file.exists(f_mpl)) {
        msg("rq2_stake_vs_r (", tr, "/", tag, "): missing input(s), skipping -- ",
            "need ", f_part, " and ", f_mpl, ".")
        next
      }
      
      part <- fread(f_part)[, .(pid = as.character(pid), mu_a_median, calib_label)]
      mpl  <- fread(f_mpl)[,  .(pid = as.character(pid), r_mean)]
      
      dt <- merge(part, mpl, by = "pid")
      if (nrow(dt) < 3L) {
        msg("rq2_stake_vs_r (", tr, "/", tag, "): fewer than 3 matched participants, skipping.")
        next
      }
      
      f_summary  <- file.path(path_out, paste0("rq2_stakes_vs_r_", tr, "_", tag, "_summary.csv"))
      f_by_label <- file.path(path_out, paste0("rq2_stakes_vs_r_", tr, "_", tag, "_by_label.csv"))
      
      if (!should_skip(f_summary, cfg, "output",
                       paste0("RQ2 stake-vs-r correlation summary (", tr, "/", tag, ")"))) {
        
        pear <- cor.test(dt$r_mean, dt$mu_a_median, method = "pearson")
        spear <- suppressWarnings(
          cor.test(dt$r_mean, dt$mu_a_median, method = "spearman", exact = FALSE)
        )
        
        summary_tbl <- data.table(
          treatment      = tr,
          tag            = tag,
          n              = nrow(dt),
          pearson_r      = unname(pear$estimate),
          pearson_p      = pear$p.value,
          spearman_rho   = unname(spear$estimate),
          spearman_p     = spear$p.value
        )
        fwrite(summary_tbl, f_summary)
        msg("Saved: ", f_summary)
      }
      
      if (!should_skip(f_by_label, cfg, "output",
                       paste0("RQ2 stake-vs-r by calib_label (", tr, "/", tag, ")"))) {
        
        by_label <- dt[, .(
          n           = .N,
          r_mean_avg  = mean(r_mean),
          r_mean_sd   = sd(r_mean),
          mu_a_avg    = mean(mu_a_median)
        ), by = calib_label]
        setorder(by_label, calib_label)
        
        fwrite(by_label, f_by_label)
        msg("Saved: ", f_by_label)
      }
    }
  }
  
  invisible(TRUE)
}