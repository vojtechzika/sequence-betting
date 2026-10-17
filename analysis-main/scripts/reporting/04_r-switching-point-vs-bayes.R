# ============================================================
# 04_r-switching-point-vs-bayes.R
#
# PURPOSE
#   Compares the classical Holt-Laury switching-point estimate of the
#   CRRA parameter r_i with the posterior mean of r_i from the Bayesian
#   MPL model. The comparison covers consistent participants, for whom
#   a switching point exists, and is reported for each treatment fit
#   and for the pooled fit, in two sets: all consistent participants
#   and interval-identified participants only (excluding those who
#   always choose the safe or the risky lottery, whose switching-point
#   estimate is a boundary value rather than an interval midpoint).
#
# INPUT
#   path_der/switching-point-r.csv
#       -- switching-point estimates (22_score-switching-point-r.R)
#   models/<selected "mpl" folder>/scored-<tag>.csv
#       -- posterior means of r_i per fit tag (01_mpl-stan-bernoulli.R);
#          the MPL folder is resolved through selected-models.csv
#
# OUTPUT
#   path_csv/r-switching-point-vs-bayes.csv
#       -- n, Pearson (with 95% CI) and Spearman correlations, mean of
#          each estimate, mean and SD of the difference (Bayes minus
#          switching point) and, for interval-identified participants, the
#          share whose posterior mean lies inside the switching-point
#          interval, by fit and participant set
#   path_fig/r-switching-point-vs-bayes.png
#       -- scatter plot of the two estimates by treatment
#
# NOTES
#   - The switching-point estimate of participants who always choose
#     the risky (safe) lottery is the boundary of the implied r
#     interval shifted by -0.5 (+0.5); these participants are marked in
#     the figure
#   - Consistency is taken from the flag in the scored file of the
#     respective fit, which must agree with the switching-point file
#   - Treatments are those in cfg$design$seq$treatments
# ============================================================

r_switching_point_vs_bayes <- function(cfg) {
  
  treatments <- names(cfg$design$seq$treatments)
  fits       <- c(treatments, "pooled")
  
  f_sw <- file.path(path_der, "switching-point-r.csv")
  stopifnot(file.exists(f_sw))
  
  sw <- fread(f_sw, encoding = "UTF-8")
  sw[, pid := as.character(pid)]
  sw <- sw[inconsistent == 0L, .(pid, treat, r_sw = r_hl, r_lo, r_hi, note)]
  
  # ---- Merge posterior means with switching-point estimates ----
  dt <- rbindlist(lapply(fits, function(tag) {
    f_sc <- model_path("mpl", paste0("scored-", gsub("_", "-", tag), ".csv"))
    sc   <- fread(f_sc, encoding = "UTF-8")
    sc[, pid := as.character(pid)]
    sc   <- sc[inconsistent == 0L, .(pid, r_bayes = r_mean)]
    
    # the fit of a treatment contains the participants of that treatment only
    sw_fit <- if (tag %in% treatments) sw[treat == tag] else sw
    d <- merge(sc, sw_fit[, !"treat"], by = "pid", all = TRUE)
    if (anyNA(d$r_bayes) || anyNA(d$r_sw)) {
      stop("Consistent participants differ between the scored file of '", tag,
           "' and switching-point-r.csv.")
    }
    d[, fit := tag]
    d
  }))
  
  # ---- Agreement statistics ----
  agreement <- function(d, fit_label, set_label) {
    ct <- cor.test(d$r_sw, d$r_bayes, method = "pearson")
    data.table(
      fit      = fit_label,
      set      = set_label,
      n        = nrow(d),
      pearson  = unname(ct$estimate),
      p_lo     = ct$conf.int[1],
      p_hi     = ct$conf.int[2],
      spearman = cor(d$r_sw, d$r_bayes, method = "spearman"),
      mean_sw    = mean(d$r_sw),
      mean_bayes = mean(d$r_bayes),
      mean_diff  = mean(d$r_bayes - d$r_sw),
      sd_diff    = sd(d$r_bayes - d$r_sw),
      # only for interval-identified participants, whose switching point
      # implies an interval [r_lo, r_hi]
      share_inside_interval = if (identical(set_label, "interval-identified"))
        mean(d$r_bayes >= d$r_lo & d$r_bayes <= d$r_hi) else NA_real_
    )
  }
  
  res <- rbindlist(lapply(fits, function(tag) {
    d <- dt[fit == tag]
    rbind(
      agreement(d,                 tag, "consistent"),
      agreement(d[note == "ok"],   tag, "interval-identified")
    )
  }))
  
  f_csv <- file.path(path_csv, "r-switching-point-vs-bayes.csv")
  fwrite(res, f_csv)
  msg("Saved: ", f_csv)
  
  # ---- Figure: both estimates by treatment ----
  pd <- dt[fit %in% treatments]
  pd[, fit  := factor(fit, levels = treatments)]
  pd[, type := fifelse(note == "ok", "Interval-identified",
                       fifelse(note == "always safe", "Always safe",
                               "Always risky"))]
  
  # Participants who always choose one lottery have no switching point; they
  # are drawn at the boundary values assigned in 22_score-switching-point-r.R,
  # which are labelled on the axis, and the number of participants is printed
  ext <- pd[type != "Interval-identified",
            .(r_sw = mean(r_sw), r_bayes = mean(r_bayes), n = .N), by = .(fit, type)]
  x_risky <- pd[type == "Always risky", unique(r_sw)]
  x_safe  <- pd[type == "Always safe",  unique(r_sw)]
  stopifnot(length(x_risky) == 1L, length(x_safe) == 1L)
  x_mid   <- pretty(range(pd[type == "Interval-identified", r_sw]), n = 4)
  x_mid   <- x_mid[x_mid > x_risky + 0.5 & x_mid < x_safe - 0.5]
  
  g <- ggplot(pd, aes(r_sw, r_bayes, colour = type, shape = type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_text(data = ext, aes(label = paste0("n = ", n)), nudge_y = 0.18, size = 2.6,
              show.legend = FALSE) +
    facet_wrap(~ fit) +
    scale_x_continuous(breaks = c(x_risky, x_mid, x_safe),
                       labels = c("always\nrisky", x_mid, "always\nsafe")) +
    scale_colour_manual(values = c("Interval-identified" = "grey20",
                                   "Always safe"         = "#1b7837",
                                   "Always risky"        = "#b2182b")) +
    labs(x = "Switching-point r", y = "Posterior mean of r (Bayesian MPL)",
         colour = NULL, shape = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  f_fig <- file.path(path_fig, "r-switching-point-vs-bayes.png")
  ggsave(f_fig, g, width = 7, height = 3.6, dpi = 300)
  msg("Saved: ", f_fig)
  
  invisible(res)
}
