# ============================================================
# 62_ex2-figures.R
#
# PURPOSE
#   Figure for EX2: forest plot of the outcome-specific slopes of optimism
#   and of response time on the participant-level outcomes (betting
#   probability, stake deviation, welfare loss, Heads probability), with the
#   posterior median and the 95% credible interval, for each treatment. The
#   figure is drawn from the summary tables of 61_ex2-tables.R and nothing is
#   recomputed from the fit. The pooled slopes are not drawn; they are in the
#   summary table.
#
# SAMPLE
#   full           all participants of each treatment in cfg$run$treatment
#
# INPUT (data/output/csv/)
#   ex2-<tr>-full-summary.csv
#
# OUTPUT
#   data/output/figures/ex2-slopes-forest.png
#   data/output/csv/ex2-slopes-forest-data.csv   plotted slopes
#
# NOTES
#   - Slopes are in standard deviations of the outcome per standard
#     deviation of the predictor; the dashed line marks zero
#   - Treatment labels and colours come from treatment_style()
#   - Output files are overwritten at every run
# ============================================================

ex2_figures <- function(cfg) {

  style  <- treatment_style(cfg)
  tr_vec <- style$treat

  out_lab  <- c(b = "Betting probability", a = "Stake deviation",
                c = "Welfare loss", h = "Heads probability")
  term_lab <- c(beta_opt = "Optimism", beta_rt = "Response time")

  dat <- rbindlist(lapply(tr_vec, function(tr) {
    f <- file.path(path_csv, paste0("ex2-", tr, "-full-summary.csv"))
    stopifnot(file.exists(f))
    x <- fread(f, encoding = "UTF-8")
    x <- x[block == "outcome" & term %in% names(term_lab)]
    x[, treat := tr]
    x
  }))

  dat[, `:=`(outcome_lab = factor(out_lab[as.character(outcome)], levels = rev(out_lab)),
             predictor   = factor(term_lab[term], levels = term_lab),
             treatment   = factor(style$label[match(treat, style$treat)], levels = style$label))]

  p <- ggplot(dat, aes(x = median, y = outcome_lab, colour = treatment)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(xmin = q025, xmax = q975),
                    position = position_dodge(width = 0.55), size = 0.45) +
    facet_wrap(~ predictor, nrow = 1) +
    scale_colour_manual(values = setNames(style$colour, style$label), name = NULL) +
    labs(x = "Slope (SD of the outcome per SD of the predictor)", y = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", strip.background = element_blank(),
          strip.text = element_text(face = "bold"))

  f_data <- file.path(path_csv, "ex2-slopes-forest-data.csv")
  fwrite(dat[, .(treatment = treat, predictor = as.character(predictor),
                 outcome = as.character(outcome), median, q025, q975, p_gt0)], f_data)
  msg("Saved: ", f_data)

  f_fig <- file.path(path_fig, "ex2-slopes-forest.png")
  ggsave(f_fig, p, width = 7.5, height = 3.2, dpi = 300)
  msg("Saved: ", f_fig)

  invisible(TRUE)
}
