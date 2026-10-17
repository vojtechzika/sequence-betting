# ============================================================
# 72_ex3-figures.R
#
# PURPOSE
#   Figure for EX3: the participant-level posterior medians of betting
#   probability (RQ1), proportional stake deviation (RQ2) and welfare loss
#   (RQ3) for the full samples of both treatments, one panel per quantity,
#   treatments overlaid as raincloud distributions with their grand medians.
#   Drawn from the participant tables of the RQ scripts; nothing is
#   recomputed from the fits.
#
# SAMPLE
#   The full sample of each treatment in cfg$run$treatment.
#
# INPUT (data/output/csv/)
#   rq<k>-<tr>-full-participants.csv     k = 1, 2, 3 (mu_i_median, mu_a_median,
#                                        mu_c_median)
#
# OUTPUT
#   data/output/figures/ex3-participants-raincloud.png
#
# NOTES
#   - Treatment labels and colours come from _reporting-style.R
#   - The density of the welfare loss is bounded at 0 (ggdist >= 3.3.0)
#   - The RQ3 participant table is that of the model selected for "rq3";
#     the contrast in 71_ex3-contrasts.R uses the hurdle-Gamma model
#   - Output files are overwritten at every run
# ============================================================

ex3_figures <- function(cfg) {

  if (utils::packageVersion("ggdist") < "3.3.0")
    stop("72_ex3-figures.R requires ggdist >= 3.3.0 (bounded densities).", call. = FALSE)

  tr_vec <- unique(as.character(cfg$run$treatment))
  style  <- treatment_style(cfg)
  cols   <- setNames(style$colour, style$label)

  spec <- list(
    list(rq = "rq1", x = "mu_i_median", bounded = FALSE,
         lab = expression(hat(mu)[i]^b ~ "(betting probability)")),
    list(rq = "rq2", x = "mu_a_median", bounded = FALSE,
         lab = expression(hat(mu)[i]^a ~ "(proportional stake deviation)")),
    list(rq = "rq3", x = "mu_c_median", bounded = TRUE,
         lab = expression(hat(mu)[i]^c ~ "(welfare loss, share of endowment)"))
  )

  panels <- lapply(spec, function(s) {
    d <- rbindlist(lapply(tr_vec, function(tr) {
      f <- file.path(path_csv, paste0(s$rq, "-", tr, "-full-participants.csv"))
      stopifnot(file.exists(f))
      t <- fread(f, encoding = "UTF-8")[, .(value = get(s$x))]
      t[, treat_label := style$label[style$treat == tr]]
      t
    }))
    d[, treat_label := factor(treat_label, levels = style$label)]
    gm <- d[, .(m = median(value)), by = treat_label]
    dens <- if (s$bounded) ggdist::density_bounded(bounds = c(0, Inf)) else ggdist::density_unbounded()

    ggplot(d, aes(x = value, y = treat_label, fill = treat_label, colour = treat_label)) +
      stat_halfeye(density = dens, adjust = 0.8, width = 0.6, .width = 0, point_colour = NA,
                   alpha = 0.55, colour = NA, position = position_nudge(y = 0.12)) +
      geom_jitter(shape = 21, size = 1.5, alpha = 0.7, stroke = 0.25, colour = "white",
                  height = 0.07, width = 0) +
      geom_point(data = gm, aes(x = m, y = treat_label), inherit.aes = FALSE,
                 shape = 124, size = 6, colour = "grey15") +
      scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
      scale_x_continuous(name = s$lab, labels = percent_format(accuracy = 1)) +
      scale_y_discrete(name = NULL) +
      theme_classic(base_size = 10) +
      theme(legend.position = "none",
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
  })

  p <- (panels[[1]] + panels[[2]]) / panels[[3]]
  f <- file.path(path_fig, "ex3-participants-raincloud.png")
  ggsave(f, p, width = 8, height = 6, dpi = 600);  msg("Saved: ", f)

  invisible(TRUE)
}
