# ============================================================
# 51_ex1-weights-scatter.R
#
# PURPOSE
#   Documents the anchor-based similarity weights of EX1 and shows why the
#   framework does not discriminate between sequences. For each sequence s
#   the posterior similarity to the anchors H = HHHHHH, T = TTTTTT and to the
#   neutral baseline 0 (the population baseline hbar) is
#       d_a(s) = P( |theta_s - theta_a| < eps | data ),   a in {H, T, 0},
#   with theta_s = mu_h[s] of the RQ4 model, and the normalised weights are
#       w_a(s) = (d_a(s) + eta) / (d_H(s) + d_T(s) + d_0(s) + 3 eta).
#   These are the collapsed (across draws) quantities of the supplement; the
#   draw-level version used for the indices is in 51_ex1-chi-draws.R.
#   The figure plots w_H against w_T for every sequence, one panel per
#   treatment, and marks the anchors and the most extreme sequences.
#
# SAMPLE
#   rq4-full of each treatment in cfg$run$treatment (as the RQ4 model).
#
# INPUT
#   model "rq4": fit-<tr>-full.rds, seq-levels-<tr>-full.rds   (mu_h, hbar)
#   cfg$design$ex1 (eps, eta), cfg$design$seq$anchor_labels
#
# OUTPUT
#   data/output/csv/ex1-<tr>-weights.csv        d, w and w_H - w_T per sequence
#                                               and tolerance (long format)
#   data/output/csv/ex1-<tr>-anchor-check.csv   the two anchors at the main
#                                               tolerance
#   data/output/figures/ex1-weights-scatter.png
#
# NOTES
#   - The first value of cfg$design$ex1$eps is the main tolerance and is the
#     one drawn; the others are reported in the csv as sensitivity
#   - Labelled sequences: the two anchors and the five sequences with the
#     largest and the five with the smallest w_H - w_T
#   - Output files are overwritten at every run
# ============================================================

ex1_weights_scatter <- function(cfg) {

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))
  style  <- treatment_style(cfg)

  eps_vec  <- as.numeric(design$ex1$eps)
  eps_main <- eps_vec[1]
  eta      <- as.numeric(design$ex1$eta)
  stopifnot(length(eps_vec) >= 1L, all(eps_vec > 0), length(eta) == 1L, eta > 0)

  lab_H <- as.character(design$seq$anchor_labels$pure_heads)
  lab_T <- as.character(design$seq$anchor_labels$pure_tails)

  panel_data <- list()

  for (tr in tr_vec) {

    fit        <- readRDS(model_path("rq4", paste0("fit-", tr, "-full.rds")))
    seq_levels <- as.character(readRDS(model_path("rq4", paste0("seq-levels-", tr, "-full.rds"))))
    post       <- rstan::extract(fit, pars = c("mu_h", "hbar"))
    theta_s    <- post$mu_h
    theta_0    <- as.numeric(post$hbar)

    iH <- match(lab_H, seq_levels); iT <- match(lab_T, seq_levels)
    if (is.na(iH) || is.na(iT)) stop("EX1: anchor missing in the RQ4 fit (", tr, ").", call. = FALSE)
    theta_H <- theta_s[, iH]
    theta_T <- theta_s[, iT]

    out <- rbindlist(lapply(eps_vec, function(eps) {
      d_H <- colMeans(abs(theta_s - theta_H) < eps)
      d_T <- colMeans(abs(theta_s - theta_T) < eps)
      d_0 <- colMeans(abs(theta_s - theta_0) < eps)
      den <- d_H + d_T + d_0 + 3 * eta
      data.table(treatment = tr, sample = "full", eps = eps, sequence = seq_levels,
                 d_H = d_H, d_T = d_T, d_0 = d_0,
                 w_H = (d_H + eta) / den, w_T = (d_T + eta) / den, w_0 = (d_0 + eta) / den,
                 direction = (d_H - d_T) / den)
    }))
    f_w <- file.path(path_csv, paste0("ex1-", tr, "-weights.csv"))
    fwrite(out, f_w);  msg("Saved: ", f_w)

    main <- out[eps == eps_main]
    f_a <- file.path(path_csv, paste0("ex1-", tr, "-anchor-check.csv"))
    fwrite(main[sequence %in% c(lab_H, lab_T)], f_a);  msg("Saved: ", f_a)

    # point roles: anchors, five most extreme on each side, others
    non_anc <- main[!sequence %in% c(lab_H, lab_T)]
    top_H <- non_anc[order(-direction)][seq_len(5L), sequence]
    top_T <- non_anc[order(direction)][seq_len(5L), sequence]
    main[, role := fcase(sequence %in% c(lab_H, lab_T), "Anchor",
                         sequence %in% c(top_H, top_T), "Extreme",
                         default = "Other")]
    main[, label := fifelse(role != "Other", sequence, NA_character_)]
    panel_label <- style$label[style$treat == tr]
    main[, panel := panel_label]
    panel_data[[tr]] <- main
  }

  all_dt <- rbindlist(panel_data)
  all_dt[, panel := factor(panel, levels = style$label)]
  all_dt[, role := factor(role, levels = c("Anchor", "Extreme", "Other"))]
  pad <- 0.02
  xl <- range(all_dt$w_H) + c(-pad, pad)
  yl <- range(all_dt$w_T) + c(-pad, pad)

  p <- ggplot(all_dt, aes(w_H, w_T)) +
    geom_point(data = all_dt[role == "Other"], colour = "#8494A8", fill = "#AAB4C4",
               shape = 21, size = 1.6, alpha = 0.6) +
    geom_point(data = all_dt[role == "Extreme"], colour = "#9B2819", fill = "#C0392B",
               shape = 21, size = 2.4, alpha = 0.9) +
    geom_point(data = all_dt[role == "Anchor"], colour = "#1A1A2E", fill = "#1A1A2E",
               shape = 23, size = 3.2) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.6, colour = "grey10",
                             segment.colour = "grey55", segment.size = 0.3,
                             box.padding = 0.4, min.segment.length = 0.15,
                             max.overlaps = 30, seed = 42, na.rm = TRUE) +
    coord_cartesian(xlim = xl, ylim = yl) +
    facet_wrap(~ panel) +
    labs(x = expression(italic(w)[H] ~ "(similarity to HHHHHH)"),
         y = expression(italic(w)[T] ~ "(similarity to TTTTTT)")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey96", colour = "grey70"),
          strip.text = element_text(face = "bold"))

  f_png <- file.path(path_fig, "ex1-weights-scatter.png")
  ggsave(f_png, p, width = 11, height = 5.5, dpi = 300);  msg("Saved: ", f_png)

  invisible(TRUE)
}
