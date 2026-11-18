# ============================================================
# 33_rq1-3-figures.R
#
# PURPOSE
#   Combined figures for RQ1 (betting), RQ2 (stake deviation) and RQ3
#   (welfare loss), drawn from the tables of 11_rq1-tables.R,
#   21_rq2-tables.R and 31_rq3-tables.R; nothing is recomputed from the fits.
#     1. Sequence forest plot with three panels (RQ3 | RQ1 | RQ2), sequences
#        ordered by the posterior mean welfare loss. A table of the five
#        sequences with the highest welfare loss reports their posterior
#        means and 95% credible intervals for the three quantities. Colours
#        give the classification of each sequence relative to the EU benchmark
#        (strong / moderate / weak / neutral departure, as in the RQ tables);
#        the vertical line and band mark the grand mean of the model.
#     2. Scatter companion: the deviations of the RQ1 and RQ2 sequence means
#        from their grand means against the rank of the sequence in welfare
#        loss, with an OLS trend fitted to the non-pure sequences. The OLS
#        slopes are exported.
#     3. Distributions of the participant-level posterior medians of the three
#        quantities (RQ1 betting probability, RQ2 stake deviation, RQ3 welfare
#        loss) as kernel densities, with the quantities in columns and three samples
#        in rows: the confirmatory sample of the FN treatment, all participants
#        of the FN treatment, and all participants of the FP treatment. A
#        summary table gives the distribution of each quantity in each sample.
#     4. Extremity of the pure sequences (one symbol repeated): their value
#        and rank among the 64 sequences in the betting rate, stake deviation
#        and welfare loss (posterior means of the models) and in the mean
#        stake of a bet and the share of all-in stakes among bets (raw data),
#        in the FN confirmatory, FN full and FP full samples. A histogram of
#        the stakes of bets compares the two pure sequences with the other
#        sequences in the same three samples.
#
# SAMPLES
#   confirmatory   treatments with cfg$design$a_flags$betting_normative = TRUE
#                  (forest plot and scatter)
#   full           all participants of a treatment; the participant
#                  distributions use the full FN and FP samples and the FN
#                  confirmatory sample
#
# INPUT (data/derived/, models)
#   master-sequences.csv                     trials (raw stakes by sequence)
#   scored-<tr>.csv (model "mpl")            Holt-Laury consistency of participants
# INPUT (data/output/csv/)
#   rq<k>-<tr>-<sample>-forest-data.csv      k = 1, 2, 3
#   rq<k>-<tr>-<sample>-participants.csv     k = 1, 2, 3 (also FP, full)
#   rq<k>-<tr>-<sample>-model-summary.csv    grand means (grand_mean_bet,
#                                            grand_mean_a, grand_mean_c)
#
# OUTPUT
#   data/output/figures/rq1-3-<tr>-confirmatory-sequences-forest.png
#   data/output/figures/rq1-3-<tr>-confirmatory-sequences-scatter.png
#   data/output/figures/rq1-3-participants-distribution.png
#   data/output/csv/rq1-3-<tr>-confirmatory-scatter-trend.csv
#   data/output/csv/rq1-3-<tr>-full-scatter-trend.csv
#                               the same trends for the full samples of each
#                               treatment, with the sequences ranked by their
#                               own welfare loss in the full sample (no figure)
#   data/output/csv/rq1-3-<tr>-confirmatory-top-sequences.csv
#                               the five sequences with the highest posterior
#                               mean welfare loss (rank, sequence, posterior
#                               mean and 95% credible interval of welfare
#                               loss, betting rate and stake deviation, and
#                               the grand-mean class of the welfare loss),
#                               followed by a row with the grand means
#   data/output/csv/rq1-3-pure-sequences.csv
#                               one row per sample, measure and pure sequence:
#                               value, rank among the 64 sequences (1 = highest
#                               value) and the mean over all 64 sequences
#   data/output/csv/rq1-3-pure-stakes-summary.csv
#                               number of bets, mean stake, and shares of stakes
#                               of 50 ECU and of the whole endowment, for the
#                               pure sequences and for the other sequences,
#                               in each of the three samples
#   data/output/figures/rq1-3-pure-stakes.png
#                               histograms of the stakes of bets, pure sequences
#                               against the other sequences, by sample
#   data/output/csv/rq1-3-participants-distribution-summary.csv
#                               one row per quantity and sample: number of
#                               participants, minimum, quartiles, maximum, mean
#                               and standard deviation of the participant
#                               posterior medians, and the grand mean of the model
#
# NOTES
#   - The sequences of the three tables must coincide; the script stops
#     otherwise
#   - Pure sequences (one symbol repeated) are labelled in the scatter plot
#     and excluded from the trend lines
#   - The participant densities are Gaussian kernel density estimates of the
#     posterior medians, reflected at the bounds of the support (0 and 1 for
#     the betting probability, -1 and 1 for the stake deviation, 0 for the
#     welfare loss). The vertical line marks the grand mean of the model
#     (posterior median)
#   - The stake deviation of the FP treatment is mechanical, since the
#     EU-optimal stake is near zero there
#   - The row labels are written vertically (FN Full, FN Confirmatory, FP
#     Full); treatment colours: _reporting-style.R
#   - Output files are overwritten at every run
# ============================================================

rq1_3_figures <- function(cfg) {

  suppressPackageStartupMessages({
    library(gridExtra); library(grid); library(scales); library(patchwork)
  })

  design <- cfg$design
  tr_vec <- unique(as.character(cfg$run$treatment))
  tr_vec <- tr_vec[vapply(tr_vec, function(tr) isTRUE(design$a_flags$betting_normative[[tr]]), logical(1))]

  read_csv <- function(rq, tr, sample, what, cls = NULL) {
    f <- file.path(path_csv, paste0(rq, "-", tr, "-", sample, "-", what, ".csv"))
    stopifnot(file.exists(f))
    fread(f, encoding = "UTF-8", colClasses = cls)
  }

  grand_levels  <- c("above", "likely_above", "neutral", "likely_below", "below")
  grand_colours <- c(above = "#C0392B", likely_above = "#E8A090", neutral = "#999999",
                     likely_below = "#90B8D8", below = "#2166AC")
  grand_names   <- c(above = "Above grand mean", likely_above = "Likely above",
                     neutral = "Equal", likely_below = "Likely below", below = "Below grand mean")

  # Classification of the sequences relative to the EU benchmark (preregistered rules of
  # the RQ tables): strong / moderate / weak departure in the direction opposite to the
  # EU prediction (underbetting, welfare loss, over-staking), or neutral. Under-staking
  # (RQ2 only) is shown in blue.
  rq_levels  <- c("strong", "moderate", "weak", "neutral",
                  "under_weak", "under_moderate", "under_strong")
  rq_colours <- c(strong = "#C0392B", moderate = "#E07B6A", weak = "#F0BFB5", neutral = "#999999",
                  under_weak = "#B8D4E8", under_moderate = "#6FA3CC", under_strong = "#2166AC")
  rq_names   <- c(strong = "Strong", moderate = "Moderate", weak = "Weak", neutral = "Neutral",
                  under_weak = "Under-staking, weak", under_moderate = "Under-staking, moderate",
                  under_strong = "Under-staking, strong")
  rq_class <- function(x) factor(sub("^over_", "", x), levels = rq_levels)

  # ============================================================
  # Scatter data and trend
  # ============================================================
  # Deviations of the RQ1 and RQ2 sequence means from their grand means, with the
  # sequences ranked by welfare loss (seq_order, highest loss first), and an OLS
  # trend of each deviation on the rank, fitted to the non-pure sequences.
  scatter_data <- function(rq1, rq2, seq_order) {
    dt <- data.table(
      sequence  = seq_order, rank_loss = seq_along(seq_order),
      bet_dev   = rq1$mu_mean[match(seq_order, as.character(rq1$sequence))] - rq1$grand_mean[1],
      stake_dev = rq2$mu_mean[match(seq_order, as.character(rq2$sequence))] - rq2$grand_mean[1]
    )
    dl <- melt(dt, id.vars = c("sequence", "rank_loss"), measure.vars = c("bet_dev", "stake_dev"),
               variable.name = "margin", value.name = "deviation")
    dl[, margin := factor(margin, levels = c("bet_dev", "stake_dev"),
                          labels = c("Betting-rate deviation", "Stake deviation"))]
    dl[, is_pure := grepl("^(.)\\1*$", sequence)]

    dtrend <- dl[is_pure == FALSE]
    trend <- rbindlist(lapply(levels(dtrend$margin), function(mg) {
      fit <- lm(deviation ~ rank_loss, data = dtrend[margin == mg])
      co  <- summary(fit)$coefficients
      data.table(margin = mg, n = nobs(fit), slope = co["rank_loss", 1], se = co["rank_loss", 2],
                 t_value = co["rank_loss", 3], p_value = co["rank_loss", 4],
                 r_squared = summary(fit)$r.squared)
    }))
    list(dl = dl, dtrend = dtrend, trend = trend)
  }

  # ============================================================
  # Forest plot and scatter companion
  # ============================================================
  forest_set <- function(tr, sample) {

    cls <- c(sequence = "character")
    rq1 <- read_csv("rq1", tr, sample, "forest-data", cls)
    rq2 <- read_csv("rq2", tr, sample, "forest-data", cls)
    rq3 <- read_csv("rq3", tr, sample, "forest-data", cls)
    stopifnot(setequal(rq1$sequence, rq3$sequence), setequal(rq2$sequence, rq3$sequence))

    flagged <- function(x) x %in% c("above", "likely_above", "below", "likely_below")
    flag_tbl <- data.table(
      sequence = rq3$sequence, mu_mean = rq3$mu_mean,
      n_flagged = flagged(rq1$grand_label[match(rq3$sequence, rq1$sequence)]) +
                  flagged(rq2$grand_label[match(rq3$sequence, rq2$sequence)])
    )
    setorder(flag_tbl, -mu_mean)
    seq_order  <- flag_tbl$sequence
    seq_levels <- rev(seq_order)                       # highest loss at the top

    prep <- function(dt, class_col) {
      dt <- copy(dt)
      dt[, sequence    := factor(sequence, levels = seq_levels)]
      dt[, grand_label := factor(grand_label, levels = grand_levels)]
      dt[, rq_label    := rq_class(get(class_col))]
      dt
    }
    rq1 <- prep(rq1, "underbet_label"); rq2 <- prep(rq2, "calib_label"); rq3 <- prep(rq3, "loss_label")
    rq_present <- rq_levels[rq_levels %in% as.character(c(rq1$rq_label, rq2$rq_label, rq3$rq_label))]

    # ---- Table: the five sequences with the highest welfare loss ----
    top_n <- 5L
    pick  <- function(dt, s_) dt[match(s_, as.character(dt$sequence))]
    top_seq <- seq_order[seq_len(top_n)]
    top <- data.table(
      rank = seq_len(top_n), sequence = top_seq,
      welfare_mean = pick(rq3, top_seq)$mu_mean,
      welfare_lo   = pick(rq3, top_seq)$mu_q025,
      welfare_hi   = pick(rq3, top_seq)$mu_q975,
      welfare_class = as.character(pick(rq3, top_seq)$grand_label),
      betting_mean = pick(rq1, top_seq)$mu_mean,
      betting_lo   = pick(rq1, top_seq)$mu_q025,
      betting_hi   = pick(rq1, top_seq)$mu_q975,
      stake_mean   = pick(rq2, top_seq)$mu_mean,
      stake_lo     = pick(rq2, top_seq)$mu_q025,
      stake_hi     = pick(rq2, top_seq)$mu_q975
    )
    gm_row <- data.table(
      rank = NA_integer_, sequence = "Grand mean",
      welfare_mean = rq3$grand_mean[1], welfare_lo = rq3$grand_lo[1], welfare_hi = rq3$grand_hi[1],
      welfare_class = NA_character_,
      betting_mean = rq1$grand_mean[1], betting_lo = rq1$grand_lo[1], betting_hi = rq1$grand_hi[1],
      stake_mean   = rq2$grand_mean[1], stake_lo   = rq2$grand_lo[1], stake_hi   = rq2$grand_hi[1]
    )
    f_top <- file.path(path_csv, paste0("rq1-3-", tr, "-", sample, "-top-sequences.csv"))
    fwrite(rbind(top, gm_row), f_top)
    msg("Saved: ", f_top)

    make_panel <- function(dt, x_label, show_y = FALSE) {
      gm <- dt$grand_mean[1]; lo <- dt$grand_lo[1]; hi <- dt$grand_hi[1]
      p <- ggplot(dt, aes(y = sequence)) +
        annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
                 fill = "steelblue", alpha = 0.08) +
        geom_vline(xintercept = gm, colour = "steelblue", linewidth = 0.55) +
        geom_segment(aes(x = mu_q025, xend = mu_q975, yend = sequence, colour = rq_label),
                     linewidth = 0.45, alpha = 0.7) +
        geom_point(aes(x = mu_mean, colour = rq_label), size = 1.8) +
        scale_colour_manual(values = rq_colours, labels = rq_names, drop = FALSE) +
        scale_x_continuous(name = x_label, labels = percent_format(accuracy = 1)) +
        annotate("text", x = gm, y = nlevels(dt$sequence) + 0.5, label = sprintf("%.1f%%", 100 * gm),
                 hjust = -0.1, size = 2.5, colour = "steelblue") +
        theme_classic(base_size = 10) +
        theme(axis.title.x = element_text(size = 8.5), axis.title.y = element_blank(),
              axis.ticks.y = element_blank(), legend.position = "none",
              panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
              plot.margin = margin(4, 8, 4, 4))
      if (show_y) {
        p + scale_y_discrete(limits = seq_levels) +
          theme(axis.text.y = element_text(family = "mono", face = "bold", size = 7, hjust = 1,
                                           colour = "black", margin = margin(r = 2)))
      } else p + theme(axis.text.y = element_blank())
    }

    p3 <- make_panel(rq3, expression(hat(mu)[s]^c ~ "(Welfare loss, share of endowment)"), TRUE)
    p1 <- make_panel(rq1, expression(hat(mu)[s]^b ~ "(Betting rate)"))
    p2 <- make_panel(rq2, expression(hat(mu)[s]^a ~ "(Stake deviation, share of endowment)"))

    legend_plot <- ggplot(data.table(x = seq_along(rq_present), y = 0,
                                     label = factor(rq_present, levels = rq_present)),
                          aes(x = x, y = y, colour = label)) +
      geom_point(size = 4) +
      scale_colour_manual(values = rq_colours[rq_present], labels = rq_names[rq_present],
                          name = "Classification relative to the EU benchmark") +
      guides(colour = guide_legend(nrow = 1, title.hjust = 0.5, override.aes = list(size = 4))) +
      coord_cartesian(xlim = c(10, 20)) + theme_void() +
      theme(legend.position = "bottom", legend.direction = "horizontal",
            legend.text = element_text(size = 10), legend.title = element_text(size = 10.5, face = "bold"),
            plot.margin = margin(0, 0, 0, 0))

    aligned <- cbind(ggplotGrob(p3), ggplotGrob(p1), ggplotGrob(p2), size = "max")   # cbind.gtable
    final <- arrangeGrob(aligned, ggplotGrob(legend_plot), ncol = 1,
                         heights = unit(c(1, 0.06), c("null", "npc")))
    stem <- paste0("rq1-3-", tr, "-", sample)
    f_forest <- file.path(path_fig, paste0(stem, "-sequences-forest.png"))
    ggsave(f_forest, plot = final, width = 17, height = 18, dpi = 300)
    msg("Saved: ", f_forest)

    # ---- Scatter companion: deviations from the grand means against loss rank ----
    sc <- scatter_data(rq1, rq2, seq_order)
    dl <- sc$dl; dtrend <- sc$dtrend; trend <- sc$trend
    margin_colours <- c("Betting-rate deviation" = "#2166AC", "Stake deviation" = "#C0392B")

    f_trend <- file.path(path_csv, paste0(stem, "-scatter-trend.csv"))
    fwrite(trend, f_trend); msg("Saved: ", f_trend)

    lab_data <- dl[is_pure == TRUE]
    lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = lab_data, aes(label = sequence), colour = "black", size = 2.8,
                               family = "mono", fontface = "bold", show.legend = FALSE,
                               min.segment.length = 0, seed = 1)
    } else geom_text(data = lab_data, aes(label = sequence), colour = "black", size = 2.8,
                     family = "mono", fontface = "bold", show.legend = FALSE, vjust = -0.8)

    ps <- ggplot(dl, aes(x = rank_loss, y = deviation, colour = margin)) +
      geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
      geom_point(aes(shape = is_pure, size = is_pure), alpha = 0.75) +
      geom_smooth(data = dtrend, aes(fill = margin), method = "lm", formula = y ~ x, se = TRUE,
                  linewidth = 0.9, alpha = 0.12) +
      lab_layer +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
      scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3.2), guide = "none") +
      scale_colour_manual(values = margin_colours, name = NULL) +
      scale_fill_manual(values = margin_colours, guide = "none") +
      scale_x_continuous(name = "Sequence rank (1 = highest welfare loss)") +
      scale_y_continuous(name = "Deviation from grand mean", labels = percent_format(accuracy = 1)) +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", axis.title = element_text(size = 9),
            plot.margin = margin(4, 8, 4, 4))
    f_sc <- file.path(path_fig, paste0(stem, "-sequences-scatter.png"))
    ggsave(f_sc, ps, width = 20, height = 10, units = "cm", dpi = 300)
    msg("Saved: ", f_sc)
  }

  # ============================================================
  # Participant distributions
  # ============================================================
  participants_distribution <- function() {

    style <- treatment_style(cfg)
    tr_fn <- style$treat[style$type == "FN"]
    tr_fp <- style$treat[style$type == "FP"]
    stopifnot(length(tr_fn) == 1L, length(tr_fp) == 1L)

    samples <- data.table(
      sample = c("FN, confirmatory", "FN, full", "FP, full"),
      tr     = c(tr_fn, tr_fn, tr_fp),
      smp    = c("confirmatory", "full", "full"),
      fill   = c(style[treat == tr_fn, colour_light], style[treat == tr_fn, colour],
                 style[treat == tr_fp, colour])
    )

    spec <- list(
      list(rq = "rq1", x = "mu_i_median", par = "grand_mean_bet", name = "betting probability",
           bw = 0.05, lim = c(0, 1), brk = seq(0, 1, by = 0.25),
           xlab = "Betting probability of a participant"),
      list(rq = "rq2", x = "mu_a_median", par = "grand_mean_a", name = "stake deviation",
           bw = 0.1, lim = c(-1, 1), brk = seq(-1, 1, by = 0.5),
           xlab = "Stake deviation of a participant (share of endowment)"),
      list(rq = "rq3", x = "mu_c_median", par = "grand_mean_c", name = "welfare loss",
           bw = 0.025, lim = c(0, 1), brk = seq(0, 1, by = 0.25),
           xlab = "Welfare loss of a participant (share of endowment)")
    )

    fills <- setNames(samples$fill, samples$sample)
    describe <- function(x) {
      data.table(n = length(x), min = min(x), q25 = unname(quantile(x, 0.25)), median = median(x),
                 q75 = unname(quantile(x, 0.75)), max = max(x), mean = mean(x),
                 sd = if (length(x) >= 2L) sd(x) else NA_real_)
    }

    # data of one quantity: participant posterior medians and model grand mean per sample
    quantity_data <- function(s) {
      rbindlist(lapply(seq_len(nrow(samples)), function(k) {
        d <- read_csv(s$rq, samples$tr[k], samples$smp[k], "participants", c(pid = "character"))
        m <- read_csv(s$rq, samples$tr[k], samples$smp[k], "model-summary")[parameter == s$par]
        stopifnot(nrow(m) == 1L)
        v <- d[[s$x]]; v <- v[is.finite(v)]
        data.table(sample = samples$sample[k], value = v, grand_mean = m$median)
      }))
    }

    # kernel densities of the participant posterior medians, reflected at the bounds of the support
    plot_col <- function(d, s, first) {
      gm <- d[, .(m = grand_mean[1], n = .N), by = sample]
      gm[, lab := sprintf("mean %.1f%%", 100 * m)]
      gm[, hj  := fifelse(m > mean(s$lim), 1.08, -0.08)]
      d[,  sample := factor(sample, levels = samples$sample)]
      gm[, sample := factor(sample, levels = samples$sample)]
      bnd <- c(s$lim[1], if (s$rq == "rq3") Inf else s$lim[2])

      ggplot(d, aes(x = value, fill = sample)) +
        geom_density(bounds = bnd, colour = "grey20", linewidth = 0.4, alpha = 0.85) +
        geom_vline(data = gm, aes(xintercept = m), colour = "grey20", linewidth = 0.7) +
        geom_text(data = gm, aes(x = m, y = Inf, label = lab, hjust = hj), inherit.aes = FALSE,
                  size = 3, colour = "grey20", vjust = 1.6) +
        geom_text(data = gm, aes(x = -Inf, y = Inf, label = paste0("n = ", n)), inherit.aes = FALSE,
                  size = 3, colour = "grey35", hjust = -0.15, vjust = 1.6) +
        scale_fill_manual(values = fills, guide = "none") +
        scale_x_continuous(name = s$xlab, breaks = s$brk,
                           labels = percent_format(accuracy = 1), expand = expansion(mult = 0.01)) +
        scale_y_continuous(name = "Density", expand = expansion(mult = c(0, 0.25))) +
        coord_cartesian(xlim = s$lim) +
        facet_grid(sample ~ ., switch = "y",
                   labeller = as_labeller(c("FN, full" = "FN Full",
                                            "FN, confirmatory" = "FN Confirmatory",
                                            "FP, full" = "FP Full"))) +
        theme_classic(base_size = 11) +
        theme(strip.placement = "outside", strip.background = element_blank(),
              strip.text.y.left = if (first) element_text(angle = 90, hjust = 0.5, size = 10)
                                  else element_blank(),
              axis.title = element_text(size = 10),
              panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
    }

    dat   <- lapply(spec, quantity_data)
    panels <- lapply(seq_along(spec), function(i) plot_col(copy(dat[[i]]), spec[[i]], first = i == 1L))

    summ <- rbindlist(lapply(seq_along(spec), function(i) {
      dat[[i]][, c(list(quantity = spec[[i]]$name), describe(value), list(grand_mean = grand_mean[1])),
               by = sample]
    }))
    summ[, sample := factor(sample, levels = samples$sample)]
    summ[, quantity := factor(quantity, levels = vapply(spec, `[[`, "", "name"))]
    setorder(summ, quantity, sample)
    f_csv <- file.path(path_csv, "rq1-3-participants-distribution-summary.csv")
    fwrite(summ, f_csv)
    msg("Saved: ", f_csv)

    fig <- panels[[1]] + panels[[2]] + panels[[3]] + plot_layout(ncol = 3)
    f <- file.path(path_fig, "rq1-3-participants-distribution.png")
    ggsave(f, fig, width = 13, height = 7.5, dpi = 300)
    msg("Saved: ", f)
  }

  # ============================================================
  # Scatter trends of the full samples (no figure)
  # ============================================================
  trend_full <- function(tr, sample) {
    cls <- c(sequence = "character")
    rq1 <- read_csv("rq1", tr, sample, "forest-data", cls)
    rq2 <- read_csv("rq2", tr, sample, "forest-data", cls)
    rq3 <- read_csv("rq3", tr, sample, "forest-data", cls)
    stopifnot(setequal(rq1$sequence, rq3$sequence), setequal(rq2$sequence, rq3$sequence))
    seq_order <- rq3[order(-mu_mean)]$sequence
    trend <- scatter_data(rq1, rq2, seq_order)$trend
    f_trend <- file.path(path_csv, paste0("rq1-3-", tr, "-", sample, "-scatter-trend.csv"))
    fwrite(trend, f_trend); msg("Saved: ", f_trend)
  }

  # ============================================================
  # Extremity of the pure sequences
  # ============================================================
  # Value and rank of the pure sequences (one symbol repeated) among all 64
  # sequences. Betting rate, stake deviation and welfare loss are the posterior
  # means of the sequence effects; the mean stake of a bet and the share of
  # all-in stakes among bets are computed from the raw trials. In the raw data,
  # the confirmatory sample consists of the Holt-Laury-consistent participants
  # of the FN treatment.
  pure_sequences <- function() {

    style <- treatment_style(cfg)
    tr_fn <- style$treat[style$type == "FN"]
    tr_fp <- style$treat[style$type == "FP"]
    stopifnot(length(tr_fn) == 1L, length(tr_fp) == 1L)
    e <- as.numeric(design$seq$endowment)

    samples <- data.table(
      sample     = c("FN, confirmatory", "FN, full", "FP, full"),
      tr         = c(tr_fn, tr_fn, tr_fp),
      smp        = c("confirmatory", "full", "full"),
      consistent = c(TRUE, FALSE, FALSE),
      fill       = c(style[treat == tr_fn, colour_light], style[treat == tr_fn, colour],
                     style[treat == tr_fp, colour])
    )

    master <- fread(file.path(path_der, "master-sequences.csv"), encoding = "UTF-8")
    master[, `:=`(pid = as.character(pid), treat = as.character(treat), stake = as.numeric(stake))]
    master[is.na(stake), stake := 0]
    cons <- rbindlist(lapply(style$treat, function(tr) {
      sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
      sc[, .(pid = as.character(pid), consistent = as.integer(inconsistent) == 0L)]
    }))
    master <- merge(master, cons, by = "pid", all.x = TRUE)

    # value and rank of the pure sequences in one measure
    pure_rows <- function(d, sample, measure, source) {
      d <- copy(d)[, .(sequence = as.character(sequence), value)]
      d[, rank := frank(-value, ties.method = "min")]
      d[, mean_all := mean(value)]
      d <- d[grepl("^(.)\\1*$", sequence)]
      d[, `:=`(sample = sample, measure = measure, source = source)]
      d
    }

    out <- rbindlist(lapply(seq_len(nrow(samples)), function(k) {
      tr <- samples$tr[k]; smp <- samples$smp[k]; nm <- samples$sample[k]
      cls <- c(sequence = "character")

      model <- list(
        list(rq = "rq1", measure = "betting rate"),
        list(rq = "rq2", measure = "stake deviation"),
        list(rq = "rq3", measure = "welfare loss")
      )
      m_rows <- rbindlist(lapply(model, function(m) {
        f <- read_csv(m$rq, tr, smp, "forest-data", cls)
        pure_rows(f[, .(sequence, value = mu_mean)], nm, m$measure, "model")
      }))

      d <- master[treat == tr]
      if (samples$consistent[k]) d <- d[consistent == TRUE]
      d <- d[stake > 0]
      raw <- d[, .(stake = mean(stake), allin = mean(stake >= e)), by = .(sequence = seq)]
      r_rows <- rbind(
        pure_rows(raw[, .(sequence, value = stake)], nm, "mean stake of a bet (ECU)", "raw"),
        pure_rows(raw[, .(sequence, value = allin)], nm, "all-in share of bets", "raw")
      )
      rbind(m_rows, r_rows)
    }))

    out[, sample := factor(sample, levels = samples$sample)]
    setorder(out, measure, sample, sequence)
    setcolorder(out, c("sample", "measure", "source", "sequence", "value", "rank", "mean_all"))
    f_csv <- file.path(path_csv, "rq1-3-pure-sequences.csv")
    fwrite(out, f_csv)
    msg("Saved: ", f_csv)

    # ---- Stakes of a bet: pure sequences against the other sequences ----
    # Histogram of the stakes of all bets in each of the three samples, in bins of
    # 5 ECU (right-closed, so that stakes of 50 and 100 fall in the bins ending there),
    # shown as the share of the bets of each group
    bw <- 5
    d <- rbindlist(lapply(seq_len(nrow(samples)), function(k) {
      x <- master[treat == samples$tr[k] & stake > 0]
      if (samples$consistent[k]) x <- x[consistent == TRUE]
      x[, sample := samples$sample[k]]
      x
    }))
    d[, group := factor(fifelse(grepl("^(.)\\1*$", seq), "Pure sequences", "Other sequences"),
                        levels = c("Pure sequences", "Other sequences"))]
    d[, sample := factor(sample, levels = samples$sample)]
    d[, bin := ceiling(stake / bw) * bw - bw / 2]

    h <- d[, .N, by = .(sample, group, bin)]
    h[, share := N / sum(N), by = .(sample, group)]

    stats <- d[, .(bets = .N, mean_stake = mean(stake),
                   share_at_50 = mean(stake == 50), share_at_100 = mean(stake >= e)),
               by = .(sample, group)]
    setorder(stats, sample, group)
    f_stats <- file.path(path_csv, "rq1-3-pure-stakes-summary.csv")
    fwrite(stats, f_stats)
    msg("Saved: ", f_stats)

    fills <- setNames(samples$fill, samples$sample)
    lab <- copy(stats)
    lab[, text := sprintf("n = %s bets, mean %.1f ECU", format(bets, big.mark = ","), mean_stake)]

    g <- ggplot(h, aes(x = bin, y = share, fill = sample)) +
      geom_col(width = bw, colour = "grey20", linewidth = 0.25, alpha = 0.9) +
      geom_vline(data = stats, aes(xintercept = mean_stake), colour = "grey20", linewidth = 0.6) +
      geom_text(data = lab, aes(x = -Inf, y = Inf, label = text), inherit.aes = FALSE,
                size = 2.8, colour = "grey25", hjust = -0.05, vjust = 1.6) +
      scale_fill_manual(values = fills, guide = "none") +
      scale_x_continuous(name = "Stake of a bet (ECU)", breaks = seq(0, e, by = e / 5),
                         expand = expansion(mult = 0.01)) +
      scale_y_continuous(name = "Share of bets", labels = percent_format(accuracy = 1),
                         expand = expansion(mult = c(0, 0.15))) +
      facet_grid(group ~ sample, switch = "y",
                 labeller = labeller(sample = c("FN, confirmatory" = "FN Confirmatory",
                                                "FN, full" = "FN Full", "FP, full" = "FP Full"))) +
      theme_classic(base_size = 11) +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 90, hjust = 0.5, size = 10),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
    f_png <- file.path(path_fig, "rq1-3-pure-stakes.png")
    ggsave(f_png, g, width = 13, height = 6.5, dpi = 300)
    msg("Saved: ", f_png)
  }

  for (tr in tr_vec) {
    forest_set(tr, "confirmatory")
  }

  for (tr in treatment_style(cfg)$treat) {
    trend_full(tr, "full")
  }

  participants_distribution()
  pure_sequences()

  invisible(NULL)
}
