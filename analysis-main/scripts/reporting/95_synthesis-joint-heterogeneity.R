# ============================================================
# 95_synthesis-joint-heterogeneity.R
#
# PURPOSE
#   Participant-level heterogeneity in the joint model of betting and side
#   choice (91_synthesis-stan-joint.R), per treatment.
#   (1) Personal slopes. For each participant, the personal slope is the
#       population fixed effect plus the participant's random effect, for
#       the four cue slopes (betting: run, imbalance; side: run,
#       imbalance). Each participant is classified by the 95% credible
#       interval of the personal slope:
#         side run:        against_run / with_run / uncertain
#         side imbalance:  against_imb (bets the minority symbol) /
#                          with_imb (bets the majority symbol) / uncertain
#         betting slopes:  engaged (credibly positive) / disengaged
#                          (credibly negative) / uncertain
#       "Uncertain" means that the data do not pin down the direction, not
#       that the point estimate is zero. Counts and shares are reported.
#   (2) Run length by imbalance class. Share of bets against the terminal
#       run, by side-imbalance class and run length (2 to 6 repeats), on bets
#       made on sequences whose terminal run consists of the majority symbol.
#       Beyond three repeats the run symbol is always the majority, so run
#       length is confounded with imbalance from a run of four onward and
#       only the contrast between two and three repeats is clean. The
#       second table gives paired within-participant steps between adjacent
#       run lengths, pooled over imbalance levels and within the 4-2
#       sequences only (imbalance held fixed). The classes were estimated on
#       all trials, including those used here, so the shares describe how
#       each class behaves and are not a test of the classification.
#   (3) A figure of the distributions of the personal run and imbalance
#       slopes (posterior medians), with the population fixed effects, by
#       treatment.
#
# SAMPLE
#   synthesis-full of each treatment in cfg$run$treatment.
#
# INPUT
#   model "synthesis": fit-, pid-levels-<tr>-full.rds
#   model "mpl": sample/analysis-trials.csv (through analysis_trials())
#   data/output/csv/synthesis-<tr>-joint-fixed-effects.csv (94_synthesis-joint-effects.R)
#
# OUTPUT (data/output/csv/)
#   synthesis-<tr>-joint-run-slope-participants.csv       personal side run slope and class
#   synthesis-<tr>-joint-imb-slope-participants.csv       personal side imbalance slope and class
#   synthesis-<tr>-joint-bet-run-slope-participants.csv   personal betting run slope and class
#   synthesis-<tr>-joint-bet-imb-slope-participants.csv   personal betting imbalance slope and class
#   synthesis-<tr>-joint-class-counts.csv                 counts and shares by class
#   synthesis-<tr>-joint-run-length-by-class.csv          share against the run, t-based 95% CI
#   synthesis-<tr>-joint-run-length-steps.csv             paired steps between run lengths
# OUTPUT (data/output/figures/)
#   synthesis-slope-distributions.png
#
# NOTES
#   - Confidence intervals in the run-length tables are t-based over
#     participants (NA where participants have no variance, for example at a
#     ceiling of 100% against the run); no seed is needed
#   - Steps based on fewer than three participants are omitted
#   - Output files are overwritten at every run
# ============================================================

synthesis_joint_heterogeneity <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  smp    <- "full"
  heads  <- as.character(cfg$design$seq$side_labels$heads)
  tails  <- as.character(cfg$design$seq$side_labels$tails)

  # mean and t-based 95% CI over participants (v = one value per
  # participant). Bounds are NA when the participant values have no variance
  # (e.g. everyone at 100% against the run, a ceiling) or there is only one
  # value, where a t-interval is undefined; report such cells without an
  # interval.
  t_ci <- function(v) {
    n <- length(v)
    m <- mean(v)
    s <- if (n > 1L) sd(v) else NA_real_
    h <- if (is.finite(s) && s > 0) qt(.975, n - 1L) * s / sqrt(n) else NA_real_
    c(mean = m, lo = m - h, hi = m + h)
  }
  
  # paired steps between adjacent run lengths (2->3, 3->4, 4->5, 5->6), per
  # class, on participant-level mean shares of bets against the run; steps
  # with fewer than 3 participants are omitted
  run_steps <- function(dd, scope) {
    w <- dcast(dd[, .(against = mean(against)), by = .(imb_class, run_len, pid)],
               imb_class + pid ~ run_len, value.var = "against")
    rows <- list()
    for (cl in unique(w$imb_class)) {
      wc <- w[imb_class == cl]
      for (p in list(c(2L, 3L), c(3L, 4L), c(4L, 5L), c(5L, 6L))) {
        a <- as.character(p[1]); b <- as.character(p[2])
        if (!all(c(a, b) %in% names(wc))) next
        df <- wc[[b]] - wc[[a]]
        df <- df[is.finite(df)]
        if (length(df) < 3L) next
        ci <- t_ci(df)
        rows[[length(rows) + 1L]] <- data.table(
          scope = scope, imb_class = cl, from = p[1], to = p[2],
          n_participants = length(df), step_pp = 100 * unname(ci["mean"]),
          ci_low = 100 * unname(ci["lo"]), ci_high = 100 * unname(ci["hi"]))
      }
    }
    rbindlist(rows)
  }

  for (tr in tr_vec) {

    save_csv <- function(x, name) {
      f <- file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, ".csv"))
      fwrite(x, f);  msg("Saved: ", f)
    }

    fit <- readRDS(model_path("synthesis", paste0("fit-", tr, "-", smp, ".rds")))
    pid_levels <- readRDS(model_path("synthesis", paste0("pid-levels-", tr, "-", smp, ".rds")))
    draws <- rstan::extract(fit, pars = c("beta_bet", "beta_side", "u"))

    run_slope <- draws$beta_side[, 2] + draws$u[, 5, ]   # [iterations, P]
    imb_slope <- draws$beta_side[, 3] + draws$u[, 6, ]
    
    summarize_participant <- function(m) {
      data.table(pid    = pid_levels,
                 median = apply(m, 2, median),
                 q025   = apply(m, 2, quantile, .025),
                 q975   = apply(m, 2, quantile, .975))
    }
    
    run_tbl <- summarize_participant(run_slope)
    imb_tbl <- summarize_participant(imb_slope)
    
    # classification is w.r.t. each participant's own 95% interval on
    # their personal slope (beta_side + u_i), not a point estimate --
    # "uncertain" means the data don't credibly pin down that person's
    # direction, not that their point estimate is exactly zero.
    run_tbl[, class := fifelse(q025 > 0, "with_run",
                               fifelse(q975 < 0, "against_run", "uncertain"))]
    imb_tbl[, class := fifelse(q025 > 0, "with_imb",
                               fifelse(q975 < 0, "against_imb", "uncertain"))]
    
    
    save_csv(run_tbl, "run-slope-participants")
    save_csv(imb_tbl, "imb-slope-participants")
    
    # ---- same thing for the BETTING equation's two random slopes ----
    # draws$u row 2 = bet run-slope RE, row 3 = bet imbalance-slope RE (see
    # synthesis.stan header for the full u[] index map). Both population
    # fixed effects (beta_bet[,2], beta_bet[,3]) are positive, so a
    # participant's personal slope being classified "engaged" means that
    # cue credibly INCREASES their probability of betting; "disengaged"
    # means it credibly DECREASES it (a real, if perhaps counterintuitive,
    # possibility -- e.g. someone who avoids betting on trials they read as
    # more predictable). Same 95%-own-interval logic as the side equation.
    bet_run_slope <- draws$beta_bet[, 2] + draws$u[, 2, ]   # [iterations, P]
    bet_imb_slope <- draws$beta_bet[, 3] + draws$u[, 3, ]
    
    bet_run_tbl <- summarize_participant(bet_run_slope)
    bet_imb_tbl <- summarize_participant(bet_imb_slope)
    
    bet_run_tbl[, class := fifelse(q025 > 0, "engaged", fifelse(q975 < 0, "disengaged", "uncertain"))]
    bet_imb_tbl[, class := fifelse(q025 > 0, "engaged", fifelse(q975 < 0, "disengaged", "uncertain"))]
    
    
    save_csv(bet_run_tbl, "bet-run-slope-participants")
    save_csv(bet_imb_tbl, "bet-imb-slope-participants")
    
    # class_counts: the same table() summaries,
    # above, saved to disk -- the most direct "how many participants
    # reliably go which way" numbers for a heterogeneity write-up (counts
    # + share of P, per cue, per equation-role). Now covers all four
    # random slopes (bet-run, bet-imb, side-run, side-imb), not just the
    # side equation's two.
    class_counts <- rbind(
      data.table(regressor = "bet_run_slope", class = names(table(bet_run_tbl$class)),
                 n = as.integer(table(bet_run_tbl$class))),
      data.table(regressor = "bet_imb_slope", class = names(table(bet_imb_tbl$class)),
                 n = as.integer(table(bet_imb_tbl$class))),
      data.table(regressor = "side_run_slope", class = names(table(run_tbl$class)),
                 n = as.integer(table(run_tbl$class))),
      data.table(regressor = "side_imb_slope", class = names(table(imb_tbl$class)),
                 n = as.integer(table(imb_tbl$class)))
    )
    class_counts[, share := n / length(pid_levels)]
    save_csv(class_counts, "class-counts")

    # ---- Trial data for the run-length tables ----
    mst <- analysis_trials(tr, paste0("synthesis-", smp))
    mst <- mst[is.finite(stake) & stake > 0 & side %in% c(heads, tails)]
    mst[, `:=`(pid = as.character(pid), seq = as.character(seq))]
    seq_lk <- data.table(seq = unique(mst$seq))
    seq_lk[, last_sym := substr(seq, nchar(seq), nchar(seq))]
    seq_lk[, n_h := nchar(gsub(tails, "", seq, fixed = TRUE))]
    seq_lk[, n_t := nchar(seq) - n_h]
    seq_lk[, maj := fifelse(n_h > n_t, heads, fifelse(n_t > n_h, tails, "="))]
    seq_lk[, run_len := vapply(strsplit(seq, ""), function(x) {
      n <- length(x); k <- 1L
      while (k < n && x[n - k] == x[n]) k <- k + 1L
      k
    }, integer(1))]
    mst <- merge(mst, seq_lk, by = "seq", all.x = TRUE)
    mst[, against := as.integer(side != last_sym)]
    run_trials <- mst[run_len >= 2L & last_sym == maj]

    rl_d <- merge(run_trials,
                  imb_tbl[, .(pid = as.character(pid), imb_class = class)],
                  by = "pid")
    if (nrow(rl_d) > 0L) {
      
      rl_pm <- rl_d[, .(against = mean(against), n_bets = .N),
                    by = .(imb_class, run_len, pid)]
      rl_lev <- rl_pm[, {
        ci <- t_ci(against)
        list(n_participants = .N, n_bets = sum(n_bets),
             against_pct = 100 * unname(ci["mean"]),
             ci_low = 100 * unname(ci["lo"]), ci_high = 100 * unname(ci["hi"]))
      }, by = .(imb_class, run_len)]
      setorder(rl_lev, imb_class, run_len)
      
      rl_stp <- rbind(run_steps(rl_d, "all imbalance levels"),
                      run_steps(rl_d[n_h %in% c(2L, 4L)], "4-2 sequences only"))
      
      save_csv(rl_lev, "run-length-by-class")
      save_csv(rl_stp, "run-length-steps")
    }
  }

  # ---- Figure: personal run and imbalance slopes by treatment ----
  style <- treatment_style(cfg)
  cols  <- setNames(style$colour, style$label)
  mk_long <- function(name) rbindlist(lapply(tr_vec, function(tr) {
    d <- fread(file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, ".csv")),
               encoding = "UTF-8", colClasses = c(pid = "character"))
    data.table(label = style$label[style$treat == tr], slope = d$median)
  }))
  mk_pop <- function(pname) rbindlist(lapply(tr_vec, function(tr) {
    d <- fread(file.path(path_csv, paste0("synthesis-", tr, "-joint-fixed-effects.csv")))
    data.table(label = style$label[style$treat == tr], pop = d$median[d$param == pname])
  }))
  mk_panel <- function(long, pop, xlab) {
    long[, label := factor(label, levels = style$label)]
    pop[,  label := factor(label, levels = style$label)]
    ggplot() +
      geom_vline(xintercept = 0, colour = "grey40", linetype = "dotted") +
      geom_density(data = long, aes(x = slope, colour = label, fill = label), alpha = 0.2, linewidth = 0) +
      geom_vline(data = pop, aes(xintercept = pop, colour = label), linetype = "dashed", linewidth = 0.6) +
      scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
      labs(x = xlab, y = "Density", colour = NULL, fill = NULL) +
      theme_classic(base_size = 10) +
      theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  }
  p_run <- mk_panel(mk_long("run-slope-participants"), mk_pop("side_run_slope"),
                    "Participant run slope (posterior median)")
  p_imb <- mk_panel(mk_long("imb-slope-participants"), mk_pop("side_imb_slope"),
                    "Participant imbalance slope (posterior median)")
  fig <- p_run + p_imb + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
  f_fig <- file.path(path_fig, "synthesis-slope-distributions.png")
  ggsave(f_fig, fig, width = 9, height = 3, dpi = 300);  msg("Saved: ", f_fig)

  invisible(TRUE)
}
