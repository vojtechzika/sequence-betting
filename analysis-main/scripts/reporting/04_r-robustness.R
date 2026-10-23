# ============================================================
# 04_r-robustness.R
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
#   A second table describes the choice patterns of all participants by
#   Holt-Laury consistency: the number of safe choices and of switches
#   between the lotteries, the share who choose the risky lottery in the
#   first row and the safe lottery in the last row, and the mean posterior
#   r. In the first row the safe lottery is better for every r above a
#   threshold, and in the last row the risky lottery is better for every r
#   below one; both thresholds are computed from the design of the task and
#   reported, so that the shares can be read against them.
#
#   A third table reports the choice precision lambda_i of the Bayesian MPL
#   model by treatment and Holt-Laury consistency: the average posterior
#   mean of log(lambda_i), next to the prior value log(30). A low precision
#   means that the model attributes the choices of a participant to noise.
#
#   A fourth table compares the posterior mean of r_i of the primary
#   (individual-level) model, by treatment and Holt-Laury consistency, with
#   three alternative estimates: the switching-point estimate (consistent
#   participants only), the exact grid posterior mean of the same model,
#   computed here for every participant, and the posterior mean of the
#   hierarchical model (02_mpl-stan-hierarchical.R). The table also reports
#   the mean posterior SD of r_i of the primary and the hierarchical model.
#
# INPUT
#   path_der/switching-point-r.csv
#       -- switching-point estimates (22_score-switching-point-r.R)
#   models/<selected "mpl" folder>/scored-<tag>.csv
#       -- posterior means of r_i per fit tag (01_mpl-stan-bernoulli.R);
#          the MPL folder is resolved through selected-models.csv
#   path_src/holt-laury.csv
#       -- Holt-Laury choices c1..cK of every participant
#   cfg$design$mpl
#       -- labels of the safe and risky lottery, payoffs and probabilities
#   models/<selected "mpl" folder>/fit-<tr>.rds, r-draws-<tr>.rds
#       -- Stan fit of each treatment (draws of lambda_i) and the participants
#          in the order of its columns (01_mpl-stan-bernoulli.R)
#   models/<selected "mpl" folder>/exact-draw-participants.csv
#       -- participants whose draws of r_i were replaced (01_mpl-stan-bernoulli.R)
#   models/mpl-hierarchical/scored-<tr>.csv, scored-<tr>-consistent.csv,
#   r-draws-<tr>.rds, r-draws-<tr>-consistent.rds
#       -- posterior means and draws of r_i of the hierarchical model, fitted
#          to the full treatment sample and to its consistent participants
#          (02_mpl-stan-hierarchical.R)
#   models/<selected "mpl" folder>/r-draws-<tr>.rds
#       -- draws of r_i of the primary model (01_mpl-stan-bernoulli.R)
#   cfg$model$stan$mpl_bernoulli$prior
#       -- priors of the primary model, used for the grid posterior
#
# OUTPUT
#   path_csv/r-switching-point-vs-bayes.csv
#       -- n, Pearson (with 95% CI) and Spearman correlations, mean of
#          each estimate, mean and SD of the difference (Bayes minus
#          switching point) and, for interval-identified participants, the
#          share whose posterior mean lies inside the switching-point
#          interval, by fit and participant set
#   path_csv/r-choice-patterns-by-consistency.csv
#       -- by treatment and group (all, consistent, inconsistent): n, mean
#          posterior r, mean number of safe choices, mean number of
#          switches between the lotteries, share choosing the risky
#          lottery in the first row, share choosing the safe lottery in the
#          last row, and the values of r at which the participant is
#          indifferent in the first and in the last row
#   path_csv/r-choice-precision-by-consistency.csv
#       -- by treatment and group (all, consistent, inconsistent): n, number
#          of participants left out, average posterior mean of log(lambda_i),
#          its SD across participants, its exponential, and the prior value
#   path_csv/r-model-comparison.csv
#       -- by treatment and group (all, consistent, inconsistent): n, mean r_i
#          of the primary model, the grid, the hierarchical model and (for
#          consistent participants) the switching point; correlations of each
#          alternative with the primary model; largest absolute difference
#          to the grid, mean absolute difference to the hierarchical model, and
#          the mean posterior SD of r_i of the primary and the hierarchical model
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
#   - The choice-pattern table uses all participants of each treatment, not
#     only the consistent ones; the group is the flag in the scored file
#   - The precision table leaves out the participants listed in
#     exact-draw-participants.csv: for them the sampler's lambda_i lies in a
#     local mode that does not belong to the posterior (the same failure that
#     requires the exact draws of r_i); their number is in n_left_out
#   - In the model comparison, the grid posterior mean is computed on a grid
#     of r with step 0.005 and of log(lambda) with step 0.1; the hierarchical
#     estimate of the consistent group comes from the fit on the consistent
#     participants and that of the other groups from the full-sample fit
# ============================================================

r_robustness <- function(cfg) {
  
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
  
  # ---- Choice patterns by consistency ----
  # How the choices respond to the design of the task, for all participants
  # and by Holt-Laury consistency. A participant who responds to the payoffs
  # chooses the safe lottery in the first row (low probability of the high
  # payoff) and the risky lottery in the last row; r_first and r_last are the
  # values of r at which the participant is indifferent in those rows.
  mpl_design  <- cfg$design$mpl
  label_safe  <- toupper(trimws(as.character(mpl_design$label_safe)))
  label_risky <- toupper(trimws(as.character(mpl_design$label_risky)))
  p_sched     <- as.numeric(mpl_design$p)
  K_rows      <- length(p_sched)
  A_high <- as.numeric(mpl_design$A_high); A_low <- as.numeric(mpl_design$A_low)
  B_high <- as.numeric(mpl_design$B_high); B_low <- as.numeric(mpl_design$B_low)
  
  crra_u  <- function(x, r) ifelse(abs(r - 1) < 1e-8, log(x), x^(1 - r) / (1 - r))
  eu_diff <- function(r, p) {
    p * (crra_u(A_high, r) - crra_u(B_high, r)) +
      (1 - p) * (crra_u(A_low, r) - crra_u(B_low, r))
  }
  r_indiff <- function(p) uniroot(eu_diff, interval = c(-10, 10), p = p, tol = 1e-10)$root
  r_first  <- r_indiff(p_sched[1])        # safe lottery preferred above this r
  r_last   <- r_indiff(p_sched[K_rows])   # risky lottery preferred below this r
  
  f_hl <- file.path(path_src, "holt-laury.csv")
  stopifnot(file.exists(f_hl))
  hl <- fread(f_hl, encoding = "UTF-8")
  choice_cols <- grep("^c\\d+$", names(hl), value = TRUE)
  choice_cols <- choice_cols[order(as.integer(sub("^c", "", choice_cols)))]
  stopifnot(length(choice_cols) == K_rows)
  ch <- matrix(toupper(trimws(as.vector(as.matrix(hl[, ..choice_cols])))), nrow = nrow(hl))
  stopifnot(all(ch %in% c(label_safe, label_risky)))
  safe <- ch == label_safe
  
  patterns <- data.table(
    pid         = as.character(hl$pid),
    n_safe      = rowSums(safe),
    n_switch    = rowSums(safe[, -1L, drop = FALSE] != safe[, -K_rows, drop = FALSE]),
    risky_first = !safe[, 1L],
    safe_last   = safe[, K_rows]
  )
  
  scored_all <- rbindlist(lapply(treatments, function(tr) {
    sc <- fread(model_path("mpl", paste0("scored-", tr, ".csv")), encoding = "UTF-8")
    sc[, .(pid = as.character(pid), treat = tr, r_mean,
           inconsistent = as.integer(inconsistent))]
  }))
  pat <- merge(scored_all, patterns, by = "pid")
  if (nrow(pat) != nrow(scored_all)) stop("Participants of the scored files are missing in holt-laury.csv.")
  
  summarise_patterns <- function(d, tr, grp) {
    if (nrow(d) == 0L) return(NULL)
    data.table(
      treatment = tr, group = grp, n = nrow(d),
      mean_r_bayes          = mean(d$r_mean),
      mean_safe_choices     = mean(d$n_safe),
      mean_switches         = mean(d$n_switch),
      share_risky_first_row = mean(d$risky_first),
      share_safe_last_row   = mean(d$safe_last),
      r_indiff_first_row    = r_first,
      r_indiff_last_row     = r_last
    )
  }
  res_pat <- rbindlist(lapply(treatments, function(tr) {
    d <- pat[treat == tr]
    rbind(summarise_patterns(d,                    tr, "all"),
          summarise_patterns(d[inconsistent == 0L], tr, "consistent"),
          summarise_patterns(d[inconsistent == 1L], tr, "inconsistent"))
  }))
  
  f_pat <- file.path(path_csv, "r-choice-patterns-by-consistency.csv")
  fwrite(res_pat, f_pat)
  msg("Saved: ", f_pat)
  
  # ---- Choice precision by consistency ----
  # The Bayesian MPL model fits erratic choice patterns with a low choice
  # precision lambda_i. For each participant, the posterior mean of
  # log(lambda_i) is computed from the draws of the treatment fit; the table
  # averages it by treatment and group and reports the prior value log(30)
  # for reference. Participants with replaced draws of r_i are left out.
  f_ex <- model_path("mpl", "exact-draw-participants.csv")
  stopifnot(file.exists(f_ex))
  exact_pid <- unique(as.character(fread(f_ex, encoding = "UTF-8")$pid))
  prior_ll  <- log(cfg$model$stan$mpl_bernoulli$prior$lambda_mean)
  
  prec_all <- rbindlist(lapply(treatments, function(tr) {
    fit <- readRDS(model_path("mpl", paste0("fit-", tr, ".rds")))
    pid <- as.character(readRDS(model_path("mpl", paste0("r-draws-", tr, ".rds")))$pid)
    ll  <- colMeans(log(rstan::extract(fit, pars = "lambda")$lambda))
    stopifnot(length(ll) == length(pid))
    data.table(pid = pid, treat = tr, log_lambda = ll)
  }))
  prec_all <- merge(prec_all, scored_all[, .(pid, treat, inconsistent)],
                    by = c("pid", "treat"))
  if (nrow(prec_all) != nrow(scored_all)) stop("Participants of the fits and of the scored files differ.")
  
  summarise_precision <- function(d, tr, grp) {
    keep <- d[!pid %in% exact_pid]
    if (nrow(keep) == 0L) return(NULL)
    data.table(
      treatment = tr, group = grp, n = nrow(keep), n_left_out = nrow(d) - nrow(keep),
      mean_log_lambda  = mean(keep$log_lambda),
      sd_log_lambda    = sd(keep$log_lambda),
      exp_mean_log_lambda = exp(mean(keep$log_lambda)),
      prior_log_lambda = prior_ll
    )
  }
  res_prec <- rbindlist(lapply(treatments, function(tr) {
    d <- prec_all[treat == tr]
    rbind(summarise_precision(d,                     tr, "all"),
          summarise_precision(d[inconsistent == 0L], tr, "consistent"),
          summarise_precision(d[inconsistent == 1L], tr, "inconsistent"))
  }))
  
  f_prec <- file.path(path_csv, "r-choice-precision-by-consistency.csv")
  fwrite(res_prec, f_prec)
  msg("Saved: ", f_prec)
  
  # ---- Comparison of the risk estimates ----
  # Posterior mean of r_i of the primary model against (i) the switching
  # point, (ii) the exact posterior mean of the same model on a grid, computed
  # here for every participant, and (iii) the hierarchical model.
  grid_posterior_mean <- function(ymat, p, prior,
                                  r_grid  = seq(-8, 4, by = 0.005),
                                  ll_grid = seq(-2, 8, by = 0.1)) {
    d1  <- crra_u(A_high, r_grid) - crra_u(B_high, r_grid)
    d2  <- crra_u(A_low,  r_grid) - crra_u(B_low,  r_grid)
    lam <- exp(ll_grid)
    lprior <- outer(dnorm(r_grid,  prior$r_mean, prior$r_sd, log = TRUE),
                    dnorm(ll_grid, log(prior$lambda_mean), prior$lambda_sd, log = TRUE),
                    "+")
    K <- length(p)
    logS <- logR <- vector("list", K)
    for (k in seq_len(K)) {
      z <- outer(p[k] * d1 + (1 - p[k]) * d2, lam)
      logS[[k]] <- -(pmax(-z, 0) + log1p(exp(-abs(z))))   # log logistic(z)
      logR[[k]] <- logS[[k]] - z
    }
    vapply(seq_len(nrow(ymat)), function(i) {
      lp <- lprior
      for (k in seq_len(K)) lp <- lp + (if (ymat[i, k] == 1L) logS[[k]] else logR[[k]])
      w <- rowSums(exp(lp - max(lp)))
      sum(w / sum(w) * r_grid)
    }, numeric(1))
  }
  
  grid_mean <- data.table(
    pid    = as.character(hl$pid),
    r_grid = grid_posterior_mean(matrix(as.integer(safe), nrow = nrow(safe)), p_sched,
                                 cfg$model$stan$mpl_bernoulli$prior)
  )
  
  dir_hier <- file.path(path_mod, "mpl-hierarchical")
  read_hier <- function(tag) {
    f <- file.path(dir_hier, paste0("scored-", gsub("_", "-", tag), ".csv"))
    if (!file.exists(f)) stop("Hierarchical model not found: ", f, call. = FALSE)
    sc <- fread(f, encoding = "UTF-8")[, .(pid = as.character(pid), r_hier = r_mean)]
    # posterior SD of r_i from the draws of the same fit
    f_dr <- file.path(dir_hier, paste0("r-draws-", gsub("_", "-", tag), ".rds"))
    if (!file.exists(f_dr)) stop("Hierarchical draws not found: ", f_dr, call. = FALSE)
    dr <- readRDS(f_dr)
    merge(sc, data.table(pid = as.character(dr$pid), sd_hier = apply(dr$r_draws, 2, sd)), by = "pid")
  }
  
  # Posterior SD of r_i of the primary model (draws of the treatment fit)
  read_primary_sd <- function(tr) {
    dr <- readRDS(model_path("mpl", paste0("r-draws-", tr, ".rds")))
    data.table(pid = as.character(dr$pid), sd_prim = apply(dr$r_draws, 2, sd))
  }
  
  compare_group <- function(tr, grp) {
    d <- scored_all[treat == tr]
    d <- switch(grp, all = d, consistent = d[inconsistent == 0L], inconsistent = d[inconsistent == 1L])
    n <- nrow(d)
    d <- merge(d[, .(pid, r_prim = r_mean)], grid_mean, by = "pid")
    d <- merge(d, read_primary_sd(tr), by = "pid")
    d <- merge(d, read_hier(if (grp == "consistent") paste0(tr, "_consistent") else tr), by = "pid")
    d <- merge(d, dt[fit == tr, .(pid, r_sw)], by = "pid", all.x = TRUE)
    if (nrow(d) != n) stop("Participants of the model comparison differ (", tr, ", ", grp, ").", call. = FALSE)
    cons <- identical(grp, "consistent")
    data.table(
      treatment = tr, group = grp, n = n,
      mean_r_primary      = mean(d$r_prim),
      mean_r_grid         = mean(d$r_grid),
      mean_r_hierarchical = mean(d$r_hier),
      mean_r_switching    = if (cons) mean(d$r_sw) else NA_real_,
      cor_grid            = cor(d$r_prim, d$r_grid),
      cor_hierarchical    = cor(d$r_prim, d$r_hier),
      cor_switching       = if (cons) cor(d$r_prim, d$r_sw) else NA_real_,
      max_abs_diff_grid   = max(abs(d$r_prim - d$r_grid)),
      mean_abs_diff_hierarchical = mean(abs(d$r_hier - d$r_prim)),
      mean_post_sd_primary      = mean(d$sd_prim),
      mean_post_sd_hierarchical = mean(d$sd_hier)
    )
  }
  res_cmp <- rbindlist(lapply(treatments, function(tr) {
    rbindlist(lapply(c("all", "consistent", "inconsistent"), function(g) compare_group(tr, g)))
  }))
  
  f_cmp <- file.path(path_csv, "r-model-comparison.csv")
  fwrite(res_cmp, f_cmp)
  msg("Saved: ", f_cmp)
  
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
