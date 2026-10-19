# ============================================================
# _rq2-ppc.R
#
# PURPOSE
#   Posterior predictive checks of all-in stakes for the one-inflated RQ2
#   models (22_rq2-stan-one-inflated.R, 23_rq2-stan-one-inflated-seq.R).
#   Three statistics are computed from replicated stakes generated with the
#   fitted participant, sequence and drift effects:
#     allin          share of betting trials with an all-in stake (stake = e);
#                    preregistered check, one-sided (excess of all-in stakes)
#     pure_contrast  all-in share in pure sequences (one symbol repeated)
#                    minus the all-in share in all other sequences;
#                    not preregistered, two-sided
#     seq_sd         SD across sequences of the sequence-level all-in share;
#                    not preregistered, one-sided (excess dispersion)
#   The preregistered statistic is reproduced by construction by a model
#   with a free all-in probability per participant. The two other statistics
#   test whether the all-in probability varies with the sequence as in the
#   data, which a model without a sequence term in the all-in probability
#   cannot generate beyond chance.
#
# INPUT (arguments)
#   fit, d (columns seq, stake, block_c), pid_levels, seq_levels,
#   include_drift, e, K_ppc, seed, cut, tr, sample, run (output folder),
#   model_label (figure title); seq_effect = TRUE if the fit contains the
#   sequence effects c_s of the all-in probability; figure = FALSE skips the
#   figure (run is then not used)
#
# OUTPUT
#   data.table with one row per statistic: observed value, replicate median
#   and 95% interval, p-value, threshold, pass; and the figure
#   <run>/ppc-<tr>-<sample>.png (one panel per statistic)
#
# NOTES
#   - Pure sequences are identified from the sequence string (a single
#     repeated symbol); a sample without a pure sequence skips pure_contrast
#   - p-values: P(replicate >= observed) for one-sided statistics; twice the
#     smaller tail probability for the two-sided statistic. A check passes
#     if p >= 1 - cut (compared after rounding to 10 decimals, so that a
#     p-value equal to the cut passes)
# ============================================================

rq2_ppc_stats <- function(fit, d, pid_levels, seq_levels, include_drift, e,
                          K_ppc, seed, cut, tr, sample, run, model_label,
                          seq_effect = FALSE, figure = TRUE) {

  pid_i <- match(d$pid, pid_levels)
  sid_s <- match(d$seq, seq_levels)
  stopifnot(!anyNA(pid_i), !anyNA(sid_s))
  Tn <- nrow(d)
  S  <- length(seq_levels)

  pars <- c("alpha", "u", "b", "v", "phi", "alpha_omega", "gamma_drift",
            if (seq_effect) "c")
  post  <- rstan::extract(fit, pars = pars)
  K_all <- length(post$alpha)
  K_use <- min(K_ppc, K_all)
  set.seed(seed + 1001L)
  k_idx <- sort(sample.int(K_all, K_use))

  is_pure  <- grepl("^(.)\\1*$", seq_levels)
  trial_pure <- is_pure[sid_s]
  n_s <- tabulate(sid_s, nbins = S)
  use_pure <- any(trial_pure) && any(!trial_pure)

  stats_of <- function(allin) {
    share_s <- as.vector(rowsum(as.numeric(allin), sid_s, reorder = TRUE)) / n_s
    c(allin = mean(allin),
      pure_contrast = if (use_pure) mean(allin[trial_pure]) - mean(allin[!trial_pure]) else NA_real_,
      seq_sd = sd(share_s))
  }

  obs <- stats_of(d$stake >= e)

  set.seed(seed + 2000L)
  rep <- t(vapply(k_idx, function(k) {
    eta <- post$alpha[k] + post$u[k, pid_i] + post$b[k, sid_s]
    if (include_drift) eta <- eta + post$gamma_drift[k] * d$block_c
    mu <- plogis(eta)
    lo <- post$alpha_omega[k] + post$v[k, pid_i]
    if (seq_effect) lo <- lo + post$c[k, sid_s]
    omega <- plogis(lo)
    p     <- rbeta(Tn, pmax(mu * post$phi[k], 1e-8), pmax((1 - mu) * post$phi[k], 1e-8))
    stake <- 1L + rbinom(Tn, e - 1L, p)
    stake[runif(Tn) < omega] <- e
    stats_of(stake >= e)
  }, numeric(3L)))
  colnames(rep) <- names(obs)

  spec <- data.table(
    statistic     = c("allin", "pure_contrast", "seq_sd"),
    label         = c("Share of bets that are all-in stakes",
                      "All-in share, pure minus other sequences",
                      "SD across sequences of the all-in share"),
    preregistered = c(TRUE, FALSE, FALSE),
    two_sided     = c(FALSE, TRUE, FALSE)
  )
  if (!use_pure) spec <- spec[statistic != "pure_contrast"]

  rows <- rbindlist(lapply(seq_len(nrow(spec)), function(j) {
    s <- spec$statistic[j]
    r <- rep[, s]
    p_up <- mean(r >= obs[[s]])
    p_lo <- mean(r <= obs[[s]])
    p    <- if (spec$two_sided[j]) min(1, 2 * min(p_up, p_lo)) else p_up
    data.table(
      treatment = tr, sample = sample, statistic = s, label = spec$label[j],
      preregistered = spec$preregistered[j], observed = obs[[s]],
      rep_median = median(r), rep_q025 = unname(quantile(r, 0.025)),
      rep_q975 = unname(quantile(r, 0.975)),
      ratio = obs[[s]] / median(r), n_rep = K_use, p_value = p,
      threshold = 1 - cut, pass = round(p, 10) >= round(1 - cut, 10)
    )
  }))

  if (figure) {
    # ---- Figure: replicated distribution and observed value per statistic ----
    long <- rbindlist(lapply(spec$statistic, function(s)
      data.table(statistic = s, v = rep[, s])))
    long[, panel := factor(spec$label[match(statistic, spec$statistic)], levels = spec$label)]
    ann <- rows[, .(statistic, observed,
                    panel = factor(spec$label[match(statistic, spec$statistic)], levels = spec$label),
                    txt = ifelse(p_value < 1 / K_use, sprintf("p < %.3f", 1 / K_use),
                                 sprintf("p = %.3f", p_value)))]
    sample_lab <- paste0(c(m25 = "FN", m19 = "FP")[[tr]], " ", sample)
    g <- ggplot(long, aes(x = v)) +
      geom_density(fill = "#74ADD1", alpha = 0.5, colour = "#2166AC", linewidth = 0.7) +
      geom_vline(data = ann, aes(xintercept = observed), colour = "#C0392B",
                 linewidth = 0.8, linetype = "dashed") +
      geom_text(data = ann, aes(label = txt), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
                size = 3, colour = "grey30", inherit.aes = FALSE) +
      facet_wrap(~panel, scales = "free", nrow = 1) +
      labs(title = paste0(model_label, ", posterior predictive checks: ", sample_lab),
           x = NULL, y = "Density") +
      theme_classic(base_size = 10) +
      theme(strip.background = element_blank(), strip.text = element_text(size = 8.5),
            axis.title.y = element_text(size = 8.5),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
            plot.title = element_text(size = 11), plot.margin = margin(4, 8, 4, 4))
    ggsave(file.path(run, paste0("ppc-", tr, "-", sample, ".png")), g,
           width = 3.2 * nrow(spec) + 0.6, height = 3.6, dpi = 300)
  }

  rows
}
