# ============================================================
# 95_synthesis_stan_outputs.R
#
# PURPOSE
#   Everything downstream of the synthesis joint model fit produced
#   by 94_synthesis_stan.R (synthesis_fit_<tr>.rds): convergence
#   diagnostics + PPC, the fixed-effect/dissociation/within-corr
#   tables, per-participant heterogeneity classification (for the
#   debrief comparison), population-averaged marginal effects on the
#   probability scale, draw-level Cohen's h for the run and imbalance
#   cues, and the slope-distribution figure -- all in one function,
#   one pass over the fit, one draws extraction.
#
#   Assumes the additive (no-interaction) specification: beta_side
#   has 3 parameters (see synthesis.stan).
#
# INPUT
#   path_mod/synthesis_fit_<tr>.rds  (from 94_synthesis_stan.R)
#
# OUTPUT
#   All output files from this script share the "synthesis_stan_" prefix
#   (distinct from the "synthesis_fit_<tr>.rds" input above, which comes
#   from 94_synthesis_stan.R).
#   path_out/synthesis_stan_diagnostics_<tr>.csv             -- convergence gate + PPC (nothing is
#                                                                 printed; check passed_gate here)
#   path_out/synthesis_stan_fixed_effects_<tr>.csv
#   path_out/synthesis_stan_dissociation_<tr>.csv            -- the formal dissociation test
#   path_out/synthesis_stan_dissociation_extra_<tr>.csv      -- NEW: magnitude-based dissociation
#                                                                 rows + participant-level
#                                                                 selection checks (see the
#                                                                 EXTRA block below)
#   path_out/synthesis_stan_within_corr_<tr>.csv              -- within-equation RE correlations
#   path_out/synthesis_stan_tau_<tr>.csv                      -- participant-level RE SDs
#                                                                 (heterogeneity magnitude, one
#                                                                 row per equation x regressor)
#   path_out/synthesis_stan_ame_<tr>.csv                       -- population-averaged marginal
#                                                                 effects, probability scale,
#                                                                 plus RE-integrated baseline
#                                                                 probabilities (last 2 rows)
#   path_out/synthesis_stan_cohens_h_<tr>.csv                  -- NEW: draw-level Cohen's h (with
#                                                                 95% CrI) for the terminal-run and
#                                                                 imbalance cues, betting and side
#                                                                 choice, plus the two counterfactual
#                                                                 probabilities behind each h (see
#                                                                 the COHEN'S H block below)
#   path_out/synthesis_stan_run_slope_participants_<tr>.csv   -- per-participant side run-slope
#   path_out/synthesis_stan_imb_slope_participants_<tr>.csv   -- per-participant side imb-slope
#                                                                 (for the debrief comparison)
#   path_out/synthesis_stan_bet_run_slope_participants_<tr>.csv -- per-participant bet run-slope
#   path_out/synthesis_stan_bet_imb_slope_participants_<tr>.csv -- per-participant bet imb-slope
#   path_out/synthesis_stan_class_counts_<tr>.csv              -- counts of participants in each
#                                                                 class bucket, for all four random
#                                                                 slopes (bet-run, bet-imb via
#                                                                 engaged/disengaged/uncertain;
#                                                                 side-run via with_run/against_run/
#                                                                 uncertain; side-imb via with_imb/
#                                                                 against_imb/uncertain -- "with_imb"
#                                                                 = picks the majority side, extending
#                                                                 the imbalance; "against_imb" = picks
#                                                                 the minority side, countering it)
#   path_fig/synthesis_stan_slope_distributions.png            -- two-panel run/imb slope figure,
#                                                                 m25 and m19 density curves overlaid
#                                                                 per panel (ggplot2 + patchwork)
#
# CONSOLE OUTPUT
#   Only "Saved: <path>" lines. Diagnostics go to synthesis_stan_diagnostics_<tr>.csv.
#
# NOTE ON SIGN CODING
#   Betting-equation slopes are UNSIGNED (bet more when a run is present /
#   the sequence is imbalanced). Side-equation slopes are SIGNED toward
#   Heads: the run slope is negative for betting against the run, the
#   imbalance slope positive for going with the majority. The side
#   intercept is a baseline preference for the "Heads" label. Keep this in
#   mind when reading cross-equation correlations -- see the EXTRA block.
#
# ============================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------
# DIAGNOSTICS + TABLES + HETEROGENEITY + AME + COHEN'S H + FIGURE, one pass per treatment
# ------------------------------------------------------------
synthesis_stan_outputs <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  gate <- cfg$model$diagnostics$synthesis %||%
    list(rhat_max = 1.01, divergent_max = 0L, ess_bulk_min = 400L)
  
  col_bar <- "#2a78d6"
  col_pop <- "#e34948"
  
  qsumm <- function(x) c(median = median(x), q025 = quantile(x, .025),
                         q975 = quantile(x, .975), p_gt0 = mean(x > 0))
  
  # Cohen's h for two proportions: difference of arcsine-transformed
  # probabilities, h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p0)). Signed: positive
  # when p1 > p0. Conventional benchmarks are |h| = 0.2 small, 0.5 medium,
  # 0.8 large.
  cohen_h <- function(p1, p0) 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p0))
  
  # only console output of this script: the path of each saved file
  save_csv <- function(x, stem) {
    f <- file.path(path_out, paste0(stem, tr, ".csv"))
    fwrite(x, f)
    msg("Saved: ", f)
  }
  
  out <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(path_mod, paste0("synthesis_fit_", tr, ".rds"))
    if (!file.exists(f_fit)) next
    obj <- readRDS(f_fit)
    fit <- obj$fit
    dl  <- obj$data_list
    pid_levels <- obj$pid_levels
    
    # ---- single draws extraction, reused by diagnostics/tables/heterogeneity/AME ----
    draws <- rstan::extract(fit, pars = c(
      "beta_bet", "beta_side", "tau", "sigma_v_bet", "sigma_v_side",
      "corr_intercepts", "corr_run_slopes", "corr_imb_slopes", "Omega",
      "u", "v_bet", "v_side"
    ))
    
    # ============================================================
    # DIAGNOSTICS
    # ============================================================
    summ <- summary(fit)$summary
    summ <- summ[!grepl("^(u\\[|z\\[|L_Omega\\[|log_lik|.*_rep\\[|lp__)", rownames(summ)), , drop = FALSE]
    
    rhat_bad <- rownames(summ)[which(summ[, "Rhat"] > gate$rhat_max)]
    ess_bad  <- rownames(summ)[which(summ[, "n_eff"] < gate$ess_bulk_min)]
    
    sp <- get_sampler_params(fit, inc_warmup = FALSE)
    divergences <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
    
    ok <- divergences <= gate$divergent_max && length(rhat_bad) == 0 && length(ess_bad) == 0
    # Convergence results are written to synthesis_stan_diagnostics_<tr>.csv
    # (below, together with the PPC) instead of being printed; check the
    # `passed_gate` row before using the outputs. If it is FALSE, see the
    # NOTE ON SAMPLER SETTINGS header in 94_synthesis_stan.R (tau[6] and the
    # two slope correlations are the usual suspects).
    
    # tau diagnostics: participant-level RE SDs, one per equation x regressor.
    # Printed explicitly (not just gated on Rhat/n_eff) because a near-zero tau
    # on a given slope is itself diagnostic: it means that random slope has
    # little real participant-level variance, which is the usual reason the
    # corresponding cross-equation correlation (corr_run_slopes / corr_imb_slopes)
    # mixes poorly or fails to converge -- more iterations won't fix that, a
    # small tau there is telling you the correlation is weakly identified.
    #
    # tau IS ALSO the direct answer to "how much individual heterogeneity is
    # there in run-sensitivity / imbalance-sensitivity, for betting and for
    # side choice": tau[2]/tau[3] are the bet-equation run/imbalance
    # participant-level SDs, tau[5]/tau[6] are the same for the side equation.
    # A large tau relative to the population fixed effect (beta_bet/beta_side)
    # means participants differ a lot in whether/how strongly they use that
    # cue -- exactly the quantity needed to pivot into a heterogeneity
    # analysis (e.g., what fraction of participants have a *credibly*
    # positive vs. negative personal slope; see the classification counts
    # below). Saved to CSV here.
    tau_labels <- c("bet_intercept_sd", "bet_run_sd", "bet_imb_sd",
                    "side_intercept_sd", "side_run_sd", "side_imb_sd")
    tau_med <- apply(draws$tau, 2, median)
    tau_lo  <- apply(draws$tau, 2, quantile, probs = .025)
    tau_hi  <- apply(draws$tau, 2, quantile, probs = .975)
    tau_rhat <- summ[paste0("tau[", 1:6, "]"), "Rhat"]
    
    # fixed-effect magnitude for the same six rows, so the tau table can be
    # read alongside "how big is the population average effect this
    # heterogeneity sits on top of" without a second join -- NA for the two
    # intercept rows (no single comparable "effect" for an intercept SD).
    fe_for_tau <- c(NA_real_,
                    median(draws$beta_bet[, 2]), median(draws$beta_bet[, 3]),
                    NA_real_,
                    median(draws$beta_side[, 2]), median(draws$beta_side[, 3]))
    
    tau_tbl <- data.table(
      param             = tau_labels,
      median            = tau_med,
      q025              = tau_lo,
      q975              = tau_hi,
      rhat              = tau_rhat,
      fixed_effect_median = fe_for_tau,
      heterogeneity_ratio  = abs(tau_med / fe_for_tau)   # tau relative to |population effect|
    )
    save_csv(tau_tbl, "synthesis_stan_tau_")
    
    # posterior predictive check: observed vs. replicated bet rate and
    # H-share. Generated here in R on a subsample of draws rather than in
    # Stan at full draw count -- bet_rep/side_rep are per-trial x per-draw
    # arrays that at production draw counts (tens of thousands across 4
    # chains) blow past typical machine memory if Stan materializes one
    # replicate per retained draw. PPC only ever needs a few hundred.
    #
    # The same per-draw loop also builds the population-averaged marginal
    # effects (AME) below, and the draw-level Cohen's h: all need
    # eta_bet/eta_side computed at the *actual* data (real pid/sid, real
    # covariates, real participant- and sequence-level random effects) for
    # a subsample of draws, so it is computed once here rather than in a
    # second pass over draws.
    # Fixed seed, reset per treatment: the draw subsample (idx) and the
    # rbinom() PPC replicates below are then identical on every re-run and
    # independent of treatment order, so the AME / PPC / Cohen's h numbers
    # do not drift.
    set.seed(as.integer(cfg$run$seed %||% 1L))
    k <- 1000L
    n_draws <- nrow(draws$beta_bet)
    idx <- sample(seq_len(n_draws), min(k, n_draws))
    
    rep_bet_rate  <- numeric(length(idx))
    rep_side_rate <- numeric(length(idx))
    
    # AME accumulators. All four use the SAME definition -- instantaneous
    # dp/dx = coefficient * p * (1-p), evaluated at each observation's
    # actual covariate profile and that draw's actual participant/sequence
    # random effects, then averaged over the N / N_bet observations -- so
    # "per unit" means the same thing in every row of the table. (Earlier
    # versions used a discrete 0->1 / full -1->+1 swing for the two run
    # regressors, which mixed a "per unit" convention for the imbalance
    # cues with a "discrete change," sometimes a 2-unit change, for the
    # run cues -- inconsistent units across the table. This version:
    #   ame_bet_run  : dp/dx of run_dum (0-1 scale)
    #   ame_bet_imb  : dp/dx of imb_mag (0-1 scale)
    #   ame_side_run : dp/dx of run_h_s (-1 to 1 scale)
    #   ame_side_imb : dp/dx of imb_h_s (-1 to 1 scale)
    # NOT evaluated at random effects = 0, which is what a naive plogis(beta)
    # conversion would silently assume.
    ame_bet_run  <- numeric(length(idx))
    ame_bet_imb  <- numeric(length(idx))
    ame_side_run <- numeric(length(idx))
    ame_side_imb <- numeric(length(idx))
    
    # Baseline predicted probability (no terminal run, zero imbalance) for
    # each equation's Intercept row. NOT plogis(fixed-effect intercept) --
    # that reference-point number ignores intercept heterogeneity entirely,
    # and bet_intercept_sd/side_intercept_sd (see synthesis_stan_tau_<tr>.csv)
    # show that heterogeneity is NOT negligible, especially for the betting
    # equation (bet_intercept_sd is the single largest tau in the model).
    # Because plogis() is nonlinear, plogis(mean(eta)) != mean(plogis(eta))
    # once the random-effect spread is large -- the naive version
    # systematically overstates the population-average baseline probability.
    # Computed here by averaging plogis() over the full P x S grid of that
    # draw's actual participant- and sequence-level intercept random
    # effects (independent per the model), at run_dum = 0 / imb_mag = 0
    # (side equation: run_h_s = 0 / imb_h_s = 0).
    baseline_bet  <- numeric(length(idx))
    baseline_side <- numeric(length(idx))
    
    # COHEN'S H accumulators (one row per draw, one column per quantity).
    # Each h is a draw-level COUNTERFACTUAL contrast: for that draw, the
    # cue is switched from its reference level (p0) to its treated level
    # (p1) on the actual trials, holding everything else (the other cue,
    # the participant and sequence random effects) at its actual value,
    # the probabilities are averaged over those trials, and then
    # h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p0)). This is the exact
    # discrete-change analogue of the AME above (the AME is the
    # instantaneous dp/dx), and the CrI comes from the posterior draws, so
    # no reference-point or interval-endpoint approximation is involved.
    # Definitions (all signed; the sign follows the probability that is
    # modelled):
    #   bet_run             : P(bet), run_dum 1 vs 0, all trials
    #   bet_imb             : P(bet), imb_mag 1 vs 0 (fully imbalanced vs.
    #                         balanced), all trials
    #   side_run            : P(choose the RUN's symbol), actual run vs.
    #                         run_h_s set to 0, run-ending trials only.
    #                         Negative = betting against the run.
    #   side_run_vs_chance  : same p1, but p0 = 0.5 (H and T are symmetric
    #                         by design), the natural benchmark for the
    #                         out-of-sample accuracy in the text
    #   side_imb_observed   : P(choose the MAJORITY symbol), observed
    #                         imbalance vs. imb_h_s set to 0, imbalanced
    #                         trials only. Positive = hot-handish.
    #   side_imb_full       : same, but with |imb_h_s| = 1 (0-6 sequences)
    #                         instead of the observed mix of 2-4 / 1-5 /
    #                         0-6, i.e. the counterpart of the "balanced
    #                         to fully imbalanced" contrast in the text
    #   side_imb_*_run      : the two imbalance rows above, restricted to
    #   side_imb_*_norun      imbalanced trials that END IN a terminal run
    #                         (_run) or DO NOT (_norun). The model is
    #                         additive (no run x imbalance interaction), so
    #                         the imbalance slope on the logit scale is the
    #                         same in both strata; the two h values differ
    #                         only through the trial composition and the
    #                         baseline probability the run cue shifts
    #                         (near 0.5 without a run, pushed away from the
    #                         run's symbol with one). They are a
    #                         model-implied split, not a test of an
    #                         interaction.
    ch_rows <- c("bet_run", "bet_imb", "side_run", "side_run_vs_chance",
                 "side_imb_observed", "side_imb_full",
                 "side_imb_observed_run", "side_imb_observed_norun",
                 "side_imb_full_run", "side_imb_full_norun")
    ch_def  <- c("P(bet): run present vs. absent, all trials",
                 "P(bet): fully imbalanced (imb_mag = 1) vs. balanced, all trials",
                 "P(choose run's symbol): run vs. run cue set to 0, run-ending trials",
                 "P(choose run's symbol): run vs. chance (0.5), run-ending trials",
                 "P(choose majority symbol): observed imbalance vs. balanced, imbalanced trials",
                 "P(choose majority symbol): |imbalance| = 1 vs. balanced, imbalanced trials",
                 "P(choose majority symbol): observed imbalance vs. balanced, imbalanced trials WITH a terminal run",
                 "P(choose majority symbol): observed imbalance vs. balanced, imbalanced trials WITHOUT a terminal run",
                 "P(choose majority symbol): |imbalance| = 1 vs. balanced, imbalanced trials WITH a terminal run",
                 "P(choose majority symbol): |imbalance| = 1 vs. balanced, imbalanced trials WITHOUT a terminal run")
    ch_h  <- matrix(NA_real_, nrow = length(idx), ncol = length(ch_rows),
                    dimnames = list(NULL, ch_rows))
    ch_p1 <- ch_h
    ch_p0 <- ch_h
    # store one counterfactual pair for a draw (h + the two probabilities)
    ch_rec <- function(j, nm, p1, p0) {
      ch_h[j, nm]  <<- cohen_h(p1, p0)
      ch_p1[j, nm] <<- p1
      ch_p0[j, nm] <<- p0
    }
    
    for (j in seq_along(idx)) {
      d <- idx[j]
      
      run_coef_bet <- draws$beta_bet[d, 2] + draws$u[d, 2, dl$pid]
      imb_coef_bet <- draws$beta_bet[d, 3] + draws$u[d, 3, dl$pid]
      re_bet <- draws$beta_bet[d, 1] + draws$u[d, 1, dl$pid] + draws$v_bet[d, dl$sid]
      eta_bet <- re_bet + run_coef_bet * dl$run_dum + imb_coef_bet * dl$imb_mag
      p_bet <- plogis(eta_bet)
      rep_bet_rate[j] <- mean(rbinom(dl$N, 1, p_bet))
      
      # both AMEs: instantaneous dp/dx at each trial's actual covariate
      # profile and actual random effects -- same formula, same "per unit"
      # convention, for both regressors
      ame_bet_run[j] <- mean(run_coef_bet * p_bet * (1 - p_bet))
      ame_bet_imb[j] <- mean(imb_coef_bet * p_bet * (1 - p_bet))
      
      # Cohen's h, betting equation: switch one cue, keep the other actual
      ch_rec(j, "bet_run",
             mean(plogis(re_bet + run_coef_bet * 1 + imb_coef_bet * dl$imb_mag)),
             mean(plogis(re_bet + imb_coef_bet * dl$imb_mag)))
      ch_rec(j, "bet_imb",
             mean(plogis(re_bet + run_coef_bet * dl$run_dum + imb_coef_bet * 1)),
             mean(plogis(re_bet + run_coef_bet * dl$run_dum)))
      
      run_coef_side <- draws$beta_side[d, 2] + draws$u[d, 5, dl$pid_s]
      imb_coef_side <- draws$beta_side[d, 3] + draws$u[d, 6, dl$pid_s]
      re_side <- draws$beta_side[d, 1] + draws$u[d, 4, dl$pid_s] + draws$v_side[d, dl$sid_s]
      eta_side <- re_side + run_coef_side * dl$run_h_s + imb_coef_side * dl$imb_h_s
      p_side <- plogis(eta_side)
      rep_side_rate[j] <- mean(rbinom(dl$N_bet, 1, p_side))
      
      # both AMEs: same dp/dx formula and "per unit" convention as the bet
      # equation above
      ame_side_run[j] <- mean(run_coef_side * p_side * (1 - p_side))
      ame_side_imb[j] <- mean(imb_coef_side * p_side * (1 - p_side))
      
      # Cohen's h, side equation. The cues are signed toward Heads, so each
      # trial's probability is re-expressed as the probability of choosing
      # the RUN's symbol (run trials) or the MAJORITY symbol (imbalanced
      # trials): plogis(eta) if the cue points to Heads, 1 - plogis(eta)
      # if it points to Tails.
      k_r <- dl$run_h_s != 0
      if (any(k_r)) {
        sg_r  <- dl$run_h_s[k_r]
        base_r <- re_side[k_r] + imb_coef_side[k_r] * dl$imb_h_s[k_r]
        p_run_sym <- function(eta) ifelse(sg_r > 0, plogis(eta), 1 - plogis(eta))
        q1 <- mean(p_run_sym(base_r + run_coef_side[k_r] * sg_r))   # actual run
        q0 <- mean(p_run_sym(base_r))                               # run cue set to 0
        ch_rec(j, "side_run", q1, q0)
        ch_rec(j, "side_run_vs_chance", q1, 0.5)
      }
      
      k_i <- dl$imb_h_s != 0
      if (any(k_i)) {
        sg_i  <- sign(dl$imb_h_s[k_i])
        base_i <- re_side[k_i] + run_coef_side[k_i] * dl$run_h_s[k_i]
        p_maj_sym <- function(eta) ifelse(sg_i > 0, plogis(eta), 1 - plogis(eta))
        q0i <- mean(p_maj_sym(base_i))                                                     # balanced
        ch_rec(j, "side_imb_observed",
               mean(p_maj_sym(base_i + imb_coef_side[k_i] * dl$imb_h_s[k_i])), q0i)        # observed mix
        ch_rec(j, "side_imb_full",
               mean(p_maj_sym(base_i + imb_coef_side[k_i] * sg_i)), q0i)                   # |imb| = 1
        
        # same contrasts, split by whether the trial's sequence ends in a
        # terminal run (stratum masks are relative to the imbalanced trials)
        has_run_i <- dl$run_h_s[k_i] != 0
        eta_obs_i  <- base_i + imb_coef_side[k_i] * dl$imb_h_s[k_i]
        eta_full_i <- base_i + imb_coef_side[k_i] * sg_i
        for (st in c("run", "norun")) {
          m <- if (st == "run") has_run_i else !has_run_i
          if (!any(m)) next
          pm_st <- function(eta) mean(ifelse(sg_i[m] > 0, plogis(eta[m]), 1 - plogis(eta[m])))
          q0_st <- pm_st(base_i)
          ch_rec(j, paste0("side_imb_observed_", st), pm_st(eta_obs_i),  q0_st)
          ch_rec(j, paste0("side_imb_full_", st),     pm_st(eta_full_i), q0_st)
        }
      }
      
      # baseline probabilities: joint P x S grid of intercept random effects,
      # fixed effect intercept only (reference covariates = 0)
      baseline_bet[j]  <- mean(plogis(outer(draws$u[d, 1, ], draws$v_bet[d, ],  "+") + draws$beta_bet[d, 1]))
      baseline_side[j] <- mean(plogis(outer(draws$u[d, 4, ], draws$v_side[d, ], "+") + draws$beta_side[d, 1]))
    }
    
    obs_bet_rate  <- mean(dl$bet)
    obs_side_rate <- mean(dl$side)
    
    diag_tbl <- data.table(
      check = c("divergences", "rhat_above_max", "n_eff_below_min", "passed_gate",
                "ppc_bet_rate_observed", "ppc_bet_rate_rep_q025", "ppc_bet_rate_rep_q975",
                "ppc_heads_share_observed", "ppc_heads_share_rep_q025", "ppc_heads_share_rep_q975"),
      value = c(as.character(divergences),
                if (length(rhat_bad)) paste0(rhat_bad, "=", sprintf("%.3f", summ[rhat_bad, "Rhat"]), collapse = "; ") else "none",
                if (length(ess_bad))  paste0(ess_bad,  "=", round(summ[ess_bad, "n_eff"]),          collapse = "; ") else "none",
                as.character(ok),
                sprintf("%.4f", c(obs_bet_rate, quantile(rep_bet_rate, c(.025, .975)),
                                  obs_side_rate, quantile(rep_side_rate, c(.025, .975)))))
    )
    save_csv(diag_tbl, "synthesis_stan_diagnostics_")
    
    # ============================================================
    # AVERAGE MARGINAL EFFECTS (probability scale)
    # ============================================================
    # Population-averaged AMEs, correctly integrating over the empirical
    # distribution of participants, sequences, and covariates (via the
    # actual N/N_bet trial-level data) and over posterior uncertainty (via
    # the k subsampled draws above) -- NOT a reference-point conversion at
    # random effects = 0. See ame_* definitions in the comment above the
    # loop. p_gt0 here means "share of subsampled draws with a positive
    # average marginal effect," directly comparable in spirit to the p_gt0
    # column in the other tables in this file.
    ame_tbl <- data.table(
      quantity = c("bet_run (dp/dx per unit of run_dum, 0-1 scale -- run vs. no run)",
                   "bet_imb (dp/dx per unit of imbalance magnitude, 0-1 scale -- full range)",
                   "side_run (dp/dx per unit of signed run cue, -1 to 1 scale -- neutral to one direction)",
                   "side_imb (dp/dx per unit of signed imbalance, -1 to 1 scale -- neutral to one direction)",
                   "bet_baseline (P(bet), no run, zero imbalance, RE-integrated)",
                   "side_baseline (P(heads), no run, zero imbalance, RE-integrated)"),
      t(sapply(list(ame_bet_run, ame_bet_imb, ame_side_run, ame_side_imb,
                    baseline_bet, baseline_side), qsumm))
    )
    save_csv(ame_tbl, "synthesis_stan_ame_")
    
    # ============================================================
    # COHEN'S H (draw-level, effect size for the run and imbalance cues)
    # ============================================================
    # One row per contrast defined above the loop. Columns: the posterior
    # median / 95% CrI of the signed h, p_gt0 (share of draws with h > 0),
    # and the posterior medians of the two counterfactual probabilities
    # (p1 = cue at its treated level, p0 = reference), so the h can be
    # read next to the percentage-point difference p1 - p0 it summarises.
    # For the run rows h is NEGATIVE when participants bet against the
    # run (probability of choosing the run's symbol falls); report |h|
    # for the size of the effect. Benchmarks: 0.2 small, 0.5 medium,
    # 0.8 large.
    cohen_tbl <- data.table(
      quantity   = ch_rows,
      definition = ch_def,
      t(sapply(ch_rows, function(nm)
        c(qsumm(ch_h[, nm]),
          p1_median = median(ch_p1[, nm]),
          p0_median = median(ch_p0[, nm]),
          diff_pp_median = 100 * median(ch_p1[, nm] - ch_p0[, nm]))))
    )
    save_csv(cohen_tbl, "synthesis_stan_cohens_h_")
    
    # ============================================================
    # TABLES -- fixed effects + the dissociation-test correlations
    # ============================================================
    fixed_tbl <- rbind(
      data.table(param = c("bet_intercept", "bet_run_slope", "bet_imb_slope"),
                 t(sapply(1:3, function(j) qsumm(draws$beta_bet[, j])))),
      data.table(param = c("side_intercept", "side_run_slope", "side_imb_slope"),
                 t(sapply(1:3, function(j) qsumm(draws$beta_side[, j]))))
    )
    
    dissoc_tbl <- data.table(
      correlation = c("bet_intercept <-> side_intercept",
                      "bet_run_slope <-> side_run_slope",
                      "bet_imb_slope <-> side_imb_slope"),
      t(sapply(list(draws$corr_intercepts, draws$corr_run_slopes, draws$corr_imb_slopes), qsumm))
    )
    
    # within-equation correlations (Omega[1,2]/[1,3]/[2,3] for bet;
    # [4,5]/[4,6]/[5,6] for side) -- needed for the exact (non-diagonal)
    # sequence-level participant-disagreement calculation: Var(eta_s) =
    # x' Sigma x where Sigma = diag(tau) %*% Omega_block %*% diag(tau) and
    # x = c(1, cue1_s, cue2_s). The dissociation table above only carries
    # the three CROSS-equation entries; this table carries the six
    # WITHIN-equation ones the diagonal-only approximation was ignoring.
    within_tbl <- data.table(
      correlation = c("bet_intercept <-> bet_run_slope",
                      "bet_intercept <-> bet_imb_slope",
                      "bet_run_slope <-> bet_imb_slope",
                      "side_intercept <-> side_run_slope",
                      "side_intercept <-> side_imb_slope",
                      "side_run_slope <-> side_imb_slope"),
      t(sapply(list(draws$Omega[, 1, 2], draws$Omega[, 1, 3], draws$Omega[, 2, 3],
                    draws$Omega[, 4, 5], draws$Omega[, 4, 6], draws$Omega[, 5, 6]),
               qsumm))
    )
    
    save_csv(fixed_tbl, "synthesis_stan_fixed_effects_")
    save_csv(dissoc_tbl, "synthesis_stan_dissociation_")
    save_csv(within_tbl, "synthesis_stan_within_corr_")
    
    # ============================================================
    # EXTRA DISSOCIATION + SELECTION CHECKS
    # ============================================================
    # Why: the side-choice slopes are SIGNED (toward Heads / toward the
    # majority), while the betting slopes are unsigned. For the run cue this
    # is harmless -- nearly every participant's side run-slope is negative, so
    # the signed slope is essentially minus its magnitude. For the imbalance
    # cue it is not: participants split into with- and against-majority, so a
    # participant who responds strongly in EITHER direction has a large
    # |side imb-slope|, and a linear correlation between the unsigned betting
    # slope and the SIGNED side slope cannot detect a shared "responsiveness"
    # disposition (a U-shaped link averages to ~0). These rows therefore
    # correlate against the MAGNITUDE of the personal side slopes.
    #
    # The two "SELECTION" rows test the participant-level version of the
    # mechanical explanation in the synthesis footnote: if selective bettors
    # (low betting intercept) are more directional, the betting intercept
    # should correlate POSITIVELY with the (negative) side run-slope and
    # NEGATIVELY with |side imb-slope|.
    #
    # Rows marked "Omega" are the model's population correlation parameter;
    # rows marked "derived" are computed per draw across participants from the
    # personal effects (fixed effect + random effect), which are shrunk, so
    # they are conservative (attenuated toward zero).
    rowcor <- function(A, B) {                       # per-draw Pearson across participants
      A <- A - rowMeans(A); B <- B - rowMeans(B)
      rowSums(A * B) / sqrt(rowSums(A^2) * rowSums(B^2))
    }
    bet_int_i  <- draws$beta_bet[, 1]  + draws$u[, 1, ]   # [iterations, P]
    bet_run_i  <- draws$beta_bet[, 2]  + draws$u[, 2, ]
    bet_imb_i  <- draws$beta_bet[, 3]  + draws$u[, 3, ]
    side_run_i <- draws$beta_side[, 2] + draws$u[, 5, ]
    side_imb_i <- draws$beta_side[, 3] + draws$u[, 6, ]
    
    extra_tbl <- data.table(
      correlation = c("bet_run_slope <-> |side_run_slope| (derived, magnitude)",
                      "bet_imb_slope <-> |side_imb_slope| (derived, magnitude)",
                      "SELECTION: bet_intercept <-> side_run_slope (Omega[1,5])",
                      "SELECTION: bet_intercept <-> |side_imb_slope| (derived, magnitude)"),
      t(sapply(list(rowcor(bet_run_i, abs(side_run_i)),
                    rowcor(bet_imb_i, abs(side_imb_i)),
                    draws$Omega[, 1, 5],
                    rowcor(bet_int_i, abs(side_imb_i))), qsumm))
    )
    save_csv(extra_tbl, "synthesis_stan_dissociation_extra_")
    rm(bet_int_i, bet_run_i, bet_imb_i, side_run_i, side_imb_i)
    
    # ============================================================
    # HETEROGENEITY -- per-participant side-choice slopes, for the
    # debrief comparison (self-reported strategy vs. model-inferred type)
    # ============================================================
    # draws$u has dims [iterations, 6, P]; row 5 = side run-slope RE,
    # row 6 = side imbalance-slope RE (see synthesis.stan header for the
    # full u[] index map)
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
    
    
    save_csv(run_tbl, "synthesis_stan_run_slope_participants_")
    save_csv(imb_tbl, "synthesis_stan_imb_slope_participants_")
    
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
    
    
    save_csv(bet_run_tbl, "synthesis_stan_bet_run_slope_participants_")
    save_csv(bet_imb_tbl, "synthesis_stan_bet_imb_slope_participants_")
    
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
    save_csv(class_counts, "synthesis_stan_class_counts_")
    
    out[[tr]] <- list(fixed = fixed_tbl, dissociation = dissoc_tbl, dissociation_extra = extra_tbl,
                      within = within_tbl, tau = tau_tbl, ame = ame_tbl, cohens_h = cohen_tbl,
                      class_counts = class_counts,
                      run = run_tbl, imb = imb_tbl, bet_run = bet_run_tbl, bet_imb = bet_imb_tbl)
  }
  
  # ============================================================
  # FIGURE -- participant-level run-slope / imbalance-slope distributions,
  # FN (m25) vs FP (m19), overlaid in each panel for direct comparison
  # ============================================================
  # Two panels (run slope | imbalance slope), each showing a density curve
  # per treatment rather than a histogram -- overlaying bars from two
  # differently-sized samples (n=149 vs n=53) is hard to read side by
  # side, a density curve on a common x-axis is not. FN is drawn on top of
  # FP and listed first in the legend; the grey dotted line marks 0, the
  # colored dashed lines mark each treatment's population fixed effect.
  # Hardcoded to the two treatments this pipeline actually runs (unlike an
  # earlier, more generic version of this block, which supported an
  # arbitrary set of treatment codes via a fallback palette -- never
  # exercised in practice, dropped here for simplicity).
  #
  # Z-order is forced by literal ggplot layer-addition order (mk_panel adds
  # one geom_density + one geom_vline PER treatment, in draw_order, back to
  # front) rather than by factor-level sorting -- verified empirically that
  # relying on level order alone does not reliably put the last level on
  # top once a `breaks`-based legend reorder is also in play. Legend order
  # is set independently of draw order via explicit `breaks`.
  #
  # ggplot2 + patchwork -- both assumed loaded centrally (same convention
  # as data.table/rstan elsewhere in this pipeline: no library() calls in
  # this file). Add library(ggplot2) / library(patchwork) to 00_setup.R
  # if they aren't already there.
  if (all(c("m19", "m25") %in% names(out))) {
    
    tr_label <- c(m25 = "FN Treatment (m = 2.5)", m19 = "FP Treatment (m = 1.9)")
    tr_color <- c(m25 = "#2166AC", m19 = "#C0392B")
    draw_order   <- c("m19", "m25")   # back -> front: FP behind, FN on top
    legend_order <- c("m25", "m19")   # FN listed first in the legend
    
    level_order   <- tr_label[draw_order]     # factor levels follow draw order
    legend_breaks <- tr_label[legend_order]   # legend order set independently, via `breaks`
    pal <- setNames(tr_color[draw_order], level_order)
    
    mk_long <- function(field) rbindlist(lapply(draw_order, function(tr)
      data.table(label = factor(tr_label[[tr]], levels = level_order), slope = out[[tr]][[field]]$median)
    ))
    mk_pop <- function(pname) data.table(
      label = factor(tr_label[draw_order], levels = level_order),
      pop   = sapply(draw_order, function(tr) out[[tr]]$fixed[param == pname, median])
    )
    
    run_long <- mk_long("run");  pop_run_df <- mk_pop("side_run_slope")
    imb_long <- mk_long("imb");  pop_imb_df <- mk_pop("side_imb_slope")
    
    mk_panel <- function(long, pop_df, xlab) {
      p <- ggplot() + geom_vline(xintercept = 0, color = "grey40", linetype = "dotted")
      for (tr in draw_order) {
        lbl <- tr_label[[tr]]
        p <- p +
          geom_density(data = long[label == lbl], aes(x = slope, color = label, fill = label),
                       alpha = 0.2, linewidth = 0) +
          geom_vline(data = pop_df[label == lbl], aes(xintercept = pop, color = label),
                     linetype = "dashed", linewidth = 0.6, show.legend = TRUE)
      }
      p +
        scale_color_manual(values = pal, breaks = legend_breaks) +
        scale_fill_manual(values = pal, breaks = legend_breaks) +
        labs(x = xlab, y = "Density", color = NULL, fill = NULL) +
        theme_classic(base_size = 10) +
        theme(panel.grid.minor = element_blank(), legend.position = "bottom")
    }
    
    p_run <- mk_panel(run_long, pop_run_df, "Participant run-slope (posterior median)")
    p_imb <- mk_panel(imb_long, pop_imb_df, "Participant imbalance-slope (posterior median)")
    
    f_fig <- file.path(path_fig, "synthesis_stan_slope_distributions.png")
    if (!should_skip(f_fig, cfg, "output", "Synthesis slope-distribution figure (m25 vs m19)")) {
      fig <- p_run + p_imb + plot_layout(guides = "collect") & theme(legend.position = "bottom")
      ggsave(f_fig, fig, width = 9, height = 3, dpi = 300)
      msg("Saved: ", f_fig)
    }
  }
  
  invisible(out)
}