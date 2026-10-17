# ============================================================
# 94_synthesis-joint-effects.R
#
# PURPOSE
#   Tables of the joint model of betting and side choice
#   (91_synthesis-stan-joint.R), one set per treatment:
#     fixed effects            population intercepts and cue slopes of the
#                              betting equation (run, imbalance magnitude)
#                              and the side equation (signed run, signed
#                              imbalance), log-odds scale
#     dissociation             cross-equation correlations of the
#                              participant effects (intercepts, run slopes,
#                              imbalance slopes)
#     within-corr              the six within-equation correlations
#     dissociation-extra       correlations of the betting slopes with the
#                              magnitude of the personal side slopes, and
#                              the two selection checks (a betting
#                              intercept that correlates with the side
#                              slopes would indicate that selective bettors
#                              are more directional). Rows marked "derived"
#                              use personal effects (fixed plus random
#                              effect), which are shrunk and therefore
#                              conservative
#     tau                      participant-level SDs of the six random
#                              effects and their ratio to the absolute
#                              population effect
#     ame                      population-averaged marginal effects on the
#                              probability scale, dp/dx at each trial's
#                              actual covariates and random effects, and the
#                              random-effect-integrated baseline
#                              probabilities
#     cohens-h                 draw-level Cohen's h for the run and
#                              imbalance cues, with the two counterfactual
#                              probabilities behind each h
#   All summaries are posterior medians with 95% credible intervals and the
#   share of draws above zero.
#
#   Sign coding: betting slopes are unsigned (bet more when a run is present
#   or the sequence is imbalanced); side slopes are signed toward Heads (a
#   negative run slope means betting against the run, a positive imbalance
#   slope betting on the majority symbol). The side intercept is a baseline
#   preference for the Heads label.
#
# SAMPLE
#   synthesis-full of each treatment in cfg$run$treatment.
#
# INPUT
#   model "synthesis": fit-, pid-levels-, data-<tr>-full.rds
#
# OUTPUT (data/output/csv/)
#   synthesis-<tr>-joint-<name>.csv   name in fixed-effects, dissociation,
#                                     within-corr, dissociation-extra, tau,
#                                     ame, cohens-h
#
# NOTES
#   - AME and Cohen's h are computed for a subsample of 1,000 posterior draws
#     (seeded from cfg$run$seed), at each draw's actual participant and
#     sequence effects; benchmarks for h: 0.2 small, 0.5 medium, 0.8 large
#   - Output files are overwritten at every run
# ============================================================

synthesis_joint_effects <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  smp    <- "full"

  qsumm <- function(x) c(median = median(x), q025 = unname(quantile(x, .025)),
                         q975 = unname(quantile(x, .975)), p_gt0 = mean(x > 0))

  # Cohen's h for two proportions: difference of arcsine-transformed
  # probabilities, h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p0)). Signed: positive
  # when p1 > p0. Conventional benchmarks are |h| = 0.2 small, 0.5 medium,
  # 0.8 large.
  cohen_h <- function(p1, p0) 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p0))
  

  for (tr in tr_vec) {

    save_csv <- function(x, name) {
      f <- file.path(path_csv, paste0("synthesis-", tr, "-joint-", name, ".csv"))
      fwrite(x, f);  msg("Saved: ", f)
    }

    fit <- readRDS(model_path("synthesis", paste0("fit-", tr, "-", smp, ".rds")))
    dl  <- readRDS(model_path("synthesis", paste0("data-", tr, "-", smp, ".rds")))
    pid_levels <- readRDS(model_path("synthesis", paste0("pid-levels-", tr, "-", smp, ".rds")))

    draws <- rstan::extract(fit, pars = c(
      "beta_bet", "beta_side", "tau", "sigma_v_bet", "sigma_v_side",
      "corr_intercepts", "corr_run_slopes", "corr_imb_slopes", "Omega",
      "u", "v_bet", "v_side"
    ))

    tau_labels <- c("bet_intercept_sd", "bet_run_sd", "bet_imb_sd",
                    "side_intercept_sd", "side_run_sd", "side_imb_sd")
    tau_med <- apply(draws$tau, 2, median)
    tau_lo  <- apply(draws$tau, 2, quantile, probs = .025)
    tau_hi  <- apply(draws$tau, 2, quantile, probs = .975)
    tau_rhat <- summary(fit, pars = "tau")$summary[, "Rhat"]
    
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
    save_csv(tau_tbl, "tau")

    set.seed(as.integer(cfg$run$seed) + 3001L)
    k <- 1000L
    n_draws <- nrow(draws$beta_bet)
    idx <- sample(seq_len(n_draws), min(k, n_draws))
    
    
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
    # and bet_intercept_sd/side_intercept_sd (see the tau output file) show that
    # the heterogeneity is not negligible, especially for the betting
    # equation.
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
    save_csv(ame_tbl, "ame")
    
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
    save_csv(cohen_tbl, "cohens-h")
    
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
    
    save_csv(fixed_tbl, "fixed-effects")
    save_csv(dissoc_tbl, "dissociation")
    save_csv(within_tbl, "within-corr")
    
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
    save_csv(extra_tbl, "dissociation-extra")
    rm(bet_int_i, bet_run_i, bet_imb_i, side_run_i, side_imb_i)
  }

  invisible(TRUE)
}
