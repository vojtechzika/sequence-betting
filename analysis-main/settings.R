# ============================================================
# settings.R
#
# Settings a researcher may change between runs. Read once by
# run-analysis.R and bundled into cfg (cfg$run, cfg$model);
# the scripts read their settings from cfg.
#
# NOT HERE
#   scripts/config/01_design.R   experimental design and inferential thresholds
#                                (the study as preregistered)
#   selected-models.csv          model folder used for each model key; kept
#                                separate because it is edited between the
#                                steps of a run and is read at the moment a
#                                script needs it
#
# NOTES
#   - Convergence criteria are those of the preregistration. Sampler
#     settings (iter, warmup, adapt_delta, treedepth) may be raised until the
#     criteria are met; priors and likelihoods are not changed for that purpose.
#   - Settings of each Stan model are recorded in its model folder and in
#     run-environment.txt.
# ============================================================

settings <- list(
  
  # ----------------------------------------------------------
  # Run controls
  # ----------------------------------------------------------
  run = list(
    
    # Treatments to analyze (first = confirmatory)
    treatment = c("m25", "m19"),
    
    # Restrict RQ2/RQ3 to MPL-consistent participants only
    consistent_only = TRUE,
    
    # Global seed for reproducibility
    seed = 12345,
    
    # Prior-sensitivity analysis: models refitted with the
    # broadened priors of their prior_broad list, as an additional pass written
    # to <model folder>-broad (never selected in selected-models.csv).
    # Names are the blocks of model$stan below; c() switches the analysis off.
    prior_sensitivity = c("mpl_hierarchical", "rq1_bernoulli", "rq2_gaussian", "rq2_one_inflated",
                          "rq2_one_inflated_seq", "rq3_hurdle_gamma", "rq3_gamma",
                          "rq3_gaussian", "rq4_bernoulli")
  ),
  
  # ----------------------------------------------------------
  # Estimation, PPC, and simulation controls
  # ----------------------------------------------------------
  model = list(
    
    # --------------------------------------------------------
    # Stan sampling settings
    # --------------------------------------------------------
    stan = list(
      
      mpl_bernoulli = list(
        # partial_pooling = FALSE: individual-level model. The priors below are fixed, no population
        # parameters are estimated, so each participant's posterior depends only on that
        # participant's own choices. Estimates for any subset of participants (e.g. the
        # Holt-Laury-consistent ones) are therefore obtained by filtering the full-sample fit,
        # and the a* step relies on this. A model that pools across participants requires
        # a separate script and a separate consistent-sample path in the a* step.
        partial_pooling = FALSE,
        # Bland (2023), lognormal mean on natural scale, lognormal SD on log scale
        prior = list(r_mean = 0.27, r_sd = 0.36, lambda_mean = 30, lambda_sd = 0.5), 
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      mpl_hierarchical = list(
        # Population model: r_i ~ Normal(mu_r, sigma_r), log(lambda_i) ~ Normal(mu_ll, sigma_ll).
        # Hyperpriors: mu_r ~ Normal(mu_r_mean, mu_r_sd), sigma_r ~ HalfNormal(0, sigma_r_sd),
        # mu_ll ~ Normal(log(lambda_mean), mu_ll_sd), sigma_ll ~ HalfNormal(0, sigma_ll_sd).
        # Centred on the individual-level prior of mpl_bernoulli.
        prior = list(mu_r_mean = 0.27, mu_r_sd = 0.36, sigma_r_sd = 1,
                     lambda_mean = 30, mu_ll_sd = 0.5, sigma_ll_sd = 1),
        # broadened priors of the prior-sensitivity pass: the weakly informative scales of the
        # population model; the calibrated Bland (2023) centres and spreads are kept
        prior_broad = list(mu_r_mean = 0.27, mu_r_sd = 0.36, sigma_r_sd = 2,
                           lambda_mean = 30, mu_ll_sd = 0.5, sigma_ll_sd = 2),
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 10
      ),
      
      rq1_bernoulli = list(
        # intercept: Normal(0, alpha_sd), participant SD: HalfNormal(0, sigma_u_sd), sequence SD: HalfNormal(0, sigma_s_sd)
        prior = list(alpha_sd = 1.5, sigma_u_sd = 1, sigma_s_sd = 1),
        # broadened priors of the prior-sensitivity pass (preregistered: Normal(0, 2), HalfNormal(0, 2))
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2),
        iter = 8000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq2_gaussian = list(
        # intercept: Normal(0, alpha_sd); participant, sequence and residual SDs: HalfNormal(0, *_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2, sigma_sd = 2),
        iter = 8000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq2_one_inflated = list(
        # intercept: Normal(0, alpha_sd); participant and sequence SDs: HalfNormal(0, *_sd);
        # all-in propensity: intercept Normal(0, alpha_omega_sd), participant SD HalfNormal(0, sigma_v_sd);
        # Beta-Binomial concentration: HalfNormal(0, phi_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_v_sd = 1,
                     alpha_omega_sd = 1.5, phi_sd = 5),
        # broadened priors of the prior-sensitivity pass; phi_sd is already wider than 2 and is kept
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2, sigma_v_sd = 2,
                           alpha_omega_sd = 2, phi_sd = 5),
        iter = 12000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq2_one_inflated_seq = list(
        # as rq2_one_inflated, plus the SD of the sequence effect on the all-in
        # probability: HalfNormal(0, sigma_c_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_v_sd = 1,
                     sigma_c_sd = 1, alpha_omega_sd = 1.5, phi_sd = 5),
        # broadened priors of the prior-sensitivity pass; phi_sd is already wider than 2 and is kept
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2, sigma_v_sd = 2,
                           sigma_c_sd = 2, alpha_omega_sd = 2, phi_sd = 5),
        iter = 12000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq3_hurdle_gamma = list(
        # hurdle and positive-part intercepts (a0, ap): Normal(0, *_sd); participant and
        # sequence SDs of both parts (su0, sb0, sup, sbp): HalfNormal(0, *_sd)
        prior = list(a0_sd = 1, ap_sd = 1, su0_sd = 1, sb0_sd = 1, sup_sd = 1, sbp_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(a0_sd = 2, ap_sd = 2, su0_sd = 2, sb0_sd = 2, sup_sd = 2, sbp_sd = 2),
        iter = 12000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq3_gamma = list(
        # intercept: Normal(0, ap_sd); participant and sequence SDs: HalfNormal(0, *_sd)
        prior = list(ap_sd = 1, sup_sd = 1, sbp_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(ap_sd = 2, sup_sd = 2, sbp_sd = 2),
        iter = 8000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq3_gaussian = list(
        # intercept: Normal(0, alpha_sd); participant, sequence and residual SDs: HalfNormal(0, *_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2, sigma_sd = 2),
        iter = 20000, warmup = 2000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      
      rq4_bernoulli = list(
        # intercept: Normal(0, alpha_sd); participant and sequence SDs: HalfNormal(0, *_sd)
        prior = list(alpha_sd = 1.5, sigma_u_sd = 1, sigma_s_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(alpha_sd = 2, sigma_u_sd = 2, sigma_s_sd = 2),
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      ex1_regression = list(
        # intercept: Normal(0, alpha_sd); slopes: Normal(0, b_sd); residual SD: HalfNormal(0, sigma_sd)
        prior = list(alpha_sd = 1, b_sd = 1, sigma_sd = 1),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(alpha_sd = 2, b_sd = 2, sigma_sd = 2),
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.99, treedepth = 13
      ),
      
      ex2 = list(
        # outcome intercepts: Normal(0, alpha_sd); residual SDs: HalfNormal(0, sigma_sd);
        # pooled slopes: Normal(0, beta_bar_sd); SDs of the slopes across outcomes: HalfNormal(0, tau_sd)
        prior = list(alpha_sd = 1, sigma_sd = 1, beta_bar_sd = 1, tau_sd = 0.5),
        # broadened priors of the prior-sensitivity pass
        prior_broad = list(alpha_sd = 2, sigma_sd = 2, beta_bar_sd = 2, tau_sd = 2),
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.995, treedepth = 14
      ),
      
      synthesis = list(
        # fixed effects: Normal(0, beta_sd); participant and sequence-effect scales:
        # Student-t(3, 0, *_scale); participant-effect correlation: LKJ(lkj_eta)
        prior = list(beta_sd = 1.5, tau_scale = 1, lkj_eta = 2, sigma_v_scale = 1),
        # broadened priors of the prior-sensitivity pass; lkj_eta is a shape parameter and is kept
        prior_broad = list(beta_sd = 2, tau_scale = 2, lkj_eta = 2, sigma_v_scale = 2),
        iter = 12000, warmup = 2000, chains = 4, adapt_delta = 0.97, treedepth = 13
      )
    ),
    
    # --------------------------------------------------------
    # Convergence criteria (preregistered), applied to every fit:
    # Rhat_max <= rhat_max, ESS_min >= ess_min, divergences <= divergences
    # --------------------------------------------------------
    convergence = list(
      rhat_max    = 1.01,
      ess_min     = 400,
      divergences = 0L
    ),
    
    # --------------------------------------------------------
    # Posterior predictive replicates
    # --------------------------------------------------------
    ppc = list(
      rq1_interval = c(0.05, 0.95),
      rq1_k        = 1000L,
      rq2_k        = 1000L,
      rq3_k        = 1000L,
      rq4_k        = 1000L,
      rq4_mc_M     = 2000L,
      synthesis_k  = 1000L
    ),
    
    # --------------------------------------------------------
    # Simulation / resampling reps
    # --------------------------------------------------------
    simulation = list(
      ex1_trep = 200L,
      ex2_trep = 200L
    )
  ),

  # ----------------------------------------------------------
  # Discussion
  # ----------------------------------------------------------
  discussion = list(
    alpha = 0.05                    # level of the binomial classification test
  ),

  # ----------------------------------------------------------
  # Synthesis: descriptive checks outside the Stan model
  # ----------------------------------------------------------
  synthesis = list(
    debrief = list(
      alternative  = "two.sided",   # alternative of the directional Fisher tests
      min_coverage = 0.30,          # warn if share of codable strategies is lower
      n_boot       = 5000L          # bootstrap draws for the Likert contrasts
    ),
    oos = list(
      n_boot = 10000L               # bootstrap draws for the out-of-sample comparison
    ),
    costdir = list(
      n_draw   = 1000L,             # posterior draws used (subsampled from each fit)
      min_bets = 30L                # minimum number of bets in the robustness subset
    )
  )
)
