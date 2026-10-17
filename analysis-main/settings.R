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
    seed = 12345
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
      
      rq1_bernoulli = list(
        # intercept: Normal(0, alpha_sd), participant SD: HalfNormal(0, sigma_u_sd), sequence SD: HalfNormal(0, sigma_s_sd)
        prior = list(alpha_sd = 1.5, sigma_u_sd = 1, sigma_s_sd = 1),
        iter = 8000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq2_gaussian = list(
        # intercept: Normal(0, alpha_sd); participant, sequence and residual SDs: HalfNormal(0, *_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_sd = 1),
        iter = 8000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq2_one_inflated = list(
        # intercept: Normal(0, alpha_sd); participant and sequence SDs: HalfNormal(0, *_sd);
        # all-in propensity: intercept Normal(0, alpha_omega_sd), participant SD HalfNormal(0, sigma_v_sd);
        # Beta-Binomial concentration: HalfNormal(0, phi_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_v_sd = 1,
                     alpha_omega_sd = 1.5, phi_sd = 5),
        iter = 12000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq2_one_inflated_seq = list(
        # as rq2_one_inflated, plus the SD of the sequence effect on the all-in
        # probability: HalfNormal(0, sigma_c_sd)
        prior = list(alpha_sd = 1, sigma_u_sd = 1, sigma_s_sd = 1, sigma_v_sd = 1,
                     sigma_c_sd = 1, alpha_omega_sd = 1.5, phi_sd = 5),
        iter = 8000, warmup = 2000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq3_hurdle_gamma = list(
        iter = 12000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq3_gamma = list(
        iter = 8000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      rq3_gaussian = list(
        iter = 8000, warmup = 1000, chains = 4, adapt_delta = 0.97, treedepth = 13
      ),
      
      
      rq4_bernoulli = list(
        # intercept: Normal(0, alpha_sd); participant and sequence SDs: HalfNormal(0, *_sd)
        prior = list(alpha_sd = 1.5, sigma_u_sd = 1, sigma_s_sd = 1),
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      ex1_regression = list(
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.99, treedepth = 13
      ),
      
      ex2 = list(
        iter = 3000, warmup = 1000, chains = 4, adapt_delta = 0.995, treedepth = 14
      ),
      
      synthesis = list(
        # fixed effects: Normal(0, beta_sd); participant and sequence-effect scales:
        # Student-t(3, 0, *_scale); participant-effect correlation: LKJ(lkj_eta)
        prior = list(beta_sd = 1.5, tau_scale = 1, lkj_eta = 2, sigma_v_scale = 1),
        iter = 12000, warmup = 1500, chains = 4, adapt_delta = 0.97, treedepth = 13
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
