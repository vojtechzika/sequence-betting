// ============================================================
// stan/ex2_associations.stan
//
// EX2: Participant-level associations with optimism and
// response time across RQ1--RQ4 outcomes, with partial
// pooling of outcome-specific slopes.
//
// Posterior uncertainty in participant-level outcomes from
// the upstream RQ models is propagated via a multiple-
// imputation (MI) likelihood: for each stacked observation n
// the likelihood is the average over Trep posterior draws of
// the standardized outcome.
//
// Model:
//   y_{t,n} ~ Normal(mu_n, sigma_{kid[n]})
//   mu_n = alpha_{kid[n]}
//          + beta_opt_{kid[n]} * Zopt_n
//          + beta_rt_{kid[n]}  * Zrt_n
//
// Partial pooling on slopes (non-centred parametrisation):
//   beta_opt_k = beta_opt_bar + tau_opt * beta_opt_raw[k]
//   beta_rt_k  = beta_rt_bar  + tau_rt  * beta_rt_raw[k]
//
// Priors (weakly informative, per preregistration):
//   alpha_k      ~ Normal(0, 1)
//   sigma_k      ~ HalfNormal(0, 1)    [implied by lower=0 + normal(0,1)]
//   beta_opt_bar ~ Normal(0, 1)
//   beta_rt_bar  ~ Normal(0, 1)
//   tau_opt      ~ HalfNormal(0, 0.5)
//   tau_rt       ~ HalfNormal(0, 0.5)
//   beta_opt_raw ~ Normal(0, 1)
//   beta_rt_raw  ~ Normal(0, 1)
//
// Data are stacked across outcomes: row n has outcome index
// kid[n] in {1,...,K}, predictors Zopt[n] and Zrt[n], and
// Trep standardized outcome draws y_rep[1..Trep, n].
// K = 4 is fixed (b=1, a=2, c=3, h=4) regardless of which
// outcomes have observations, so Stan parameter indexing is
// stable across treatments and filter configurations.
// ============================================================

data {
  int<lower=1>                    K;       // number of outcomes (4)
  int<lower=1>                    Nobs;    // total stacked participants
  int<lower=1>                    Trep;    // MI posterior draws per participant
  matrix[Trep, Nobs]              y_rep;   // standardized outcome draws
  array[Nobs] int<lower=1,upper=K> kid;   // outcome id per observation
  vector[Nobs]                    Zopt;    // z-scored LOT-R optimism
  vector[Nobs]                    Zrt;     // z-scored log response time
}

parameters {
  vector[K]          alpha_k;       // outcome-specific intercepts
  vector<lower=0>[K] sigma_k;       // outcome-specific residual SDs
  real               beta_opt_bar;  // pooled optimism slope
  real               beta_rt_bar;   // pooled response-time slope
  real<lower=0>      tau_opt;       // between-outcome SD of optimism slopes
  real<lower=0>      tau_rt;        // between-outcome SD of RT slopes
  vector[K]          beta_opt_raw;  // non-centred raw slopes (optimism)
  vector[K]          beta_rt_raw;   // non-centred raw slopes (RT)
}

transformed parameters {
  // Outcome-specific slopes (non-centred parametrisation)
  vector[K] beta_opt_k = beta_opt_bar + tau_opt * beta_opt_raw;
  vector[K] beta_rt_k  = beta_rt_bar  + tau_rt  * beta_rt_raw;
}

model {
  // ---- Priors ----
  alpha_k      ~ normal(0, 1);
  sigma_k      ~ normal(0, 1);    // HalfNormal(0,1) given lower=0 constraint
  beta_opt_bar ~ normal(0, 1);
  beta_rt_bar  ~ normal(0, 1);
  tau_opt      ~ normal(0, 0.5);  // HalfNormal(0, 0.5)
  tau_rt       ~ normal(0, 0.5);
  beta_opt_raw ~ normal(0, 1);
  beta_rt_raw  ~ normal(0, 1);

  // ---- MI likelihood ----
  // For each observation n, the likelihood marginalises over Trep
  // posterior draws of the standardized outcome:
  //   log p += log[ (1/Trep) * sum_t Normal(y_rep[t,n] | mu_n, sigma_k) ]
  //          = log_sum_exp_t Normal_lpdf - log(Trep)
  for (n in 1:Nobs) {
    real mu_n    = alpha_k[kid[n]]
                   + beta_opt_k[kid[n]] * Zopt[n]
                   + beta_rt_k[kid[n]]  * Zrt[n];
    real sigma_n = sigma_k[kid[n]];
    vector[Trep] lp;
    for (t in 1:Trep)
      lp[t] = normal_lpdf(y_rep[t, n] | mu_n, sigma_n);
    target += log_sum_exp(lp) - log(Trep);
  }
}

generated quantities {
  // Posterior probability that each pooled slope is positive
  int<lower=0,upper=1> p_beta_opt_bar_pos = (beta_opt_bar > 0);
  int<lower=0,upper=1> p_beta_rt_bar_pos  = (beta_rt_bar  > 0);

  // Outcome-specific direction probabilities
  array[K] int<lower=0,upper=1> p_beta_opt_k_pos;
  array[K] int<lower=0,upper=1> p_beta_rt_k_pos;
  for (k in 1:K) {
    p_beta_opt_k_pos[k] = (beta_opt_k[k] > 0);
    p_beta_rt_k_pos[k]  = (beta_rt_k[k]  > 0);
  }
}
