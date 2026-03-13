// ============================================================
// EX2: Participant-level associations with optimism and response time
// across RQ1--RQ4 outcomes, with partial pooling across outcomes.
// Uncertainty propagation via MI likelihood over posterior draws.
//
// Data are stacked across outcomes k=1..K with per-row outcome id kid[n].
// For each n, y_rep[t,n] provides the t-th posterior draw of standardized outcome.
//
// Model:
//   y_{t,n} ~ Normal(mu_n, sigma_{kid[n]})
//   mu_n = alpha_{kid[n]} 
//          + beta_opt_{kid[n]}*Zopt_n 
//          + beta_rt_{kid[n]}*Zrt_n
//
// Partial pooling on slopes:
//   beta_*_k = beta_*_bar + tau_* * beta_*_raw[k]
//
// Priors:
//   beta_*_bar ~ Normal(0,1)
//   tau_* ~ HalfNormal(0,0.5)
//   alpha_k ~ Normal(0,1)
//   sigma_k ~ HalfNormal(0,1)
// ============================================================

data {
  int<lower=1> K;                      // number of outcomes (4)
  int<lower=1> Nobs;                   // stacked rows
  int<lower=1> Trep;                   // MI posterior draws

  matrix[Trep, Nobs] y_rep;            // outcome draws

  int<lower=1,upper=K> kid[Nobs];      // outcome id

  vector[Nobs] Zopt;                   // optimism
  vector[Nobs] Zrt;                    // log response time
}

parameters {
  vector[K] alpha_k;
  vector<lower=0>[K] sigma_k;

  real beta_opt_bar;
  real beta_rt_bar;

  real<lower=0> tau_opt;
  real<lower=0> tau_rt;

  vector[K] beta_opt_raw;
  vector[K] beta_rt_raw;
}

transformed parameters {
  vector[K] beta_opt_k = beta_opt_bar + tau_opt * beta_opt_raw;
  vector[K] beta_rt_k  = beta_rt_bar  + tau_rt  * beta_rt_raw;
}

model {

  // priors
  alpha_k ~ normal(0, 1);
  sigma_k ~ normal(0, 1);

  beta_opt_bar ~ normal(0, 1);
  beta_rt_bar  ~ normal(0, 1);

  tau_opt ~ normal(0, 0.5);
  tau_rt  ~ normal(0, 0.5);

  beta_opt_raw ~ normal(0, 1);
  beta_rt_raw  ~ normal(0, 1);

  // MI likelihood
  for (n in 1:Nobs) {

    real mu = alpha_k[kid[n]]
              + beta_opt_k[kid[n]] * Zopt[n]
              + beta_rt_k[kid[n]]  * Zrt[n];

    vector[Trep] lp;

    for (t in 1:Trep)
      lp[t] = normal_lpdf(y_rep[t, n] | mu, sigma_k[kid[n]]);

    target += log_sum_exp(lp) - log(Trep);
  }
}

generated quantities {

  // pooled slope direction probabilities
  real p_beta_opt_bar_pos = beta_opt_bar > 0;
  real p_beta_rt_bar_pos  = beta_rt_bar  > 0;

}
