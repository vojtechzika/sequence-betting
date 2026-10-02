// EX1.2: Predictors of GHI with uncertainty propagation via MI likelihood
// y_rep[t,i] are draws of chi_i^(t). We integrate over t with equal weights.
//
// Model:
//   chi_i^(t) ~ Normal(mu_i, sigma)
//   mu_i = alpha + b_opt * Zopt_i + b_rt * Zrt_i + b_r * Zr_i
//
// Priors:
//   alpha, b_* ~ Normal(0,1)
//   sigma ~ HalfNormal(0,1)

data {
  int<lower=1> N;                 // participants
  int<lower=1> Trep;              // number of retained chi-draws used for MI
  matrix[Trep, N] y_rep;          // chi draws (Trep x N)

  vector[N] Zopt;                 // standardized LOT-R
  vector[N] Zrt;                  // standardized log RT
  vector[N] Zr;                   // standardized risk (r)
}

parameters {
  real alpha;
  real b_opt;
  real b_rt;
  real b_r;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu = alpha + b_opt * Zopt + b_rt * Zrt + b_r * Zr;
}

model {
  // priors
  alpha ~ normal(0, 1);
  b_opt ~ normal(0, 1);
  b_rt  ~ normal(0, 1);
  b_r   ~ normal(0, 1);
  sigma ~ normal(0, 1);

  // MI likelihood: for each participant, average over Trep imputations (draws)
  for (i in 1:N) {
    vector[Trep] lp;
    for (t in 1:Trep) {
      lp[t] = normal_lpdf(y_rep[t, i] | mu[i], sigma);
    }
    target += log_sum_exp(lp) - log(Trep);
  }
}

generated quantities {
  // Useful for reporting
  real p_b_opt_pos = b_opt > 0;
  real p_b_rt_pos  = b_rt  > 0;
  real p_b_r_pos   = b_r   > 0;
}
