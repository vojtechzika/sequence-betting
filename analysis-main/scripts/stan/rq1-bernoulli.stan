// ------------------------------------------------------------
// RQ1: Probability of betting
// y_t = 1(stake > 0)
//
// Likelihood: Bernoulli-logit
//   y_t ~ Bernoulli(pi_ist)
//   logit(pi_ist) = alpha + u_i + b_s [+ gamma_drift * block_c_t]
//
// Convergence and posterior predictive checks are computed in
// 11_rq1-stan-bernoulli.R. A replacement likelihood, if ever required,
// would be a separate Stan file with its own fitting script; this file
// is not modified.
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s), sum-to-zero via centering
//
// Optional linear drift (include_drift = 1):
//   gamma_drift * block_c_t, block_c in {-1.5, -0.5, 0.5, 1.5}
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//   When include_drift = 0, gamma_drift is sampled but not used in
//   the likelihood; the prior shrinks it toward zero.
//
// Priors (scales are passed as data):
//   alpha       ~ Normal(0, prior_alpha_sd)
//   sigma_u     ~ HalfNormal(0, prior_sigma_u_sd)
//   sigma_s     ~ HalfNormal(0, prior_sigma_s_sd)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//
// Generated quantities (all evaluated at block_c = 0):
//   mu_b[s]        = mean_i inv_logit(alpha + u_i + b_s)
//   mu_b_i[i]      = mean_s inv_logit(alpha + u_i + b_s)
//   grand_mean_bet = mean_{i,s} inv_logit(alpha + u_i + b_s)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> y[T];
  int<lower=0,upper=1> include_drift;   // 1 = include linear drift, 0 = no drift
  vector[T] block_c;                    // centered block index {-1.5,-0.5,0.5,1.5}
  real<lower=0> prior_gamma_sd;         // prior SD for drift coefficient
  real<lower=0> prior_alpha_sd;         // prior SD for the intercept
  real<lower=0> prior_sigma_u_sd;       // prior scale for sigma_u (half-normal)
  real<lower=0> prior_sigma_s_sd;       // prior scale for sigma_s (half-normal)
}
parameters {
  real alpha;
  vector[N] u_raw;
  real<lower=0> sigma_u;
  vector[S] b_raw;
  real<lower=0> sigma_s;
  real gamma_drift;                     // always sampled; prior shrinks to 0 when not needed
}
transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}
model {
  alpha       ~ normal(0, prior_alpha_sd);
  sigma_u     ~ normal(0, prior_sigma_u_sd);
  sigma_s     ~ normal(0, prior_sigma_s_sd);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  gamma_drift ~ normal(0, prior_gamma_sd);

  for (t in 1:T) {
    real eta = alpha + u[pid[t]] + b[sid[t]];
    if (include_drift) eta += gamma_drift * block_c[t];
    y[t] ~ bernoulli_logit(eta);
  }
}
generated quantities {
  vector[S] mu_b;
  vector[N] mu_b_i;
  real grand_mean_bet;

  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b[s] = acc / N;
  }
  for (i in 1:N) {
    real acc = 0;
    for (s in 1:S) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b_i[i] = acc / S;
  }
  {
    real acc = 0;
    for (i in 1:N)
      for (s in 1:S)
        acc += inv_logit(alpha + u[i] + b[s]);
    grand_mean_bet = acc / (N * S);
  }
}
