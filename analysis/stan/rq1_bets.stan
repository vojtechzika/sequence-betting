// ------------------------------------------------------------
// RQ1: Probability of betting
// y_t = 1(stake > 0)
//
// Primary likelihood (use_bb = 0): Bernoulli-logit
//   y_t ~ Bernoulli(pi_ist)
//
// Robustness likelihood (use_bb = 1): Beta-Binomial
//   y_t ~ BetaBinomial(1, mu_t * phi, (1 - mu_t) * phi)
//   where mu_t = inv_logit(eta_t)
//   phi ~ HalfNormal(0, 5)  -- overdispersion parameter
//
// Linear predictor:
//   logit(pi_ist) = alpha + u_i + b_s [+ gamma_drift * block_c_t]
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Optional linear drift (include_drift = 1):
//   gamma_drift * block_c_t
//   block_c in {-1.5, -0.5, 0.5, 1.5} (centered block index)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//   When include_drift = 0, gamma_drift is sampled but not used
//   in the likelihood; the prior shrinks it toward zero.
//
// Priors:
//   alpha       ~ Normal(0, 1.5)
//   sigma_u     ~ HalfNormal(0, 1)
//   sigma_s     ~ HalfNormal(0, 1)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//   phi         ~ HalfNormal(0, 5)  [only used when use_bb = 1]
//
// Data inputs:
//   use_bb         : 0 = Bernoulli (primary), 1 = Beta-Binomial (robustness)
//   include_drift  : 0/1 flag for drift adjustment
//   block_c        : centered block index for each observation
//   prior_gamma_sd : prior SD for drift coefficient (from design_cfg)
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
  int<lower=0,upper=1> use_bb;          // 0 = Bernoulli, 1 = Beta-Binomial
  int<lower=0,upper=1> include_drift;   // 1 = include linear drift, 0 = no drift
  vector[T] block_c;                    // centered block index {-1.5,-0.5,0.5,1.5}
  real<lower=0> prior_gamma_sd;         // prior SD for drift coefficient
}
parameters {
  real alpha;
  vector[N] u_raw;
  real<lower=0> sigma_u;
  vector[S] b_raw;
  real<lower=0> sigma_s;
  real gamma_drift;                     // always sampled; prior shrinks to 0 when not needed
  real<lower=0> phi;                    // BB overdispersion; sampled always, used only if use_bb
}
transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}
model {
  // priors
  alpha       ~ normal(0, 1.5);
  sigma_u     ~ normal(0, 1);
  sigma_s     ~ normal(0, 1);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  gamma_drift ~ normal(0, prior_gamma_sd);
  phi         ~ normal(0, 5);           // HalfNormal via <lower=0> constraint

  // likelihood
  for (t in 1:T) {
    real eta = alpha + u[pid[t]] + b[sid[t]];
    if (include_drift) eta += gamma_drift * block_c[t];
    if (use_bb) {
      real mu = inv_logit(eta);
      y[t] ~ beta_binomial(1, mu * phi, (1 - mu) * phi);
    } else {
      y[t] ~ bernoulli_logit(eta);
    }
  }
}
generated quantities {
  vector[S] mu_b;
  vector[N] mu_b_i;
  real grand_mean_bet;

  // per-sequence population mean betting prob (at block_c = 0)
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b[s] = acc / N;
  }
  // per-participant mean betting prob (at block_c = 0)
  for (i in 1:N) {
    real acc = 0;
    for (s in 1:S) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b_i[i] = acc / S;
  }
  // grand mean (at block_c = 0)
  {
    real acc = 0;
    for (i in 1:N)
      for (s in 1:S)
        acc += inv_logit(alpha + u[i] + b[s]);
    grand_mean_bet = acc / (N * S);
  }
}