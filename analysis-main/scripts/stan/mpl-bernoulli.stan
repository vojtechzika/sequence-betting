// ------------------------------------------------------------
// Holt-Laury MPL: individual CRRA risk parameter r_i and choice
// precision lambda_i (individual-level, fixed priors)
//
// Likelihood: y_t ~ Bernoulli(logistic(lambda_i * (UA_t - UB_t)))
//   y_t = 1 if the safe lottery A is chosen in row t
//   U_k = p_t * u(k1; r_i) + (1 - p_t) * u(k2; r_i),  k in {A, B}
//   u(x; r) = x^(1 - r) / (1 - r), or log(x) if r = 1
//
// The four payoffs are the same in every row; only p_t varies. The
// utility difference is therefore computed from four utilities per
// participant:
//   UA_t - UB_t = p_t * d1_i + (1 - p_t) * d2_i,
//   d1_i = u(A1; r_i) - u(B1; r_i),  d2_i = u(A2; r_i) - u(B2; r_i)
//
// Priors:
//   r      ~ Normal(prior_r_mean, prior_r_sd)
//   lambda ~ Lognormal(log(prior_lambda_mean), prior_lambda_sd)
// ------------------------------------------------------------
functions {
  real crra_u(real x, real r) {
    if (fabs(r - 1.0) < 1e-8) return log(x);
    else return pow(x, 1.0 - r) / (1.0 - r);
  }
}
data {
  int<lower=1> N;                      // participants
  int<lower=1> T;                      // observations (participant x row)
  int<lower=1,upper=N> uid[T];         // participant of observation t
  vector<lower=0,upper=1>[T] p;        // probability of the high payoff in row t
  real<lower=0> A1;                    // lottery A, high payoff
  real<lower=0> A2;                    // lottery A, low payoff
  real<lower=0> B1;                    // lottery B, high payoff
  real<lower=0> B2;                    // lottery B, low payoff
  int<lower=0,upper=1> y[T];           // 1 = safe lottery A chosen
  real prior_r_mean;
  real<lower=0> prior_r_sd;
  real<lower=0> prior_lambda_mean;     // on natural scale; log taken in model
  real<lower=0> prior_lambda_sd;       // on log scale
}
parameters {
  vector[N] r;
  vector<lower=0>[N] lambda;
}
model {
  vector[N] d1;
  vector[N] d2;

  r ~ normal(prior_r_mean, prior_r_sd);
  lambda ~ lognormal(log(prior_lambda_mean), prior_lambda_sd);

  for (i in 1:N) {
    d1[i] = crra_u(A1, r[i]) - crra_u(B1, r[i]);
    d2[i] = crra_u(A2, r[i]) - crra_u(B2, r[i]);
  }

  y ~ bernoulli_logit(lambda[uid] .* (p .* d1[uid] + (1 - p) .* d2[uid]));
}
