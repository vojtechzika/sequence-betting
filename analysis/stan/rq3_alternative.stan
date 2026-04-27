// ------------------------------------------------------------
// RQ3: Certainty-equivalent welfare loss (Gaussian diagnostic)
// y_t = Δc_is / e
//
// Diagnostic alternative to the primary hurdle-Gamma model.
// Does not accommodate boundary at zero or right skew but
// provides a useful check on sign and ordering of sequence
// effects. Conclusions rely on stability across specifications.
//
// y_t ~ Normal(mu_t, sigma)
// mu_t = alpha + u_i + b_s [+ gamma_drift * block_c_t]
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Priors:
//   alpha ~ Normal(0, 1)
//   sigma, sigma_u, sigma_s ~ HalfNormal(0, 1)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//
// Generated quantities (evaluated at block_c = 0):
//   mu_c[s]   = mean_i(alpha + u_i + b_s)
//   mu_c_i[i] = mean_s(alpha + u_i + b_s)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector[T] y;
  int<lower=0,upper=1> include_drift;
  vector[T] block_c;
  real<lower=0> prior_gamma_sd;
}
parameters {
  real alpha;
  vector[N] u_raw;
  real<lower=0> sigma_u;
  vector[S] b_raw;
  real<lower=0> sigma_s;
  real<lower=0> sigma;
  real gamma_drift;
}
transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}
model {
  alpha       ~ normal(0, 1);
  sigma_u     ~ normal(0, 1);
  sigma_s     ~ normal(0, 1);
  sigma       ~ normal(0, 1);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  gamma_drift ~ normal(0, prior_gamma_sd);

  for (t in 1:T) {
    real mu_t = alpha + u[pid[t]] + b[sid[t]];
    if (include_drift) mu_t += gamma_drift * block_c[t];
    y[t] ~ normal(mu_t, sigma);
  }
}
generated quantities {
  vector[S] mu_c;
  vector[N] mu_c_i;

  // Evaluated at block_c = 0 (drift marginalized out)
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += alpha + u[i] + b[s];
    mu_c[s] = acc / N;
  }
  for (i in 1:N) {
    real acc = 0;
    for (s in 1:S) acc += alpha + u[i] + b[s];
    mu_c_i[i] = acc / S;
  }
}