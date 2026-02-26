// ------------------------------------------------------------
// RQ1: Probability of betting (Bernoulli-logit)
// y_t = 1(stake > 0)
//
// y_t ~ Bernoulli(pi_is)
// logit(pi_is) = alpha + u_i + b_s
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Priors:
//   alpha ~ Normal(0, 1.5)
//   sigma_u, sigma_s ~ HalfNormal(0, 1) via <lower=0> normal(0,1)
//
// Generated quantities:
//   mu_b[s]   = mean_i inv_logit(alpha + u_i + b_s)
//   mu_b_i[i] = mean_s inv_logit(alpha + u_i + b_s)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> y[T];
}

parameters {
  real alpha;

  vector[N] u_raw;
  real<lower=0> sigma_u;

  vector[S] b_raw;
  real<lower=0> sigma_s;
}

transformed parameters {
  vector[N] u = sigma_u * u_raw;

  // sum-to-zero sequence effects (identified relative to grand mean)
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}

model {
  // priors
  alpha   ~ normal(0, 1.5);
  sigma_u ~ normal(0, 1);
  sigma_s ~ normal(0, 1);

  u_raw ~ normal(0, 1);
  b_raw ~ normal(0, 1);

  // likelihood
  for (t in 1:T) {
    y[t] ~ bernoulli_logit(alpha + u[pid[t]] + b[sid[t]]);
  }
}

generated quantities {
  vector[S] mu_b;
  vector[N] mu_b_i;

  // per-sequence population mean betting prob
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b[s] = acc / N;
  }

  // per-participant mean betting prob (averaged across sequences)
  for (i in 1:N) {
    real acc = 0;
    for (s in 1:S) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b_i[i] = acc / S;
  }
}
