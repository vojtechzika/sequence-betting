// ------------------------------------------------------------
// RQ4: Side choices conditional on betting (Bernoulli-logit)
// h_t = 1 if Heads, 0 if Tails, only on trials with stake > 0
//
// h_t ~ Bernoulli(pi_is)
// logit(pi_is) = alpha + u_i + beta_s
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   beta_s ~ Normal(0, sigma_s), sum-to-zero via centering
//
// Priors:
// alpha ~ Normal(0, 1.5)
// sigma_u, sigma_s ~ HalfNormal(0, 1)
//
// Generated quantities:
// mu_h[s]   = mean_i inv_logit(alpha + u_i + beta_s)
// hbar      = mean_i inv_logit(alpha + u_i)          (baseline, beta_s=0)
// mu_h_i[i] = inv_logit(alpha + u_i)                 (participant baseline)
// ------------------------------------------------------------

data {
  int<lower=1> N;                 // participants
  int<lower=1> S;                 // sequences
  int<lower=1> T;                 // trials (betting trials only)
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> h[T];      // 1=Heads, 0=Tails
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

  // sum-to-zero sequence effects with hierarchical scale
  vector[S] beta = sigma_s * (b_raw - mean(b_raw));
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
    h[t] ~ bernoulli_logit(alpha + u[pid[t]] + beta[sid[t]]);
  }
}

generated quantities {
  vector[S] mu_h;
  real hbar;
  vector[N] mu_h_i;

  // participant baseline (beta_s = 0)
  for (i in 1:N) {
    mu_h_i[i] = inv_logit(alpha + u[i]);
  }

  // baseline averaged across participants
  {
    real acc = 0;
    for (i in 1:N) acc += mu_h_i[i];
    hbar = acc / N;
  }

  // per-sequence population mean
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + beta[s]);
    mu_h[s] = acc / N;
  }
}
