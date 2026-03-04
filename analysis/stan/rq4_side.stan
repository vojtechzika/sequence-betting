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
  vector[S] mu_h;           // population mean per sequence: E_u logistic(alpha + u + beta_s)
  real hbar;                // population baseline:          E_u logistic(alpha + u)
  vector[N] mu_h_i;         // participant baselines: logistic(alpha + u_i)

  // OPTIONAL sanity/traceability:
  vector[S] mu_h_sample;    // sample mean over realized u_i: mean_i logistic(alpha + u_i + beta_s)
  real hbar_sample;         // sample baseline:               mean_i logistic(alpha + u_i)

  // participant-specific baselines
  for (i in 1:N)
    mu_h_i[i] = inv_logit(alpha + u[i]);

  // sample-based summaries (exact given u)
  {
    real acc0 = 0;
    for (i in 1:N) acc0 += inv_logit(alpha + u[i]);
    hbar_sample = acc0 / N;

    for (s in 1:S) {
      real accs = 0;
      for (i in 1:N) accs += inv_logit(alpha + u[i] + beta[s]);
      mu_h_sample[s] = accs / N;
    }
  }

  // population-based summaries via MC with common random numbers
  {
    int M = 2000;
    array[M] real u_new;
    real acc0 = 0;

    for (m in 1:M) {
      u_new[m] = normal_rng(0, sigma_u);
      acc0 += inv_logit(alpha + u_new[m]);
    }
    hbar = acc0 / M;

    for (s in 1:S) {
      real accs = 0;
      for (m in 1:M) accs += inv_logit(alpha + u_new[m] + beta[s]);
      mu_h[s] = accs / M;
    }
  }
}

