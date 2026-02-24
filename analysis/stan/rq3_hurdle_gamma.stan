// ------------------------------------------------------------
// RQ3: Certainty-equivalent welfare loss (hurdle-Gamma) -- PREREG PRIMARY
// Outcome y_t = Δc_is / e (non-negative; includes zeros)
//
// Hurdle part:
//   Pr(y_t = 0) = pi0_t
//   logit(pi0_t) = a0 + u0_i + b0_s
//
// Positive part:
//   y_t | (y_t > 0) ~ Gamma(shape, rate)
//   mean_t = exp(ap + up_i + bp_s)
//   rate_t = shape / mean_t
//
// Random effects (sum-to-zero for sequence effects via centering):
//   u0_i ~ Normal(0, su0) ; b0_s raw ~ Normal(0, sb0)
//   up_i ~ Normal(0, sup) ; bp_s raw ~ Normal(0, sbp)
//
// Priors (weakly informative):
//   a0, ap ~ Normal(0,1)
//   SDs ~ HalfNormal(0,1) via <lower=0>
//   shape ~ lognormal(0,1)
//
// Generated quantities:
//   mu_c[s] = E_i[ (1-pi0_is)*mean_pos_is ]  (expected loss share by sequence)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector<lower=0>[T] y;          // includes zeros
  int<lower=0,upper=1> is_zero[T];
}

parameters {
  // hurdle part
  real a0;
  vector[N] u0_raw;
  real<lower=0> su0;
  vector[S] b0_raw;
  real<lower=0> sb0;

  // positive mean part
  real ap;
  vector[N] up_raw;
  real<lower=0> sup;
  vector[S] bp_raw;
  real<lower=0> sbp;

  // gamma shape
  real<lower=0> shape;
}

transformed parameters {
  vector[N] u0 = su0 * u0_raw;
  vector[N] up = sup * up_raw;

  vector[S] b0 = b0_raw - mean(b0_raw);
  vector[S] bp = bp_raw - mean(bp_raw);
}

model {
  a0  ~ normal(0, 1);
  ap  ~ normal(0, 1);

  su0 ~ normal(0, 1);
  sb0 ~ normal(0, 1);
  sup ~ normal(0, 1);
  sbp ~ normal(0, 1);

  u0_raw ~ normal(0, 1);
  up_raw ~ normal(0, 1);

  b0_raw ~ normal(0, sb0);
  bp_raw ~ normal(0, sbp);

  shape ~ lognormal(0, 1);

  for (t in 1:T) {
    real logit_pi0 = a0 + u0[pid[t]] + b0[sid[t]];
    real pi0 = inv_logit(logit_pi0);

    if (is_zero[t] == 1) {
      target += bernoulli_lpmf(1 | pi0);
    } else {
      real mean_pos = exp(ap + up[pid[t]] + bp[sid[t]]);
      real rate = shape / mean_pos;

      target += bernoulli_lpmf(0 | pi0);
      target += gamma_lpdf(y[t] | shape, rate);
    }
  }
}

generated quantities {
  vector[S] mu_c;

  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      real pi0 = inv_logit(a0 + u0[i] + b0[s]);
      real mean_pos = exp(ap + up[i] + bp[s]);
      acc += (1 - pi0) * mean_pos;
    }
    mu_c[s] = acc / N;
  }
}
