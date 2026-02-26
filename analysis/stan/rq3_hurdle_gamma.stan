// ------------------------------------------------------------
// RQ3: Certainty-equivalent welfare loss (hurdle-gamma)
// y_t = Δc_is / e, can be zero
//
// hurdle:
//   P(y=0) = pi0_is
//   logit(pi0_is) = a0 + u0_i + b0_s
//
// positive part:
//   y | y>0 ~ Gamma(shape, rate = shape / mu_is)
//   log(mu_is) = ap + up_i + bp_s
//
// Random effects (each with sum-to-zero sequence effects):
//   u0_i ~ Normal(0, su0),  b0_s ~ Normal(0, sb0) centered
//   up_i ~ Normal(0, sup),  bp_s ~ Normal(0, sbp) centered
//
// Priors:
//   a0, ap ~ Normal(0, 1.5)
//   su0, sb0, sup, sbp ~ HalfNormal(0, 1) via <lower=0> normal(0,1)
//   shape ~ Gamma(2, 0.1)
//
// Generated quantities:
//   mu_c[s]   = mean_i E[y_{is} | params]   (per-sequence population mean)
//   mu_c_i[i] = mean_s E[y_{is} | params]   (per-participant mean across sequences)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector<lower=0>[T] y;
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

  real<lower=0> shape;
}

transformed parameters {
  vector[N] u0 = su0 * u0_raw;
  vector[N] up = sup * up_raw;

  // sum-to-zero sequence effects with hierarchical scales
  vector[S] b0 = sb0 * (b0_raw - mean(b0_raw));
  vector[S] bp = sbp * (bp_raw - mean(bp_raw));
}

model {
  // priors
  a0  ~ normal(0, 1.5);
  ap  ~ normal(0, 1.5);

  su0 ~ normal(0, 1);
  sb0 ~ normal(0, 1);
  sup ~ normal(0, 1);
  sbp ~ normal(0, 1);

  u0_raw ~ normal(0, 1);
  up_raw ~ normal(0, 1);
  b0_raw ~ normal(0, 1);
  bp_raw ~ normal(0, 1);

  shape ~ gamma(2, 0.1);

  // likelihood
  for (t in 1:T) {
    real pi0 = inv_logit(a0 + u0[pid[t]] + b0[sid[t]]);
    if (y[t] == 0) {
      target += bernoulli_lpmf(1 | pi0);
    } else {
      real mu_pos = exp(ap + up[pid[t]] + bp[sid[t]]);
      target += bernoulli_lpmf(0 | pi0);
      target += gamma_lpdf(y[t] | shape, shape / mu_pos);
    }
  }
}

generated quantities {
  vector[S] mu_c;
  vector[N] mu_c_i;

  // per-sequence population mean: mean_i E[y_is]
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      real pi0 = inv_logit(a0 + u0[i] + b0[s]);
      real mu_pos = exp(ap + up[i] + bp[s]);
      acc += (1 - pi0) * mu_pos;
    }
    mu_c[s] = acc / N;
  }

  // per-participant mean across sequences: mean_s E[y_is]
  for (i in 1:N) {
    real acc = 0;
    for (s in 1:S) {
      real pi0 = inv_logit(a0 + u0[i] + b0[s]);
      real mu_pos = exp(ap + up[i] + bp[s]);
      acc += (1 - pi0) * mu_pos;
    }
    mu_c_i[i] = acc / S;
  }
}
