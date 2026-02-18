functions {
  real crra_u(real x, real r) {
    if (fabs(r - 1.0) < 1e-6) return log(x);
    return pow(x, 1.0 - r) / (1.0 - r);
  }
}

data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1,upper=N> uid[T];
  vector<lower=0,upper=1>[T] p;

  real A1[T];
  real A2[T];
  real B1[T];
  real B2[T];
  int<lower=0,upper=1> y[T];
}

parameters {
  // non-centered
  vector[N] r_raw;
  vector[N] ll_raw;

  real mr;
  real<lower=0> sr;

  real mloglambda;
  real<lower=0> sloglambda;
}

transformed parameters {
  vector[N] r = mr + sr * r_raw;
  vector<lower=0>[N] lambda;

  for (i in 1:N) {
    lambda[i] = exp(mloglambda + sloglambda * ll_raw[i]);
  }
}

model {
  mr ~ normal(0.27, 0.36);
  sr ~ normal(0, 0.5) T[0,];

  mloglambda ~ normal(log(30), 0.5);
  sloglambda ~ normal(0, 0.5) T[0,];

  r_raw ~ normal(0, 1);
  ll_raw ~ normal(0, 1);

  for (t in 1:T) {
    real uA1 = crra_u(A1[t], r[uid[t]]);
    real uA2 = crra_u(A2[t], r[uid[t]]);
    real uB1 = crra_u(B1[t], r[uid[t]]);
    real uB2 = crra_u(B2[t], r[uid[t]]);

    real UA = p[t] * uA1 + (1 - p[t]) * uA2;
    real UB = p[t] * uB1 + (1 - p[t]) * uB2;

    // avoid p hitting exactly 0/1 numerically
    y[t] ~ bernoulli_logit(lambda[uid[t]] * (UA - UB));
  }
}
