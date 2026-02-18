functions {
  real crra_u(real x, real r) {
    if (fabs(r - 1.0) < 1e-8) return log(x);
    else return pow(x, 1.0 - r) / (1.0 - r);
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
  vector[N] r;
  vector<lower=0>[N] lambda;
}

model {
  r ~ normal(0.27, 0.36);
  lambda ~ lognormal(log(30), 0.5);

  for (t in 1:T) {
    real UA = p[t] * crra_u(A1[t], r[uid[t]]) + (1 - p[t]) * crra_u(A2[t], r[uid[t]]);
    real UB = p[t] * crra_u(B1[t], r[uid[t]]) + (1 - p[t]) * crra_u(B2[t], r[uid[t]]);
    y[t] ~ bernoulli_logit(lambda[uid[t]] * (UA - UB));
  }
}

