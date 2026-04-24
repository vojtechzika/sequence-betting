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
  real prior_r_mean;
  real<lower=0> prior_r_sd;
  real<lower=0> prior_lambda_mean;   // on natural scale; log taken in model
  real<lower=0> prior_lambda_sd;     // on log scale
}
parameters {
  vector[N] r;
  vector<lower=0>[N] lambda;
}
model {
  r ~ normal(prior_r_mean, prior_r_sd);
  lambda ~ lognormal(log(prior_lambda_mean), prior_lambda_sd);
  for (t in 1:T) {
    real UA = p[t] * crra_u(A1[t], r[uid[t]]) + (1 - p[t]) * crra_u(A2[t], r[uid[t]]);
    real UB = p[t] * crra_u(B1[t], r[uid[t]]) + (1 - p[t]) * crra_u(B2[t], r[uid[t]]);
    y[t] ~ bernoulli_logit(lambda[uid[t]] * (UA - UB));
  }
}