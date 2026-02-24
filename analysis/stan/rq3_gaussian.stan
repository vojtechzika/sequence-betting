// ------------------------------------------------------------
// RQ3: Certainty-equivalent welfare loss (Gaussian diagnostic)
// Outcome y_t = Δc_is / e (non-negative by construction; may include many zeros)
//
// y_t ~ Normal(mu_t, sigma)
// mu_t = alpha + u_i + b_s
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Priors:
//   alpha ~ Normal(0,1)
//   sigma, sigma_u, sigma_s ~ HalfNormal(0,1) via <lower=0>
//
// Generated quantities:
//   mu_c[s] = mean_i(alpha + u_i + b_s)
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector[T] y;
}

parameters {
  real alpha;

  vector[N] u_raw;
  real<lower=0> sigma_u;

  vector[S] b_raw;
  real<lower=0> sigma_s;

  real<lower=0> sigma;
}

transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = b_raw - mean(b_raw);
}

model {
  alpha   ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  sigma_s ~ normal(0, 1);
  sigma   ~ normal(0, 1);

  u_raw ~ normal(0, 1);
  b_raw ~ normal(0, sigma_s);

  for (t in 1:T) {
    y[t] ~ normal(alpha + u[pid[t]] + b[sid[t]], sigma);
  }
}

generated quantities {
  vector[S] mu_c;

  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      acc += alpha + u[i] + b[s];
    }
    mu_c[s] = acc / N;
  }
}
