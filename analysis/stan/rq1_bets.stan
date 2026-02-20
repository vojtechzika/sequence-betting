data {
  int<lower=1> N;                 // participants
  int<lower=1> S;                 // sequences (64)
  int<lower=1> T;                 // trials
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> y[T];      // 1=bet, 0=no bet
}

parameters {
  real alpha;                      // grand intercept
  vector[N] u_raw;                 // participant RE (non-centered)
  real<lower=0> sigma_u;

  vector[S] b_raw;                 // sequence effects (will be sum-to-zero)
  real<lower=0> sigma_s;
}

transformed parameters {
  vector[N] u = sigma_u * u_raw;

  // sum-to-zero constraint for sequence effects
  vector[S] b;
  b = b_raw - mean(b_raw);
}

model {
  // Priors per prereg (weakly informative)
  alpha   ~ normal(0, 1.5);
  sigma_u ~ normal(0, 1);          // HalfNormal(0,1) via <lower=0>
  sigma_s ~ normal(0, 1);          // HalfNormal(0,1) via <lower=0>
  u_raw   ~ normal(0, 1);
  b_raw   ~ normal(0, sigma_s);

  for (t in 1:T) {
    y[t] ~ bernoulli_logit(alpha + u[pid[t]] + b[sid[t]]);
  }
}

generated quantities {
  // sequence-level mean betting probability: mu_b_s = E_i[ logistic(alpha + u_i + b_s) ]
  vector[S] mu_b;
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b[s] = acc / N;
  }
}
