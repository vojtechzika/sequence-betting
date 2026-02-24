data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> h[T];   // 1=Heads, 0=Tails
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

  int J = 2000;           // MC draws per iteration (tune)
  real acc;

  // baseline (beta = 0)
  acc = 0;
  for (j in 1:J) {
    real u0 = normal_rng(0, sigma_u);
    acc += inv_logit(alpha + u0);
  }
  hbar = acc / J;

  // sequence means
  for (s in 1:S) {
    acc = 0;
    for (j in 1:J) {
      real u0 = normal_rng(0, sigma_u);
      acc += inv_logit(alpha + u0 + beta[s]);
    }
    mu_h[s] = acc / J;
  }
}
