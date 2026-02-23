// ------------------------------------------------------------
// RQ1 Betting Model (Hierarchical Logistic Regression)
//
// This model estimates sequence-level deviations in the probability
// of placing a bet (extensive margin) relative to a grand mean.
//
// Data structure:
//   b_is ∈ {0,1} indicates whether participant i placed a bet
//   in sequence s.
//
// Model:
//   b_is ~ Bernoulli(π_is)
//   logit(π_is) = α + u_i + β_s
//
// Hierarchical structure:
//   u_i  ~ Normal(0, σ_u)         // participant random effects
//   β_s  ~ Normal(0, σ_s)         // sequence effects
//
// Identification:
//   Sequence effects are constrained to sum to zero
//   (Σ_s β_s = 0), so α represents the population-average
//   log-odds of betting.
//
// Priors (weakly informative):
//   α        ~ Normal(0, 1.5)
//   σ_u, σ_s ~ HalfNormal(0, 1)
//
// Implementation note:
//   Random effects are implemented using a non-centered
//   parameterization for improved sampling efficiency.
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=0,upper=1> y[T];
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

  // sum-to-zero constraint for sequence effects (identified relative to grand mean)
  vector[S] b;
  b = b_raw - mean(b_raw);
}

model {
  alpha   ~ normal(0, 1.5);
  sigma_u ~ normal(0, 1);   // half-normal via <lower=0>
  sigma_s ~ normal(0, 1);   // half-normal via <lower=0>

  u_raw ~ normal(0, 1);

  // sequence effects on sigma_s scale, then sum-to-zero in transformed parameters
  b_raw ~ normal(0, sigma_s);

  for (t in 1:T)
    y[t] ~ bernoulli_logit(alpha + u[pid[t]] + b[sid[t]]);
}

generated quantities {
  vector[S] mu_b;
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) acc += inv_logit(alpha + u[i] + b[s]);
    mu_b[s] = acc / N;
  }
}
