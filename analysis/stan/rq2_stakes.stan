// ------------------------------------------------------------
// RQ2 Staking Model (Hierarchical Normal on Z^{Δa})
//
// Outcome:
//   z_is = standardized within-participant stake deviation
//          Z^{Δa}_{is} = (Δa_is - mean_i(Δa)) / s_i*(Δa)
//
// Model:
//   z_is ~ Normal(η_is, σ)
//   η_is = α + u_i + b_s
//
// Random effects:
//   u_i ~ Normal(0, σ_u)
//   b_s ~ Normal(0, σ_s) with sum-to-zero constraint Σ_s b_s = 0
//
// Priors (weakly informative):
//   α ~ Normal(0, 1)
//   σ, σ_u, σ_s ~ HalfNormal(0, 1)
//
// Generated quantities:
//   mu_a[s] = (1/e) * mean_i[ Δa_bar_i + (α + u_i + b_s) * s_star_i ]
//   (descriptive calibration contrast on absolute scale)
// ------------------------------------------------------------

data {
  int<lower=1> N;                 // participants
  int<lower=1> S;                 // sequences
  int<lower=1> T;                 // betting trials used in RQ2
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];

  vector[T] z;                    // Z^{Δa}_{is}

  // participant-level quantities used for mapping back to absolute scale
  vector[N] delta_bar;            // \bar{Δa}_i (computed across betting trials)
  vector<lower=0>[N] s_star;      // s_i*(Δa) = max(sd_i(Δa), floor)
  real<lower=1> e;                // endowment (e.g., 100)
}

parameters {
  real alpha;

  vector[N] u_raw;
  real<lower=0> sigma_u;

  vector[S] b_raw;               // sequence effects (will be centered to sum-to-zero)
  real<lower=0> sigma_s;

  real<lower=0> sigma;           // residual SD
}

transformed parameters {
  vector[N] u = sigma_u * u_raw;

  // sum-to-zero constraint for sequence effects (identified relative to grand mean)
  vector[S] b;
  b = b_raw - mean(b_raw);
}

model {
  // Priors
  alpha   ~ normal(0, 1);
  sigma_u ~ normal(0, 1);         // half-normal via <lower=0>
  sigma_s ~ normal(0, 1);         // half-normal via <lower=0>
  sigma   ~ normal(0, 1);         // half-normal via <lower=0>

  u_raw ~ normal(0, 1);

  // sequence effects with scale sigma_s, then sum-to-zero in transformed parameters
  b_raw ~ normal(0, sigma_s);

  // Likelihood
  for (t in 1:T) {
    z[t] ~ normal(alpha + u[pid[t]] + b[sid[t]], sigma);
  }
}

generated quantities {
  // Descriptive back-mapped sequence contrast on absolute scale (as fraction of endowment)
  vector[S] mu_a;

  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      acc += delta_bar[i] + (alpha + u[i] + b[s]) * s_star[i];
    }
    mu_a[s] = (acc / N) / e;
  }
}
