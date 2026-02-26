// ------------------------------------------------------------
// RQ2: Proportional stake deviation on betting trials (Normal)
//
// Observed (per trial t with stake>0):
//   delta_t = stake_t - a_star_{i,treat}(r_i)   (in ECU)
//   z_t = (delta_t - delta_bar_i) / s_star_i    (standardized within participant)
//
// Model on z_t:
//   z_t ~ Normal(alpha + u_i + b_s, sigma)
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Priors:
//   alpha ~ Normal(0,1.5)
//   sigma, sigma_u, sigma_s ~ HalfNormal(0,1) via <lower=0>
//
// Generated quantities:
//   mu_a[s]   = mean_i( delta_bar_i + (alpha + u_i + b_s) * s_star_i ) / e
//   mu_a_i[i] = ( delta_bar_i + (alpha + u_i) * s_star_i ) / e
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector[T] z;                 // standardized stake deviation
  vector[N] delta_bar;         // participant mean deviation (ECU)
  vector[N] s_star;            // participant sd of deviation (ECU)
  int<lower=1> e;              // endowment (ECU) for proportional scaling
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

  // sum-to-zero sequence effects with hierarchical scale
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}

model {
  // priors
  alpha   ~ normal(0, 1.5);
  sigma_u ~ normal(0, 1);
  sigma_s ~ normal(0, 1);
  sigma   ~ normal(0, 1);

  u_raw ~ normal(0, 1);
  b_raw ~ normal(0, 1);

  // likelihood
  for (t in 1:T) {
    z[t] ~ normal(alpha + u[pid[t]] + b[sid[t]], sigma);
  }
}

generated quantities {
  vector[S] mu_a;
  vector[N] mu_a_i;

  // per-sequence: population mean proportional deviation
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      acc += (delta_bar[i] + (alpha + u[i] + b[s]) * s_star[i]) / e;
    }
    mu_a[s] = acc / N;
  }

  // per-participant: mean proportional deviation (sequence effects averaged out by centering)
  for (i in 1:N) {
    mu_a_i[i] = (delta_bar[i] + (alpha + u[i]) * s_star[i]) / e;
  }
}
