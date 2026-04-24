// ------------------------------------------------------------
// RQ2: Proportional stake deviation on betting trials
//
// Observed (per trial t with stake>0):
//   delta_t = stake_t - a_star_{i,treat}(r_i)   (in ECU)
//   z_t = (delta_t - delta_bar_i) / s_star_i    (standardized within participant)
//
// Model on z_t:
//   z_t ~ Normal(alpha + u_i + b_s [+ gamma_drift * block_c_t], sigma)
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)
//   b_s ~ Normal(0, sigma_s) with sum-to-zero via centering
//
// Optional linear drift (include_drift = 1):
//   gamma_drift * block_c_t
//   block_c in {-1.5, -0.5, 0.5, 1.5} (centered block index)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//   When include_drift = 0, gamma_drift is sampled but not used;
//   prior shrinks it toward zero.
//
// Priors:
//   alpha       ~ Normal(0, 1)
//   sigma_u     ~ HalfNormal(0, 1)
//   sigma_s     ~ HalfNormal(0, 1)
//   sigma       ~ HalfNormal(0, 1)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//
// Data inputs:
//   use_bb         : reserved for BB robustness (not yet implemented)
//   include_drift  : 0/1 flag for drift adjustment
//   block_c        : centered block index for each observation
//   prior_gamma_sd : prior SD for drift coefficient
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
  vector[T] z;                       // standardized stake deviation
  vector[N] delta_bar;               // participant mean deviation (ECU)
  vector[N] s_star;                  // participant SD of deviation (ECU)
  int<lower=1> e;                    // endowment (ECU) for proportional scaling
  int<lower=0,upper=1> include_drift; // 1 = include linear drift, 0 = no drift
  vector[T] block_c;                 // centered block index {-1.5,-0.5,0.5,1.5}
  real<lower=0> prior_gamma_sd;      // prior SD for drift coefficient
}
parameters {
  real alpha;
  vector[N] u_raw;
  real<lower=0> sigma_u;
  vector[S] b_raw;
  real<lower=0> sigma_s;
  real<lower=0> sigma;
  real gamma_drift;                  // always sampled; prior shrinks to 0 when not needed
}
transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
}
model {
  // priors
  alpha       ~ normal(0, 1);
  sigma_u     ~ normal(0, 1);
  sigma_s     ~ normal(0, 1);
  sigma       ~ normal(0, 1);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  gamma_drift ~ normal(0, prior_gamma_sd);

  // likelihood
  {
    vector[T] eta = alpha + u[pid] + b[sid];
    if (include_drift) eta += gamma_drift * block_c;
    z ~ normal(eta, sigma);
  }
}
generated quantities {
  vector[S] mu_a;
  vector[N] mu_a_i;

  // per-sequence: population mean proportional deviation (at block_c = 0)
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      acc += (delta_bar[i] + (alpha + u[i] + b[s]) * s_star[i]) / e;
    }
    mu_a[s] = acc / N;
  }

  // per-participant: mean proportional deviation (sequence effects averaged out)
  for (i in 1:N) {
    mu_a_i[i] = (delta_bar[i] + (alpha + u[i]) * s_star[i]) / e;
  }
}