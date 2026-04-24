// ------------------------------------------------------------
// RQ2 Robustness: One-inflated Beta-Binomial on raw stake scale
//
// Triggered when PPC detects boundary inflation from all-in
// stakes (a = e). Models the stake distribution as a mixture:
//   - Point mass at a = e (all-in) with probability omega_i
//   - Beta-Binomial on {1,...,e-1} otherwise
//
// The expected stake under this model is:
//   E[a_ist] = omega_i * e + (1 - omega_i) * mu_ist * e
//   where mu_ist = inv_logit(alpha + u_i + b_s [+ gamma_drift * block_c_t])
//
// Generated quantities express mu_a[s] and mu_a_i[i] as mean
// proportional deviations from EU-optimal stake a_star[i],
// directly comparable to the primary Gaussian model outputs.
//
// Linear predictor for mean stake (on logit scale):
//   logit(mu_ist) = alpha + u_i + b_s [+ gamma_drift * block_c_t]
//
// One-inflation probability (participant-level):
//   logit(omega_i) = alpha_omega + v_i
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)  -- participant stake level
//   b_s ~ Normal(0, sigma_s)  -- sequence effect, sum-to-zero
//   v_i ~ Normal(0, sigma_v)  -- participant all-in propensity
//
// Priors:
//   alpha, alpha_omega ~ Normal(0, 1.5)
//   sigma_u, sigma_s, sigma_v ~ HalfNormal(0, 1)
//   phi ~ HalfNormal(0, 5)    -- Beta-Binomial overdispersion
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//
// Data inputs:
//   y              : raw stake (1..e) per trial
//   a_star         : posterior mean EU-optimal stake per participant
//   include_drift  : 0/1 flag for linear drift
//   block_c        : centered block index {-1.5,-0.5,0.5,1.5}
//   prior_gamma_sd : prior SD for drift coefficient
//
// Generated quantities (comparable to primary Gaussian model):
//   mu_a[s]   = mean_i[ (E[a_ist] - a_star[i]) / e ]
//   mu_a_i[i] = (E[a_i] - a_star[i]) / e
// ------------------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=1> y[T];                  // raw stake (1..e)
  vector[N] a_star;                   // posterior mean EU-optimal stake per participant
  int<lower=1> e;                     // endowment (ECU)
  int<lower=0,upper=1> include_drift;
  vector[T] block_c;
  real<lower=0> prior_gamma_sd;
}
parameters {
  real alpha;
  vector[N] u_raw;
  real<lower=0> sigma_u;
  vector[S] b_raw;
  real<lower=0> sigma_s;
  real<lower=0> phi;
  real alpha_omega;
  vector[N] v_raw;
  real<lower=0> sigma_v;
  real gamma_drift;
}
transformed parameters {
  vector[N] u = sigma_u * u_raw;
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
  vector[N] v = sigma_v * v_raw;
}
model {
  // priors
  alpha       ~ normal(0, 1.5);
  sigma_u     ~ normal(0, 1);
  sigma_s     ~ normal(0, 1);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  phi         ~ normal(0, 5);
  alpha_omega ~ normal(0, 1.5);
  sigma_v     ~ normal(0, 1);
  v_raw       ~ normal(0, 1);
  gamma_drift ~ normal(0, prior_gamma_sd);

  // likelihood
  for (t in 1:T) {
    real eta   = alpha + u[pid[t]] + b[sid[t]];
    if (include_drift) eta += gamma_drift * block_c[t];
    real mu    = inv_logit(eta);
    real omega = inv_logit(alpha_omega + v[pid[t]]);

    if (y[t] == e) {
      // all-in stake: mixture of one-inflation and BB mass at upper boundary
      target += log_mix(
        omega,
        0.0,
        beta_binomial_lpmf(e - 1 | e - 1, mu * phi, (1 - mu) * phi)
      );
    } else {
      // interior stake: only BB component
      target += log1m(omega) +
                beta_binomial_lpmf(y[t] - 1 | e - 1, mu * phi, (1 - mu) * phi);
    }
  }
}
generated quantities {
  vector[S] mu_a;
  vector[N] mu_a_i;

  // per-sequence population mean proportional deviation from EU-optimal
  // E[a_ist] = omega_i * e + (1 - omega_i) * mu_ist * e
  for (s in 1:S) {
    real acc = 0;
    for (i in 1:N) {
      real mu_ist  = inv_logit(alpha + u[i] + b[s]);
      real omega_i = inv_logit(alpha_omega + v[i]);
      real E_a     = omega_i * e + (1 - omega_i) * mu_ist * e;
      acc += (E_a - a_star[i]) / e;
    }
    mu_a[s] = acc / N;
  }

  // per-participant mean proportional deviation from EU-optimal
  // (sequence effects averaged out by sum-to-zero constraint on b)
  for (i in 1:N) {
    real mu_i    = inv_logit(alpha + u[i]);
    real omega_i = inv_logit(alpha_omega + v[i]);
    real E_a     = omega_i * e + (1 - omega_i) * mu_i * e;
    mu_a_i[i]   = (E_a - a_star[i]) / e;
  }
}