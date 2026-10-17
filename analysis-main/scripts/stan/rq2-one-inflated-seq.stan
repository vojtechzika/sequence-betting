// ------------------------------------------------------------
// RQ2 (extension of the preregistered alternative): one-inflated
// Beta-Binomial on the raw stake scale in which the all-in probability
// also varies by sequence; see 23_rq2-stan-one-inflated-seq.R.
//
// The model of rq2-one-inflated.stan lets sequences shift only the stakes
// that are not all-in, while the all-in probability varies by participant
// alone. If sequences act mainly on the decision to stake the whole
// endowment, that model attributes the effect to a channel in which it does
// not occur and shrinks the sequence contrasts. Here the all-in probability
// has a sequence effect c_s of the same form as b_s; with sigma_c = 0 the
// model reduces to rq2-one-inflated.stan.
//
// Stake y_t in {1,...,e} on betting trial t of participant i = pid[t]
// and sequence s = sid[t]. With probability omega_t the stake is e
// (one-inflation); otherwise y_t - 1 follows a Beta-Binomial with e - 1
// trials, so that y_t is supported on {1,...,e}:
//   y_t - 1 ~ BetaBinomial(e - 1, mu_t * phi, (1 - mu_t) * phi)
//   logit(mu_t)    = alpha + u_i + b_s [+ gamma_drift * block_c_t]
//   logit(omega_t) = alpha_omega + v_i + c_s
// An all-in stake therefore has probability
//   omega_t + (1 - omega_t) * BetaBinomial(e - 1 | e - 1, ...).
//
// Random effects:
//   u_i ~ Normal(0, sigma_u)   participant stake level, centred
//   b_s ~ Normal(0, sigma_s)   sequence effect on the stake, centred
//   v_i ~ Normal(0, sigma_v)   participant all-in propensity
//   c_s ~ Normal(0, sigma_c)   sequence effect on the all-in probability, centred
// (u and b are centred so that they sum to zero; likewise c. The
// intercepts alpha and alpha_omega then carry the means.)
//
// Priors (scales are passed as data, set in settings.R):
//   alpha       ~ Normal(0, prior_alpha_sd)
//   sigma_u     ~ HalfNormal(0, prior_sigma_u_sd)
//   sigma_s     ~ HalfNormal(0, prior_sigma_s_sd)
//   sigma_v     ~ HalfNormal(0, prior_sigma_v_sd)
//   sigma_c     ~ HalfNormal(0, prior_sigma_c_sd)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//   alpha_omega ~ Normal(0, prior_alpha_omega_sd)   logit scale
//   phi         ~ HalfNormal(0, prior_phi_sd)       Beta-Binomial concentration
//
// Data inputs:
//   y              : stake (1..e) per betting trial
//   a_star         : posterior mean EU-optimal stake per participant
//   include_drift  : 0/1 flag for linear drift
//   block_c        : centered block index {-1.5,-0.5,0.5,1.5}
//   prior_*        : prior scales listed above
//
// Generated quantities (same estimands as the other RQ2 models, at block_c = 0):
//   E[y_is]   = omega_is * e + (1 - omega_is) * (1 + (e - 1) * mu_is)
//   mu_a[s]   = mean_i( (E[y_is] - a_star_i) / e )
//   mu_a_i[i] = mean_s( (E[y_is] - a_star_i) / e )
// ------------------------------------------------------------

data {
  int<lower=2> e;                         // endowment (ECU)
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  int<lower=1,upper=e> y[T];              // stake
  vector[N] a_star;
  int<lower=0,upper=1> include_drift;
  vector[T] block_c;
  real<lower=0> prior_gamma_sd;           // prior SD for the drift coefficient
  real<lower=0> prior_alpha_sd;           // prior SD for the intercept
  real<lower=0> prior_sigma_u_sd;         // prior scale for sigma_u (half-normal)
  real<lower=0> prior_sigma_s_sd;         // prior scale for sigma_s (half-normal)
  real<lower=0> prior_sigma_v_sd;         // prior scale for sigma_v (half-normal)
  real<lower=0> prior_sigma_c_sd;         // prior scale for sigma_c (half-normal)
  real<lower=0> prior_alpha_omega_sd;     // prior SD for the all-in intercept
  real<lower=0> prior_phi_sd;             // prior scale for phi (half-normal)
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
  vector[S] c_raw;
  real<lower=0> sigma_c;
  real gamma_drift;                       // always sampled; prior shrinks to 0 when not used
}
transformed parameters {
  vector[N] u = sigma_u * (u_raw - mean(u_raw));
  vector[S] b = sigma_s * (b_raw - mean(b_raw));
  vector[N] v = sigma_v * v_raw;
  vector[S] c = sigma_c * (c_raw - mean(c_raw));
}
model {
  alpha       ~ normal(0, prior_alpha_sd);
  sigma_u     ~ normal(0, prior_sigma_u_sd);
  sigma_s     ~ normal(0, prior_sigma_s_sd);
  sigma_v     ~ normal(0, prior_sigma_v_sd);
  sigma_c     ~ normal(0, prior_sigma_c_sd);
  u_raw       ~ normal(0, 1);
  b_raw       ~ normal(0, 1);
  v_raw       ~ normal(0, 1);
  c_raw       ~ normal(0, 1);
  phi         ~ normal(0, prior_phi_sd);
  alpha_omega ~ normal(0, prior_alpha_omega_sd);
  gamma_drift ~ normal(0, prior_gamma_sd);

  for (t in 1:T) {
    real eta = alpha + u[pid[t]] + b[sid[t]];
    if (include_drift) eta += gamma_drift * block_c[t];
    {
      real mu    = inv_logit(eta);
      real omega = inv_logit(alpha_omega + v[pid[t]] + c[sid[t]]);
      real lbb   = beta_binomial_lpmf(y[t] - 1 | e - 1, mu * phi, (1 - mu) * phi);
      if (y[t] == e) {
        // all-in: one-inflation mass plus the Beta-Binomial mass at the upper bound
        target += log_mix(omega, 0.0, lbb);
      } else {
        target += log1m(omega) + lbb;
      }
    }
  }
}
generated quantities {
  vector[S] mu_a;
  vector[N] mu_a_i;
  {
    matrix[N, S] dev;
    for (i in 1:N) {
      for (s in 1:S) {
        real omega_is = inv_logit(alpha_omega + v[i] + c[s]);
        real mu_is    = inv_logit(alpha + u[i] + b[s]);
        real E_y      = omega_is * e + (1 - omega_is) * (1 + (e - 1) * mu_is);
        dev[i, s]     = (E_y - a_star[i]) / e;
      }
    }
    for (s in 1:S) mu_a[s]   = mean(col(dev, s));
    for (i in 1:N) mu_a_i[i] = mean(row(dev, i));
  }
}
