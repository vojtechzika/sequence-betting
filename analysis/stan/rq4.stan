// ------------------------------------------------------------
// rq4.stan
//
// RQ4: Side choices conditional on betting (Bernoulli-logit)
//
// MODEL
//   h_t ~ Bernoulli(pi_is)
//   logit(pi_is) = alpha + u_i + beta_s [+ gamma_drift * block_c_t]
//
// Random effects (non-centred parameterisation):
//   u_i    ~ Normal(0, sigma_u)   participant intercepts
//   beta_s ~ Normal(0, sigma_s)   sequence effects, sum-to-zero via
//                                 centering in transformed parameters
//
// Drift (optional, activated by include_drift == 1):
//   Linear:      gamma_drift * block_c_t,  block_c_t in {-1.5,-0.5,0.5,1.5}
//   Categorical: gamma_block[k],  k = 1..4, sum-to-zero constrained
//   Which variant is selected by drift_type:
//     0 = no drift
//     1 = linear drift
//     2 = categorical drift
//
// Priors:
//   alpha              ~ Normal(0, 1.5)
//   sigma_u, sigma_s   ~ HalfNormal(0, 1)
//   gamma_drift        ~ Normal(0, prior_gamma_sd)
//   gamma_block[k]     ~ Normal(0, prior_gamma_sd), sum-to-zero
//
// Generated quantities:
//   mu_h[s]       population mean per sequence: E_{u~N(0,sigma_u)} logistic(alpha + u + beta_s)
//   hbar          population baseline:          E_{u~N(0,sigma_u)} logistic(alpha + u)
//   mu_h_i[i]     participant baseline:         logistic(alpha + u_i)
//                 NOTE: under the sum-to-zero constraint on beta_s,
//                 E_s[beta_s] = 0 by construction, so mu_h_i[i] equals
//                 the sequence-averaged participant tendency
//                 E_s[E(pi_is | data)]. If the sum-to-zero constraint were
//                 relaxed this equivalence would not hold.
//   mu_h_sample[s]  sample mean over realised u_i (exact given posterior u)
//   hbar_sample     sample baseline (exact given posterior u)
//   y_rep[T]      posterior predictive replications (for PPC / overdispersion
//                 check; sequence-wise variance of y_rep vs observed used to
//                 detect overdispersion that would trigger Beta-Binomial
//                 robustness fit)
//
// MC integration for population quantities:
//   M draws from N(0, sigma_u) are used per posterior draw.
//   M is passed as data so it can be set from cfg without recompilation.
//   Recommended default: M = 2000.
//
// INPUTS (from 41_rq4_stan.R):
//   N, S, T, pid[T], sid[T], h[T]
//   include_drift, drift_type, block_c[T], prior_gamma_sd
//   mc_M
//
// OUTPUTS consumed downstream:
//   mu_h[S], hbar, mu_h_i[N]   -> rq4_tables.R (sequence/participant summaries)
//   y_rep[T]                   -> rq4_diagnostics.R (PPC, overdispersion test)
//   beta[S]                    -> rq4_tables.R (relative sequence contrasts)
// ------------------------------------------------------------

data {
  int<lower=1> N;                       // participants
  int<lower=1> S;                       // sequences
  int<lower=1> T;                       // trials (betting trials only)
  array[T] int<lower=1,upper=N> pid;
  array[T] int<lower=1,upper=S> sid;
  array[T] int<lower=0,upper=1>  h;    // 1 = Heads, 0 = Tails

  // Drift control
  int<lower=0,upper=1> include_drift;   // 0 = no drift, 1 = include drift
  int<lower=0,upper=2> drift_type;      // 0 = none, 1 = linear, 2 = categorical
  vector[T] block_c;                    // centred block index; 0 when include_drift == 0
  real<lower=0> prior_gamma_sd;         // prior SD for drift parameters (e.g. 0.3)

  // MC integration resolution for population quantities
  int<lower=100> mc_M;                  // number of MC draws; recommended 2000
}

parameters {
  real alpha;

  // Participant random effects (non-centred)
  vector[N] u_raw;
  real<lower=0> sigma_u;

  // Sequence random effects (non-centred; sum-to-zero applied in transformed parameters)
  vector[S] b_raw;
  real<lower=0> sigma_s;

  // Drift parameters (active only when include_drift == 1)
  // Linear drift coefficient
  real gamma_drift;
  // Categorical drift (4 blocks); identified via sum-to-zero in transformed parameters
  vector[4] gamma_block_raw;
}

transformed parameters {
  // Participant effects
  vector[N] u = sigma_u * u_raw;

  // Sequence effects: non-centred + sum-to-zero
  vector[S] beta = sigma_s * (b_raw - mean(b_raw));

  // Categorical drift: sum-to-zero centering
  vector[4] gamma_block = gamma_block_raw - mean(gamma_block_raw);
}

model {
  // ---- Priors ----
  alpha   ~ normal(0, 1.5);
  sigma_u ~ normal(0, 1);
  sigma_s ~ normal(0, 1);
  u_raw   ~ normal(0, 1);
  b_raw   ~ normal(0, 1);

  // Drift priors: always placed so the sampler geometry is stable;
  // the likelihood contribution is zeroed out when include_drift == 0
  // via the multiplier below.
  gamma_drift       ~ normal(0, prior_gamma_sd);
  gamma_block_raw   ~ normal(0, prior_gamma_sd);

  // ---- Likelihood ----
  for (t in 1:T) {
    real eta = alpha + u[pid[t]] + beta[sid[t]];

    if (include_drift == 1) {
      if (drift_type == 1) {
        // Linear drift
        eta += gamma_drift * block_c[t];
      } else if (drift_type == 2) {
        // Categorical drift: block index recovered from block_c
        // block_c in {-1.5, -0.5, 0.5, 1.5} -> block 1..4
        int blk = to_int(round(block_c[t] + 2.5));  // maps to 1,2,3,4
        eta += gamma_block[blk];
      }
    }

    h[t] ~ bernoulli_logit(eta);
  }
}

generated quantities {
  // ------------------------------------------------------------------
  // 1. Posterior predictive replications (for PPC and overdispersion
  //    check; sequence-wise sum of y_rep used in rq4_diagnostics.R)
  // ------------------------------------------------------------------
  array[T] int y_rep;
  for (t in 1:T) {
    real eta = alpha + u[pid[t]] + beta[sid[t]];
    if (include_drift == 1) {
      if (drift_type == 1) {
        eta += gamma_drift * block_c[t];
      } else if (drift_type == 2) {
        int blk = to_int(round(block_c[t] + 2.5));
        eta += gamma_block[blk];
      }
    }
    y_rep[t] = bernoulli_logit_rng(eta);
  }

  // ------------------------------------------------------------------
  // 2. Sample-based summaries (exact given posterior u_i)
  //    Useful for traceability and cross-checking population quantities.
  // ------------------------------------------------------------------
  vector[S] mu_h_sample;
  real      hbar_sample;
  vector[N] mu_h_i;        // participant baselines: logistic(alpha + u_i)
                            // equals sequence-averaged tendency E_s[pi_is]
                            // under the sum-to-zero constraint on beta_s

  for (i in 1:N)
    mu_h_i[i] = inv_logit(alpha + u[i]);

  {
    real acc0 = 0.0;
    for (i in 1:N) acc0 += inv_logit(alpha + u[i]);
    hbar_sample = acc0 / N;

    for (s in 1:S) {
      real accs = 0.0;
      for (i in 1:N) accs += inv_logit(alpha + u[i] + beta[s]);
      mu_h_sample[s] = accs / N;
    }
  }

  // ------------------------------------------------------------------
  // 3. Population-based summaries via MC integration
  //    Integrates over the full participant heterogeneity distribution
  //    N(0, sigma_u) rather than conditioning on realised u_i.
  //    Uses mc_M fresh draws with common random numbers across sequences
  //    so that mu_h[s] and hbar are on the same population scale and
  //    their difference is a well-defined estimand for sequence bias.
  //    (Post-pilot addendum clarification, Section A.5.)
  // ------------------------------------------------------------------
  vector[S] mu_h;
  real      hbar;

  {
    array[mc_M] real u_new;
    real acc0 = 0.0;
    for (m in 1:mc_M) {
      u_new[m] = normal_rng(0, sigma_u);
      acc0     += inv_logit(alpha + u_new[m]);
    }
    hbar = acc0 / mc_M;

    for (s in 1:S) {
      real accs = 0.0;
      for (m in 1:mc_M)
        accs += inv_logit(alpha + u_new[m] + beta[s]);
      mu_h[s] = accs / mc_M;
    }
  }
}
