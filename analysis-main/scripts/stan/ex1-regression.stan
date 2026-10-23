// EX1.2: Predictors of GHI, one regression per posterior draw of chi_i
//
// For a given draw t of the participant-level directional index,
//   chi_i^(t) = alpha + b_opt * Zopt_i + b_rt * Zrt_i + b_r * Zr_i + eps_i,
//   eps_i ~ Normal(0, sigma).
// The script 54_ex1-stan-regression.R fits this model separately for each of
// Trep draws and pools the posteriors with equal weights, as specified in
// the preregistration (EX1, part II).
//
// Priors (preregistered):
//   alpha      ~ Normal(0, prior_alpha_sd)
//   b_*        ~ Normal(0, prior_b_sd)
//   sigma      ~ HalfNormal(0, prior_sigma_sd)

data {
  int<lower=1> N;        // participants
  vector[N] y;           // chi_i^(t), one posterior draw
  vector[N] Zopt;        // standardized LOT-R
  vector[N] Zrt;         // standardized log response time
  vector[N] Zr;          // standardized risk parameter r
  real<lower=0> prior_alpha_sd;  // prior SD of the intercept
  real<lower=0> prior_b_sd;      // prior SD of the slopes
  real<lower=0> prior_sigma_sd;  // prior scale of the residual SD (half-normal)
}

parameters {
  real alpha;
  real b_opt;
  real b_rt;
  real b_r;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, prior_alpha_sd);
  b_opt ~ normal(0, prior_b_sd);
  b_rt  ~ normal(0, prior_b_sd);
  b_r   ~ normal(0, prior_b_sd);
  sigma ~ normal(0, prior_sigma_sd);   // half-normal through the lower bound

  y ~ normal(alpha + b_opt * Zopt + b_rt * Zrt + b_r * Zr, sigma);
}
