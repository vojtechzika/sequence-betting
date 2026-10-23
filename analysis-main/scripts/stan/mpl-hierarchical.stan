// ------------------------------------------------------------
// Holt-Laury MPL: hierarchical model of the CRRA risk parameter r_i and
// the choice precision lambda_i
//
// Likelihood of participant i (as in mpl-bernoulli.stan):
//   y_t ~ Bernoulli(logistic(lambda_i * (UA_t - UB_t)))
//   UA_t - UB_t = p_t * d1_i + (1 - p_t) * d2_i,
//   d1_i = u(A1; r_i) - u(B1; r_i),  d2_i = u(A2; r_i) - u(B2; r_i)
//
// Population distribution (partial pooling):
//   r_i          ~ Normal(mu_r, sigma_r)
//   log lambda_i ~ Normal(mu_ll, sigma_ll)
// Hyperpriors:
//   mu_r      ~ Normal(h_mu_r_mean, h_mu_r_sd)
//   sigma_r   ~ HalfNormal(0, h_sigma_r_sd)
//   mu_ll     ~ Normal(log(h_lambda_mean), h_mu_ll_sd)
//   sigma_ll  ~ HalfNormal(0, h_sigma_ll_sd)
//
// Estimation. The participant parameters (r_i, log lambda_i) are
// integrated out on a fixed grid, so that the sampler explores only the four
// hyperparameters. The likelihood of a choice pattern on the grid does not
// depend on the hyperparameters and is passed as data (lik, scaled to a
// maximum of one within each pattern; the scale does not affect the
// posterior). For pattern j, the marginal likelihood is
//   m_j = sum_r f_r(r) * sum_l g_l(l) * lik[j, r, l],
// where f_r and g_l are the normal densities of the population distribution
// on the grid (times the grid step). Participants with the same choice pattern
// have the same marginal likelihood and enter through the pattern count.
// Posterior draws of r_i given the hyperparameters are generated in
// generated quantities from the exact conditional posterior on the r grid.
// The grid covers r in [r_grid[1], r_grid[R]]; the lower bounds on sigma_r
// and sigma_ll keep the population densities resolved by the grid step.
//
// lik is stored as a matrix with P blocks of R rows: row (j - 1) * R + k
// holds the likelihood of pattern j at r_grid[k], with columns indexing
// the log lambda grid.
// ------------------------------------------------------------
data {
  int<lower=1> N;                       // participants
  int<lower=1> P;                       // distinct choice patterns
  int<lower=2> R;                       // grid points for r
  int<lower=2> L;                       // grid points for log lambda
  int<lower=1,upper=P> pat[N];          // pattern of participant i
  int<lower=1> n_pat[P];                // number of participants with pattern j
  vector[R] r_grid;                     // grid for r (equally spaced)
  vector[L] ll_grid;                    // grid for log lambda (equally spaced)
  real<lower=0> dr;                     // step of the r grid
  real<lower=0> dl;                     // step of the log lambda grid
  matrix<lower=0>[R * P, L] lik;        // scaled likelihood of each pattern on the grid
  real h_mu_r_mean;
  real<lower=0> h_mu_r_sd;
  real<lower=0> h_sigma_r_sd;
  real<lower=0> h_lambda_mean;          // on natural scale; log taken in model
  real<lower=0> h_mu_ll_sd;
  real<lower=0> h_sigma_ll_sd;
}
parameters {
  real mu_r;
  real<lower=0.05> sigma_r;
  real mu_ll;
  real<lower=0.2> sigma_ll;
}
model {
  vector[R] fr;
  vector[L] gl;
  vector[R * P] t;

  mu_r     ~ normal(h_mu_r_mean, h_mu_r_sd);
  sigma_r  ~ normal(0, h_sigma_r_sd);
  mu_ll    ~ normal(log(h_lambda_mean), h_mu_ll_sd);
  sigma_ll ~ normal(0, h_sigma_ll_sd);

  // population densities on the grids
  fr = exp(-0.5 * square((r_grid - mu_r) / sigma_r)) / (sigma_r * sqrt(2 * pi())) * dr;
  gl = exp(-0.5 * square((ll_grid - mu_ll) / sigma_ll)) / (sigma_ll * sqrt(2 * pi())) * dl;

  // sum over log lambda for every pattern and r, then over r
  t = lik * gl;
  for (j in 1:P)
    target += n_pat[j] * log(dot_product(segment(t, (j - 1) * R + 1, R), fr));
}
generated quantities {
  vector[N] r_draw;                     // draw of r_i from its posterior given the hyperparameters
  {
    vector[R] fr;
    vector[L] gl;
    vector[R * P] t;
    fr = exp(-0.5 * square((r_grid - mu_r) / sigma_r)) / (sigma_r * sqrt(2 * pi())) * dr;
    gl = exp(-0.5 * square((ll_grid - mu_ll) / sigma_ll)) / (sigma_ll * sqrt(2 * pi())) * dl;
    t = lik * gl;
    for (i in 1:N) {
      vector[R] w;
      w = fr .* segment(t, (pat[i] - 1) * R + 1, R);
      w = w / sum(w);
      r_draw[i] = r_grid[categorical_rng(w)] + uniform_rng(-0.5 * dr, 0.5 * dr);
    }
  }
}
