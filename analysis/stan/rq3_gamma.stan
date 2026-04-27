// ------------------------------------------------------------
  // RQ3: Certainty-equivalent welfare loss (Gamma-only)
// y_t = Δc_is / e, strictly positive
//
  // Used when hurdle probability is negligible (mean_omega < 0.05),
// i.e. essentially all welfare losses are strictly positive and
// the hurdle component is not warranted.
//
  // positive part only:
  //   y | y > 0 ~ Gamma(shape, rate = shape / mu_is)
//   log(mu_is) = ap + up_i + bp_s [+ gamma_drift * block_c_t]
//
  // Random effects (sum-to-zero sequence effects):
  //   up_i ~ Normal(0, sup)
//   bp_s ~ Normal(0, sbp) centered
//
  // Priors:
  //   ap ~ Normal(0, 1)
//   sup, sbp ~ HalfNormal(0, 1)
//   shape ~ Gamma(2, 0.1)
//   gamma_drift ~ Normal(0, prior_gamma_sd)
//
  // Generated quantities:
  //   mu_c[s]   = mean_i E[y_{is} | params]
  //   mu_c_i[i] = mean_s E[y_{is} | params]
  // ------------------------------------------------------------
    
    data {
      int<lower=1> N;
      int<lower=1> S;
      int<lower=1> T;
      int<lower=1,upper=N> pid[T];
      int<lower=1,upper=S> sid[T];
      vector<lower=0>[T] y;
      int<lower=0,upper=1> include_drift;
      vector[T] block_c;
      real<lower=0> prior_gamma_sd;
    }
  parameters {
    real ap;
    vector[N] up_raw;
    real<lower=0> sup;
    vector[S] bp_raw;
    real<lower=0> sbp;
    real<lower=0> shape;
    real gamma_drift;
  }
  transformed parameters {
    vector[N] up = sup * up_raw;
    vector[S] bp = sbp * (bp_raw - mean(bp_raw));
  }
  model {
    // priors
    ap          ~ normal(0, 1);
    sup         ~ normal(0, 1);
    sbp         ~ normal(0, 1);
    up_raw      ~ normal(0, 1);
    bp_raw      ~ normal(0, 1);
    shape       ~ gamma(2, 0.1);
    gamma_drift ~ normal(0, prior_gamma_sd);
    
    // likelihood
    for (t in 1:T) {
      real eta = ap + up[pid[t]] + bp[sid[t]];
      if (include_drift) eta += gamma_drift * block_c[t];
      real mu_pos = exp(eta);
      target += gamma_lpdf(y[t] | shape, shape / mu_pos);
    }
  }
  generated quantities {
    vector[S] mu_c;
    vector[N] mu_c_i;
    
    // per-sequence population mean: mean_i E[y_is]
    for (s in 1:S) {
      real acc = 0;
      for (i in 1:N) {
        real eta = ap + up[i] + bp[s];
        acc += exp(eta);
      }
      mu_c[s] = acc / N;
    }
    
    // per-participant mean across sequences: mean_s E[y_is]
    for (i in 1:N) {
      real acc = 0;
      for (s in 1:S) {
        real eta = ap + up[i] + bp[s];
        acc += exp(eta);
      }
      mu_c_i[i] = acc / S;
    }
  }