data {
  int<lower=1> N;              // participants
  int<lower=1> S;              // sequences
  int<lower=1> T;              // betting trials only
  int<lower=1,upper=N> pid[T];
  int<lower=1,upper=S> sid[T];
  vector[T] z;                 // normalized stake deviations
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
  
  vector[S] b;
  b = b_raw - mean(b_raw);     // sum-to-zero
}

model {
  alpha   ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  sigma_s ~ normal(0, 1);
  sigma   ~ normal(0, 1);
  
  u_raw   ~ normal(0, 1);
  b_raw   ~ normal(0, sigma_s);
  
  z ~ normal(alpha + u[pid] + b[sid], sigma);
}