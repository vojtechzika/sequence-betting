// ============================================================
// synthesis.stan
//
// Joint hierarchical model of the extensive margin (bet vs. no-bet)
// and the side-choice decision conditional on betting (H vs. T).
// Tests whether the terminal-run and imbalance cues jointly predict
// both decisions, and whether betting and side choice share one
// participant-level tendency or are "dissociated."
//
// Cue covariates (outcome-blind, sequence-level; built in
// synthesis_stan(), 96_synthesis_stan.R):
//   -- betting equation (magnitude / presence only, unsigned) --
//   imb_mag = imbalance magnitude, |H count - T count| / 6, in [0,1]
//   run_dum = 1 if seq ends in a run of >=2 identical symbols, else 0
//   -- side-choice equation (signed) --
//   run_h_s = signed terminal-run cue (master_sequences.csv: terminal_run)
//             +1 if seq ends in a run of >=2 H's, -1 if >=2 T's, 0 if none
//   imb_h_s = signed imbalance cue = (H count - T count)/6
//             (master_sequences.csv: imbalance)
//
// Both equations get participant-level random intercepts + slopes on
// their own two regressors. The 6 participant-level random effects (3
// per equation: [bet_intercept, bet_run, bet_imb, side_intercept,
// side_run, side_imb]) share one correlated multivariate-normal prior
// (LKJ), giving the cross-equation correlations that formalize the
// dissociation hypothesis:
//   Omega[1,4] = corr_intercepts : overall bet-propensity <-> side-bias
//   Omega[2,5] = corr_run_slopes : run-cue sensitivity,  bet <-> side
//   Omega[3,6] = corr_imb_slopes : imbalance-cue sensitivity, bet <-> side
// Omega[2,5]/Omega[3,6] compare engagement-sensitivity (bet eq) to
// direction-sensitivity (side eq), not the same cue in both equations;
// only Omega[1,4] compares a directly common construct.
//
// Sequence-level random intercepts (independent across equations and
// from the participant REs) absorb residual per-sequence idiosyncrasies
// not captured by the cues.
//
// Additive specification (no run x imbalance interaction).
// ============================================================

data {
  int<lower=1> P;                    // number of participants
  int<lower=1> S;                    // number of sequences

  // -------- extensive margin (bet vs. no-bet), N rows --------
  int<lower=1> N;
  array[N] int<lower=1, upper=P> pid;
  array[N] int<lower=1, upper=S> sid;
  vector<lower=0, upper=1>[N] imb_mag;   // unsigned imbalance magnitude
  vector<lower=0, upper=1>[N] run_dum;   // terminal-run presence dummy
  array[N] int<lower=0, upper=1> bet;

  // -------- side choice | bet == 1, N_bet rows --------
  int<lower=0> N_bet;
  array[N_bet] int<lower=1, upper=P> pid_s;
  array[N_bet] int<lower=1, upper=S> sid_s;
  vector[N_bet] run_h_s;
  vector[N_bet] imb_h_s;
  array[N_bet] int<lower=0, upper=1> side;   // 1 = H

  int<lower=0, upper=1> prior_only;          // 1 -> ignore likelihood (prior predictive)
}

parameters {
  vector[3] beta_bet;                // [intercept, run slope, imbalance slope]
  vector[3] beta_side;                // [intercept, run slope, imbalance slope]

  vector<lower=0>[6] tau;            // participant-level RE scales
  cholesky_factor_corr[6] L_Omega;   // participant-level RE correlation (Cholesky)
  matrix[6, P] z;                    // raw (non-centered) participant REs

  real<lower=0> sigma_v_bet;
  real<lower=0> sigma_v_side;
  vector[S] v_bet_raw;
  vector[S] v_side_raw;
}

transformed parameters {
  matrix[6, P] u = diag_pre_multiply(tau, L_Omega) * z;   // correlated participant REs
  vector[S] v_bet  = v_bet_raw  * sigma_v_bet;
  vector[S] v_side = v_side_raw * sigma_v_side;
}

model {
  // ---- priors ----
  beta_bet  ~ normal(0, 1.5);
  beta_side ~ normal(0, 1.5);
  tau       ~ student_t(3, 0, 1);
  L_Omega   ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();

  sigma_v_bet  ~ student_t(3, 0, 1);
  sigma_v_side ~ student_t(3, 0, 1);
  v_bet_raw    ~ std_normal();
  v_side_raw   ~ std_normal();

  // ---- likelihood ----
  if (!prior_only) {
    vector[N] eta_bet;
    for (n in 1:N) {
      eta_bet[n] = beta_bet[1] + u[1, pid[n]] + v_bet[sid[n]]
                 + (beta_bet[2] + u[2, pid[n]]) * run_dum[n]
                 + (beta_bet[3] + u[3, pid[n]]) * imb_mag[n];
    }
    bet ~ bernoulli_logit(eta_bet);

    if (N_bet > 0) {
      vector[N_bet] eta_side;
      for (m in 1:N_bet) {
        eta_side[m] = beta_side[1] + u[4, pid_s[m]] + v_side[sid_s[m]]
                    + (beta_side[2] + u[5, pid_s[m]]) * run_h_s[m]
                    + (beta_side[3] + u[6, pid_s[m]]) * imb_h_s[m];
      }
      side ~ bernoulli_logit(eta_side);
    }
  }
}

generated quantities {
  // log-lik and posterior-predictive replicates are not computed here;
  // PPC is done in R on a subsample of the extracted parameter draws
  // (see 97_synthesis_stan_outputs.R).
  corr_matrix[6] Omega = multiply_lower_tri_self_transpose(L_Omega);
  real corr_intercepts = Omega[1, 4];  // bet-propensity <-> side-bias (overall)
  real corr_run_slopes = Omega[2, 5];  // run-cue sensitivity: bet eq <-> side eq
  real corr_imb_slopes = Omega[3, 6];  // imbalance-cue sensitivity: bet eq <-> side eq
}