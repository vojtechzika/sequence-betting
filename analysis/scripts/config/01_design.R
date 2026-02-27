design_cfg <- function() {
  list(
    seq = list(
      endowment  = 100L,
      xmin       = 0.01,
      seq_n      = 64L,
      coin_prob  = 0.5,      # objective probability of win
      treatments = list(
        m25 = 2.5, m19 = 1.9
        ),
      side_labels = list(
        heads = "H",
        tails = "O",
        nobet = "NB"
      ),
      anchor_labels = list(
        pure_heads = "HHHHHH",
        pure_tails = "OOOOOO"
      )
    ),
    mpl = list(
      K = 10L,
      A_high = 20.0, A_low = 16.0,
      B_high = 38.5, B_low = 1.0
    ),
    
    
    a_flags = list(
      # For each treatment label: TRUE  = a* > 0 is normatively optimal (betting benchmark); FALSE = a* = 0 is normatively optimal (no-bet benchmark)
      betting_normative = list(
        m25 = TRUE,
        m19 = FALSE
      ),
      tau = c(0.90, 0.80, 0.95) # Thresholds for classification; first = main, rest = sensitivity
    ),
    
    rhos = list(
      rq1_rho = c(0.10, 0.08, 0.12),
      rq2_rho = c(0.05, 0.03, 0.08),
      rq2_sd_floor = 2,
      rq3_rho = c(0.05, 0.03, 0.08),
      ex1_eps = c(0.05, 0.03, 0.08),
      ex1_delta = c(0.05, 0.03, 0.08) 
    ),
    exclusion = list(
      rq2_min_bets = 3L,
      P0 = c(0.80, 0.90, 0.95)
    )
  )
}