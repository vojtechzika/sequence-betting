# ============================================================
# 01_design_cfg.R
# Experimental design + inferential thresholds
# ============================================================

design_cfg <- function() {
  
  
  
  list(
    
    # --------------------------------------------------------
    # Core sequence design
    # --------------------------------------------------------
    
    seq = list(
      endowment = 100L,
      xmin      = 0.01,
      seq_n     = 64L,
      coin_prob = 0.5,
      
      treatments = list(
        m25 = 2.5,
        m19 = 1.9
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
   
    
    
    
    # --------------------------------------------------------
    # MPL risk elicitation
    # --------------------------------------------------------
    mpl = list(
      K      = 19L, # number of rows/lotteries
      A_high = 20.0,
      A_low  = 16.0,
      B_high = 38.5,
      B_low  = 1.0,
      label_safe  = "A",
      label_risky = "B"
    ),
    
    # --------------------------------------------------------
    # LOT-R (scoring)
    # --------------------------------------------------------
    lotr = list(
      scale_min = 0L,
      scale_max = 4L,
      item_cols = paste0("q", 1:10),
      scored_items = c("q1", "q3", "q4", "q7", "q9", "q10"),
      rev_items    = c("q3", "q7", "q9")
    ),
    
    # --------------------------------------------------------
    # RQ1 (Extensive margin)
    # --------------------------------------------------------
    rq1 = list(
      rho = c(0.10, 0.08, 0.12)  # first = main
    ),
    
    # --------------------------------------------------------
    # RQ2 (Intensive margin)
    # --------------------------------------------------------
    rq2 = list(
      rho      = c(0.05, 0.03, 0.08),
      sd_floor = 2,
      min_bets = 3L
    ),
    
    # --------------------------------------------------------
    # RQ3
    # --------------------------------------------------------
    rq3 = list(
      rho = c(0.05, 0.03, 0.08)
    ),
    
    # --------------------------------------------------------
    # RQ4 (Side choice on betting trials)
    # --------------------------------------------------------
    rq4 = list(
      delta = c(0.05, 0.03, 0.08),   # first = main (tolerance band around hbar)
      ppc_overdisp_cut = 0.95        # prereg trigger cutoff for Beta–Binomial robustness
    ),
    
    
    # --------------------------------------------------------
    # EX2
    # --------------------------------------------------------
    ex2 = list(
      min_Nk      = 10L
    ),
    
    # --------------------------------------------------------
    # Normative a* classification
    # --------------------------------------------------------
    a_flags = list(
      betting_normative = list(
        m25 = TRUE,
        m19 = FALSE
      ),
      tau = c(0.90, 0.80, 0.95) # first value is for the primary analysis
    ),
    
    # --------------------------------------------------------
    # Additional exclusion thresholds
    # --------------------------------------------------------
    exclusion = list(
      P0 = c(0.80, 0.90, 0.95)
    )
  )
}