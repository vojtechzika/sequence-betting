# ============================================================
# 02_model_cfg.R
# Estimation, PPC, and simulation controls
# ============================================================
model_cfg <- function() {
  
  list(
    
    # --------------------------------------------------------
    # Stan sampling settings
    # --------------------------------------------------------
    stan = list(
      
      mpl = list(
        prior = list(# Bland (2023)
          r_mean      = 0.27,
          r_sd        = 0.36,   
          lambda_mean = 30,     # lognormal mean on natural scale
          lambda_sd   = 0.5     # lognormal SD on log scale
        ),
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq1 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq2 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq3 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      rq4 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      ex1_2 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      ),
      
      ex2 = list(
        iter = 2000, warmup = 1000, chains = 4, adapt_delta = 0.95, treedepth = 12
      )
    ),
    
    # --------------------------------------------------------
    # Posterior predictive replicates
    # --------------------------------------------------------
    ppc = list(
      rq1_k        = 1000L,
      rq3_k        = 1000L,
      rq4_k        = 1000L,
      rq1_interval = c(0.05, 0.95)
    ),
    
    # --------------------------------------------------------
    # Simulation / resampling reps
    # --------------------------------------------------------
    simulation = list(
      ex1_trep = 1000L,
      ex2_trep = 1000L
    ),
    
    # --------------------------------------------------------
    # RQ3 model choice
    # --------------------------------------------------------
    rq3 = list(
      model_type = "hurdle_gamma"  # or "gaussian"
    )
  )
}