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
        pilot = list(iter=1000, warmup=750,  chains=2, adapt_delta=0.90, treedepth=12),
        main  = list(iter=2000, warmup=1000, chains=4, adapt_delta=0.95, treedepth=12)
      ),
      
      rq1 = list(
        pilot = list(iter=1500, warmup=750,  chains=2, adapt_delta=0.90, treedepth=12),
        main  = list(iter=2000, warmup=1000, chains=4, adapt_delta=0.95, treedepth=12)
      ),
      
      rq2 = list(
        pilot = list(iter=1500, warmup=750,  chains=2, adapt_delta=0.90, treedepth=12),
        main  = list(iter=2000, warmup=1000, chains=4, adapt_delta=0.95, treedepth=12)
      ),
      
      rq3 = list(
        pilot = list(iter=1500, warmup=750,  chains=2, adapt_delta=0.90, treedepth=12),
        main  = list(iter=2000, warmup=1000, chains=4, adapt_delta=0.95, treedepth=12)
      )
    ),
    
    # --------------------------------------------------------
    # Posterior predictive replicates
    # --------------------------------------------------------
    ppc = list(
      rq1_k = list(pilot = 300L,  main = 1000L),
      rq3_k = list(pilot = 300L,  main = 1000L),
      
      rq1_interval = c(0.05, 0.95)
    ),
    
    # --------------------------------------------------------
    # Simulation / resampling reps
    # --------------------------------------------------------
    simulation = list(
      ex1_trep = list(pilot = 300L,  main = 1000L),
      ex2_trep = list(pilot = 300L,  main = 1000L)
    ),
    
    # --------------------------------------------------------
    # RQ3 model choice
    # --------------------------------------------------------
    rq3 = list(
      model_type = "hurdle_gamma"  # or "gaussian"
    )
  )
}