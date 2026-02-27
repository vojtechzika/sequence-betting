run_cfg <- function() {
  list(
    # define dataset - "pilot" for testing/pilot data or "main" for data from a full experiment
    dataset   = "pilot",
    
    # A list of all treatments to run the analysis on. List all treatments unless necessary to list only some.
    # Confirmatory analysis is run on the first listed treatment!!
    treatment = c("m25","m19"),
    
    # Should models and outputs be re-generated if they already exist? Set FALSE for reproducibility, TRUE for testing.
    overwrite_outputs = TRUE, 
    overwrite_models = TRUE,
    
    # Seed for reproducibility.
    seed      = 12345,
    
    rq1_ppc_k = 300L,
    rq3_model = "hurdle_gamma", # accepted: "hurdle_gamma", "gaussian"
    rq3_krep = 300L, # set to 1000L for main data
    ex1_trep = 300L,  # set to 1000L for main data
    ex2_trep = 300L,  # set to 1000L for main data
    ex2_min_Nk = 3L, # minimum number of participants
    ex2_rq2_min_bets = 3L # minimum number of bets
  )
}