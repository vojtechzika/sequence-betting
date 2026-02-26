run_cfg <- function() {
  list(
    dataset   = "pilot",
    treatment = c("m25","m19"),  # A list of all treatments to run the analysis on. Add "pooled = TRUE" when calling the function from runners to include a pooled analysis across all treatments.
    seed      = 12345,
    rq3_model = "hurdle_gamma", # accepted: "hurdle_gamma", "gaussian"
    rq3_krep = 300L, # set to 1000L for main data
    ex1_trep = 300L,  # set to 1000L for main data
    ex2_trep = 300L,  # set to 1000L for main data
    ex2_min_Nk = 3L, # minimum number of participants
    ex2_rq2_min_bets = 3L # minimum number of bets
    )
}