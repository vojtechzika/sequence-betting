run_cfg <- function() {
  list(
    dataset   = "pilot",
    treatment = c("m25","m19"),  # A list of all treatments to run the analysis on. Add "pooled = TRUE" when calling the function from runners to include a pooled analysis across all treatments.
    seed      = 12345,
    rq3_model = "hurdle_gamma", # accepted: "hurdle_gamma", "gaussian"
    rq3_krep = 500L, # set to 99999L if you want to run the full model, but beware that it will take a very long time to run. For testing purposes, set to a smaller number (e.g., 200L) to get a quick estimate of the model parameters and check that everything is working correctly.
    ex1_trep = 500L  # set to 99999L if you want to run the full model, but beware that it will take a very long time to run. For testing purposes, set to a smaller number (e.g., 200L) to get a quick estimate of the model parameters and check that everything is working correctly.

    )
}