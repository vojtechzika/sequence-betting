run_cfg <- function() {
  list(
    dataset   = "pilot",
    treatment = c("m25","m19"),  # A list of all treatments to run the analysis on. Add "pooled = TRUE" when calling the function from runners to include a pooled analysis across all treatments.
    seed      = 12345,
    rq3_model = "hurdle_gamma", # accepted: "hurdle_gamma", "gaussian"
    rq3_krep = 200L
  )
}