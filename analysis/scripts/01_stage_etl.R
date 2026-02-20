############################
# ---- ETL Stage ----
############################

run_etl <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$run$dataset))
  ds <- cfg$run$dataset
  
  source(here::here("scripts", "etl", "01_read_merge_filter.R"))
  source(here::here("scripts", "etl", "02_table_participants.R"))
  source(here::here("scripts", "etl", "03_table_sequences.R"))
  source(here::here("scripts", "etl", "04_table_mpl.R"))
  source(here::here("scripts", "etl", "05_table_lotr.R"))
  source(here::here("scripts", "etl", "06_table_debriefing.R"))
  
  run_read_merge_filter(ds)
  make_participants_table(ds)
  make_sequences_table(ds)
  make_mpl_table(ds)
  make_lotr_table(ds)
  make_debriefing_table(ds)
  
  msg("\nETL phase completed for:", ds, "\n")
  
  invisible(TRUE)
}