############################
# ---- ETL Stage ----
############################

run_etl <- function(cfg) {
  
  source(here::here("scripts", "etl", "01_read_merge_filter.R"))
  source(here::here("scripts", "etl", "02_table_participants.R"))
  source(here::here("scripts", "etl", "03_table_sequences.R"))
  source(here::here("scripts", "etl", "04_table_mpl.R"))
  source(here::here("scripts", "etl", "05_table_lotr.R"))
  source(here::here("scripts", "etl", "06_table_debriefing.R"))
  
  run_read_merge_filter(cfg)
  make_participants_table(cfg)
  make_sequences_table(cfg)
  make_mpl_table(cfg)
  make_lotr_table(cfg)
  make_debriefing_table(cfg)
  
  msg("\nETL completed for:", ds, "\n")
  
  invisible(TRUE)
}