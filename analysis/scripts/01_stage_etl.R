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
  source(here::here("scripts", "etl", "09_build_master_sequences.R"))
  
  run_read_merge_filter(cfg)
  make_participants_table(cfg)
  make_sequences_table(cfg)
  make_mpl_table(cfg)
  make_lotr_table(cfg)
  make_debriefing_table(cfg)
  build_master_sequences(cfg)
  
  msg("\nETL completed for data_folder:", cfg$run$data_folder, "\n")
  
  invisible(TRUE)
}