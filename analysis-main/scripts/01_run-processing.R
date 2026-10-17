############################
# ---- Processing Stage ----
############################

run_processing <- function(cfg) {
  
  source(here::here("scripts", "processing", "01_read-merge-filter.R"))
  source(here::here("scripts", "processing", "11_table-participants.R"))
  source(here::here("scripts", "processing", "12_table-sequences.R"))
  source(here::here("scripts", "processing", "13_table-holt-laury.R"))
  source(here::here("scripts", "processing", "14_table-lotr.R"))
  source(here::here("scripts", "processing", "15_table-debriefing.R"))
  source(here::here("scripts", "processing", "19_master-sequences.R"))
  source(here::here("scripts", "processing", "21_score-lotr.R"))
  source(here::here("scripts", "processing", "22_score-switching-point-r.R"))
  source(here::here("scripts", "processing", "31_button-order.R"))
  
  read_merge_filter(cfg)
  
  table_participants(cfg)
  table_sequences(cfg)
  table_holt_laury(cfg)
  table_lotr(cfg)
  table_debriefing(cfg)
  master_sequences(cfg)
  
  score_lotr(cfg)
  score_switching_point_r(cfg)
  
  button_order_check(cfg)
  
  msg("\nProcessing completed.\n")
  
  invisible(TRUE)
}