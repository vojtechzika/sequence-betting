# ---- ETL Runner ----
source(here::here("scripts", "00_setup.R"))

source(here::here("scripts", "etl", "01_read_merge_filter.R"))
source(here::here("scripts", "etl", "02_table_participants.R"))
source(here::here("scripts", "etl", "03_table_sequences.R"))
source(here::here("scripts", "etl", "04_table_mpl.R"))
source(here::here("scripts", "etl", "05_table_lotr.R"))
source(here::here("scripts", "etl", "06_table_debriefing.R"))


dataset <- "pilot"   # change to "main" when ready

run_read_merge_filter(dataset)
make_participants_table(dataset)
make_sequences_table(dataset)
make_mpl_table(dataset)
make_lotr_table(dataset) 
make_debriefing_table(dataset)

cat("\nETL phase completed for:", dataset, "\n")


