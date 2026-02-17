# ---- ETL Runner ----
source("scripts/01_read_merge_filter.R")
source("scripts/02_table_participants.R")
source("scripts/03_table_sequences.R")

dataset <- "pilot"   # change to "main" when ready

run_read_merge_filter(dataset)
make_participants_table(dataset)
make_sequences_table(dataset)

cat("\nETL phase completed for:", dataset, "\n")


