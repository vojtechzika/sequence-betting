# ---- Analysis Runner ----

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "analysis", "01_r_from_mpl.R"))
source(here::here("scripts", "analysis", "02_score_lotr.R"))
source(here::here("scripts", "analysis", "03_build_master_sequences.R"))

source(here::here("scripts", "analysis", "11_descriptive_participants.R"))
source(here::here("scripts", "analysis", "12_descriptive_sequences.R"))


dataset <- "pilot"   # change to "main" when ready


### Fitting MPL Data with Stan
fit_file  <- file.path(path_clean_ds(dataset), "mpl_fit.rds")
draw_file <- file.path(path_clean_ds(dataset), "mpl_r_draws.rds")
if (file.exists(fit_file) && file.exists(draw_file)) {
    msg("MPL Stan results already exist. Skipping estimation.")
  } else {
    msg("Running MPL Stan estimation...")
    r_from_mpl(dataset)
  }

#LOT-R scoring
score_lotr(dataset)

#Master sequence 
build_master_sequences(dataset)

descriptive_participants(dataset)
descriptive_sequences(dataset)

cat("\nAnalysis phase completed for:", dataset, "\n")