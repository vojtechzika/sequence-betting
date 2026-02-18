# ---- Analysis Runner ----

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "analysis", "01_r_from_mpl.R"))

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

cat("\nAnalysis phase completed for:", dataset, "\n")