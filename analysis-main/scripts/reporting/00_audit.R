# ============================================================
# 00_audit.R
#
# PURPOSE
#   First step of the reporting stage. Verifies that the models listed in
#   selected-models.csv form a consistent set and that every model was
#   fitted on the samples as they are currently defined, then exports the
#   sample definitions for the replication package. The audit takes no
#   action on convergence or posterior predictive checks (the choice of
#   models is made in selected-models.csv); it stops only when the set of
#   models is inconsistent.
#
#   1. Lineage: every upstream folder recorded in the inputs.csv of a
#      selected model must itself be selected (check_lineage()).
#   2. Samples: for each selected model, the number of participants and
#      trials recorded in diagnostics.csv for each sample is compared with
#      the current sample-definitions.csv of the MPL model. A mismatch
#      means that the model was fitted on a sample that has since changed.
#   3. Diagnostics: all checks of the selected models are collected in one
#      table with the failed checks listed in the report.
#
# INPUT
#   selected-models.csv
#   data/models/<folder>/diagnostics.csv       (every selected model)
#   model "mpl": sample/sample-definitions.csv, sample/analysis-trials.csv
#
# OUTPUT
#   data/output/csv/audit_diagnostics.csv        diagnostics of all selected models
#   data/output/csv/audit_samples.csv            recorded vs defined sample sizes
#   data/output/csv/audit_sample-definitions.csv definition of every sample
#   data/output/csv/audit_analysis-trials.csv    trial table with sample membership
#   data/output/reports/audit.txt                summary of the audit
#
# NOTES
#   - Files are written before the audit stops, so that a failure can be
#     inspected
#   - Sample sizes are read from the check names <sample>_participants and
#     <sample>_trials written by the model scripts
# ============================================================

audit <- function(cfg) {
  
  sel <- read_selected_models()
  if (is.null(sel)) stop("selected-models.csv not found.", call. = FALSE)
  
  # ---- Lineage ----
  lineage_msg <- ""
  lineage_ok <- tryCatch({ check_lineage(); TRUE },
                         error = function(e) { lineage_msg <<- conditionMessage(e); FALSE })
  
  # ---- Diagnostics of the selected models ----
  diag_all <- rbindlist(lapply(seq_len(nrow(sel)), function(k) {
    f <- file.path(path_mod, sel$folder[k], "diagnostics.csv")
    if (!file.exists(f)) stop("diagnostics.csv missing for '", sel$key[k], "' (", f, ")",
                              call. = FALSE)
    d <- fread(f, encoding = "UTF-8")
    d[, `:=`(key = sel$key[k], folder = sel$folder[k])]
    setcolorder(d, c("key", "folder"))
    d
  }), fill = TRUE)
  
  # ---- Recorded vs defined sample sizes ----
  defs <- fread(model_path("mpl", file.path("sample", "sample-definitions.csv")),
                encoding = "UTF-8")
  defs[, `:=`(n_participants = as.integer(n_participants), n_trials = as.integer(n_trials))]
  
  rec <- diag_all[grepl("_(participants|trials)$", check)]
  rec[, `:=`(sample_name = sub("_(participants|trials)$", "", check),
             what        = sub("^.*_", "", check))]
  rec_wide <- dcast(rec, key + folder + treatment + sample_name + sample ~ what,
                    value.var = "value")
  setnames(rec_wide, c("participants", "trials"), c("n_participants_model", "n_trials_model"))
  
  samples <- merge(rec_wide,
                   defs[, .(treatment, sample_name = sample,
                            n_participants_definition = n_participants,
                            n_trials_definition       = n_trials)],
                   by = c("treatment", "sample_name"), all.x = TRUE)
  samples[, agrees := !is.na(n_participants_definition) &
            n_participants_model == n_participants_definition &
            n_trials_model       == n_trials_definition]
  setcolorder(samples, c("key", "folder", "treatment", "sample_name", "sample"))
  setorder(samples, key, treatment, sample_name, sample)
  
  # ---- Exports ----
  fwrite(diag_all, file.path(path_csv, "audit_diagnostics.csv"))
  fwrite(samples,  file.path(path_csv, "audit_samples.csv"))
  fwrite(defs,     file.path(path_csv, "audit_sample-definitions.csv"))
  fwrite(fread(model_path("mpl", file.path("sample", "analysis-trials.csv")),
               encoding = "UTF-8"),
         file.path(path_csv, "audit_analysis-trials.csv"))
  
  # ---- Report ----
  failed   <- diag_all[pass == FALSE]
  mismatch <- samples[agrees == FALSE]
  
  lines <- c(
    "AUDIT OF SELECTED MODELS",
    paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Selected models",
    paste0("  ", sel$key, " -> ", sel$folder),
    "",
    paste("Lineage:", if (lineage_ok) "consistent" else paste("FAILED\n ", lineage_msg)),
    "",
    paste0("Sample sizes: ", sum(samples$agrees), " of ", nrow(samples),
           " recorded samples match the current definitions"),
    if (nrow(mismatch)) c("  Mismatches:",
                          paste0("  ", mismatch$key, " ", mismatch$treatment, " ",
                                 mismatch$sample, ": model ", mismatch$n_participants_model,
                                 "/", mismatch$n_trials_model, ", definition ",
                                 mismatch$n_participants_definition, "/",
                                 mismatch$n_trials_definition, " (participants/trials)")),
    "",
    paste0("Diagnostics: ", nrow(failed), " of ", nrow(diag_all), " checks not passed"),
    if (nrow(failed)) paste0("  ", failed$key, " ", failed$tag, " ", failed$check,
                             ": value ", signif(failed$value, 4),
                             ", threshold ", signif(failed$threshold, 4))
  )
  f_rep <- file.path(path_rep, "audit.txt")
  writeLines(lines, f_rep)
  msg("Saved: ", f_rep)
  
  if (!lineage_ok || nrow(mismatch)) {
    stop("Audit failed; see ", f_rep, call. = FALSE)
  }
  invisible(TRUE)
}
