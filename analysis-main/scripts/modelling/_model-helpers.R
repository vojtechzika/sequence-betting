# ============================================================
# _model-helpers.R
#
# PURPOSE
#   Shared helpers for the modelling stage. Each model script writes
#   all its artefacts into one folder under path_mod. The folder is
#   built under a temporary name (<folder>.tmp) and renamed to
#   <folder> only when the run is complete, i.e. when
#   diagnostics.csv has been written. A complete folder is never
#   overwritten; to rerun a model, delete its folder.
#
#   Scripts locate the upstream models they read through
#   selected-models.csv (project root), which maps a model key to the
#   folder chosen after inspecting the diagnostics. Modelling and reporting
#   therefore share one source for the version of every model.
#
# FILES READ
#   selected-models.csv (project root)   -- columns: key, folder
#                                           (optional column: note)
#
# FILES PRODUCED (inside each model folder)
#   diagnostics.csv   -- long format: tag, treatment, sample, check,
#                        value, threshold, pass (written by the model)
#   inputs.csv        -- upstream folders the run resolved via
#                        model_path() (written by model_end)
#   run-info.txt      -- date, seed, R and package versions, model
#                        settings (written by model_end)
#
# FUNCTIONS
#   model_begin(folder)      start a run; returns the .tmp directory,
#                            or NULL if the model is already complete
#   model_path(x, file)      path of a file in a selected model folder
#   write_diagnostics(run, tbl)
#   stan_diagnostics(fit, tag, ...)   shared convergence checks
#   model_end(run, cfg, settings, key)  finish a run (rename the folder)
#   a_star_pid_set(tr, sample, tau, set)
#                            participants kept/excluded/flagged by the
#                            a* classification (07_a-star-draws.R)
#   hl_sample(cfg)           Holt-Laury sample of the confirmatory analyses
#                            ("consistent" or "full", from consistent_only)
#   check_lineage()          stop if a selected model was built from an
#                            upstream folder that is no longer selected
#
# TYPICAL USE IN A MODEL SCRIPT
#   run <- model_begin("rq1-bernoulli")
#   if (is.null(run)) return(invisible(NULL))
#   draws <- readRDS(model_path("mpl", "r-draws-m25.rds"))
#   ... fit, then saveRDS(fit, file.path(run, "fit-m25.rds")) ...
#   write_diagnostics(run, rbind(stan_diagnostics(fit, "m25", cfg = cfg), ...))
#   model_end(run, cfg, settings = cfg$model$stan$rq1_bernoulli)
# ============================================================

# Upstream folders resolved during the current run
.model_log <- new.env(parent = emptyenv())
.model_log$used <- NULL

# --------------------------------------------------------
# Selected models
# --------------------------------------------------------
read_selected_models <- function() {
  f <- here::here("selected-models.csv")
  if (!file.exists(f)) return(NULL)
  # lines starting with '#' are comments
  lines <- readLines(f, encoding = "UTF-8", warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  sel   <- data.table::fread(text = lines, encoding = "UTF-8", colClasses = "character")
  if (!all(c("key", "folder") %in% names(sel))) {
    stop("selected-models.csv must have columns 'key' and 'folder'.")
  }
  if (anyDuplicated(sel$key)) stop("selected-models.csv has duplicate keys.")
  sel
}

# Path of `file` in the model folder `x`. `x` is a key in
# selected-models.csv or, if it is not a key, a folder name under
# path_mod. With must_exist = FALSE the folder and file need not exist
# (used to build output paths).
model_path <- function(x, file = NULL, must_exist = TRUE) {
  stopifnot(length(x) == 1L, nzchar(x))
  
  sel    <- read_selected_models()
  is_key <- !is.null(sel) && x %in% sel$key
  folder <- if (is_key) sel$folder[sel$key == x] else x
  dir    <- file.path(path_mod, folder)
  
  if (must_exist) {
    if (!dir.exists(dir)) {
      stop("Model folder not found for '", x, "': ", dir,
           if (!is_key) "\n  (not a key in selected-models.csv)" else "")
    }
    .model_log$used <- rbind(.model_log$used, data.frame(key = x, folder = folder))
  }
  
  if (is.null(file)) return(dir)
  p <- file.path(dir, file)
  if (must_exist && !file.exists(p)) stop("File not found in model folder: ", p)
  p
}

# --------------------------------------------------------
# Run lifecycle
# --------------------------------------------------------
# Returns the temporary directory to write into, or NULL when the
# model folder is already complete.
model_begin <- function(folder) {
  stopifnot(length(folder) == 1L, nzchar(folder), !grepl("\\.tmp$", folder))
  
  final <- file.path(path_mod, folder)
  tmp   <- paste0(final, ".tmp")
  
  if (dir.exists(final)) {
    msg("Model '", folder, "' is complete; delete ", final, " to rerun.")
    return(NULL)
  }
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE)   # stale partial run
  dir.create(tmp, recursive = TRUE)
  
  .model_log$used <- NULL
  tmp
}

# Writes diagnostics.csv in the model folder under construction.
write_diagnostics <- function(run, tbl) {
  tbl  <- data.table::as.data.table(tbl)
  cols <- c("tag", "treatment", "sample", "check", "value", "threshold", "pass")
  if (!identical(names(tbl), cols)) {
    stop("diagnostics must have columns: ", paste(cols, collapse = ", "))
  }
  if (nrow(tbl) == 0L || anyNA(tbl$pass)) {
    stop("diagnostics must be non-empty and 'pass' must not be missing.")
  }
  tbl[, pass := as.logical(pass)]
  data.table::fwrite(tbl, file.path(run, "diagnostics.csv"))
  invisible(tbl)
}

# Convergence checks for one Stan fit, in the long diagnostics format.
# Thresholds: cfg$model$convergence (rhat_max, ess_min, divergences);
# treedepth hits must be zero and E-BFMI at least 0.3.
stan_diagnostics <- function(fit, tag, treatment = NA_character_,
                             sample = NA_character_, cfg,
                             max_treedepth = 10L, ebfmi_min = 0.3) {
  conv <- cfg$model$convergence
  summ <- rstan::summary(fit)$summary
  sp   <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  # Quantities without posterior variation (for example the diagonal of a
  # correlation matrix, equal to 1 in every draw) have no defined ESS or
  # R-hat and are left out of both minima; their number is reported.
  varying  <- is.finite(summ[, "sd"]) & summ[, "sd"] > 1e-10
  rhat_max <- max(summ[varying, "Rhat"],  na.rm = TRUE)
  ess_min  <- min(summ[varying, "n_eff"], na.rm = TRUE)
  diverg   <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  tdhits   <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth),
                         numeric(1)))
  # E-BFMI per chain: squared successive energy differences relative to
  # the energy variance; the smallest chain value is reported
  ebfmi <- min(vapply(sp, function(x) {
    e <- x[, "energy__"]
    sum(diff(e)^2) / sum((e - mean(e))^2)
  }, numeric(1)))
  
  data.table::data.table(
    tag       = tag,
    treatment = as.character(treatment),
    sample    = as.character(sample),
    check     = c("rhat_max", "ess_min", "divergences", "treedepth_hits", "ebfmi_min",
                  "constant_quantities_excluded"),
    value     = c(rhat_max, ess_min, diverg, tdhits, ebfmi, sum(!varying)),
    threshold = c(conv$rhat_max, conv$ess_min, conv$divergences, 0, ebfmi_min, NA),
    pass      = c(rhat_max <= conv$rhat_max, ess_min >= conv$ess_min,
                  diverg <= conv$divergences, tdhits == 0, ebfmi >= ebfmi_min, TRUE)
  )
}

# Completes a run: requires diagnostics.csv, records the inputs and
# the run settings, and renames <folder>.tmp to <folder>. With `key`, it
# prints the line to add to selected-models.csv.
model_end <- function(run, cfg, settings = NULL, key = NULL) {
  if (!file.exists(file.path(run, "diagnostics.csv"))) {
    stop("diagnostics.csv missing; the folder stays as ", run)
  }
  
  if (!is.null(.model_log$used)) {
    data.table::fwrite(unique(.model_log$used), file.path(run, "inputs.csv"))
  }
  
  pk <- c("rstan", "StanHeaders", "data.table")
  info <- c(
    paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("seed:", cfg$run$seed),
    paste("R:", R.version.string),
    paste0(pk, ": ", vapply(pk, function(p) as.character(utils::packageVersion(p)),
                            character(1))),
    if (!is.null(settings)) c("settings:", utils::capture.output(dput(settings)))
  )
  writeLines(info, file.path(run, "run-info.txt"))
  
  final <- sub("\\.tmp$", "", run)
  if (!file.rename(run, final)) stop("Could not rename ", run, " to ", final)
  
  msg("Saved: ", final)
  if (!is.null(key)) {
    msg("To use this model downstream, set it in selected-models.csv: ",
        key, ",", basename(final))
  }
  invisible(final)
}

# Holt-Laury sample on which the confirmatory analyses are based:
# "consistent" if cfg$run$consistent_only is TRUE, otherwise "full".
# Used wherever a script looks up the a* flags of the confirmatory sample.
hl_sample <- function(cfg) {
  if (isTRUE(cfg$run$consistent_only)) "consistent" else "full"
}

# --------------------------------------------------------
# Lineage
# --------------------------------------------------------
# Checks that the models in selected-models.csv are consistent: every
# upstream folder recorded in the inputs.csv of a selected model must itself
# be a selected folder. A model built from an MPL fit that is no longer
# selected is therefore reported as stale. Stops with a list of the stale
# models; returns TRUE invisibly otherwise.
check_lineage <- function() {
  sel <- read_selected_models()
  if (is.null(sel)) stop("selected-models.csv not found; nothing to check.")
  
  problems <- character()
  for (k in seq_len(nrow(sel))) {
    dir <- file.path(path_mod, sel$folder[k])
    if (!dir.exists(dir)) {
      problems <- c(problems, sprintf("'%s': folder %s does not exist.",
                                      sel$key[k], sel$folder[k]))
      next
    }
    f <- file.path(dir, "inputs.csv")
    if (!file.exists(f)) next                       # model without upstream
    
    inp   <- data.table::fread(f, encoding = "UTF-8", colClasses = "character")
    stale <- setdiff(unique(inp$folder), sel$folder)
    if (length(stale)) {
      problems <- c(problems, sprintf(
        "'%s' (%s) was built from %s, which is not selected.",
        sel$key[k], sel$folder[k], paste(stale, collapse = ", ")))
    }
  }
  
  if (length(problems)) {
    stop("Lineage check failed:\n  ", paste(problems, collapse = "\n  "),
         "\nRerun the stale models or correct selected-models.csv.", call. = FALSE)
  }
  invisible(TRUE)
}

# --------------------------------------------------------
# a* classification
# --------------------------------------------------------
# Participant IDs from the a* classification written by a_star_draws()
# in <folder of model "mpl">/a-star/. `set` is one of
#   "keep"    / "exclude"  for treatments with betting_normative = TRUE
#                          (exclude = normative non-betters, P0 >= tau)
#   "flag"    / "unflag"   for treatments with betting_normative = FALSE
#                          (flag = normative betters, Pplus >= tau)
# `tau` must be one of cfg$design$a_flags$tau.
a_star_pid_set <- function(tr, sample, tau,
                           set = c("keep", "exclude", "flag", "unflag")) {
  set <- match.arg(set)
  f   <- model_path("mpl",
                    file.path("a-star", paste0("a-star-flags-", tr, "-", sample, ".csv")))
  tbl <- data.table::fread(f, encoding = "UTF-8", colClasses = c(pid = "character"))
  
  tau_nm <- gsub("\\.", "", sprintf("%.2f", as.numeric(tau)))
  col    <- if (set %in% c("keep", "exclude")) paste0("exclude_nonbetters_tau", tau_nm)
            else paste0("flag_betters_tau", tau_nm)
  if (!col %in% names(tbl)) {
    stop("Column '", col, "' not found in ", f,
         ". Check that '", set, "' fits the treatment and that tau was ",
         "among the thresholds used.", call. = FALSE)
  }
  
  hit <- if (set %in% c("keep", "unflag")) !tbl[[col]] else tbl[[col]]
  sort(unique(tbl$pid[hit]))
}
