# ============================================================
# 07_a-star-draws.R
#
# PURPOSE
#   Computes the expected-utility optimal stake a*(r) for every
#   posterior draw of the risk parameter r_i from the Holt-Laury model,
#   separately for each treatment (payout multiplier m). a* is the
#   integer maximizer of EU(a; r) over {0, ..., endowment}; ties are
#   broken towards the smaller stake. From the draws, the script
#   classifies participants by the posterior probability that the
#   optimal stake is zero (P0) or positive (Pplus):
#     betting_normative = TRUE   participants with P0 >= tau are
#                                normative non-betters and are excluded
#                                from the confirmatory sample
#     betting_normative = FALSE  participants with Pplus >= tau are
#                                flagged as normative betters
#   The first tau in cfg$design$a_flags$tau is the primary threshold;
#   the others are sensitivity thresholds.
#
#   The script works only on the MPL model selected under the key "mpl"
#   in selected-models.csv and writes into a subfolder of that model, so
#   a* is computed only for the MPL model in use.
#
# SAMPLES
#   full        all Holt-Laury participants of the treatment
#   consistent  Holt-Laury-consistent participants only
#               (when cfg$run$consistent_only is TRUE); the draws are those of
#               the full sample restricted to these participants, so each
#               participant has one set of a* draws in every sample
#
# INPUT
#   model "mpl": r-draws-<tr>.rds, scored-<tr>.csv (inconsistency flag),
#                model-info.csv (partial_pooling must be FALSE)
#   scripts/config/02_eu-functions.R   (a_star_discrete, R_TOL_DEFAULT)
#   cfg$design$seq (endowment, xmin, coin_prob, treatments),
#   cfg$design$a_flags (betting_normative, tau)
#
# OUTPUT (subfolder a-star/ of the selected MPL model folder)
#   a-star-draws-<tr>-<sample>.rds   list(pid, a_star_draws), draws x pid
#   a-star-flags-<tr>-<sample>.csv   pid, P0, Pplus, betting_normative and
#                                    one flag column per tau
#   diagnostics.csv                  consistency checks per treatment/sample
#   a-star-flag-summary.csv          participants and flags per threshold
#   inputs.csv, run-info.txt
#   Flag columns: exclude_nonbetters_tau<tau> (P0 >= tau) when
#   betting_normative is TRUE, flag_betters_tau<tau> (Pplus >= tau)
#   otherwise; <tau> is written without the decimal point (0.90 -> 090).
#   a_star_pid_set() (_model-helpers.R) reads these files.
#
# NOTES
#   - Computed separately per treatment (no pooling)
#   - The subfolder is built as a-star.tmp and renamed when complete;
#     an existing a-star folder is never overwritten (delete it to rerun,
#     for example after changing the design parameters or tau)
#   - A missing MPL fit stops the script
# ============================================================

a_star_draws <- function(cfg) {
  
  # a-star subfolder inside the selected MPL model folder
  mpl_dir <- model_path("mpl")
  
  # The consistent sample is obtained by filtering the full-sample draws, which
  # is valid only if r_i is not pooled across participants (model property
  # recorded by the MPL script)
  info <- data.table::fread(model_path("mpl", "model-info.csv"), encoding = "UTF-8",
                            colClasses = "character")
  pooling <- info[property == "partial_pooling", value]
  if (!identical(pooling, "FALSE")) {
    stop("The selected MPL model is not recorded as non-pooled (partial_pooling = ",
         paste(pooling, collapse = ""), "). Deriving the consistent sample by ",
         "filtering is valid only without pooling; a pooled model needs its own ",
         "consistent-sample path in this script.", call. = FALSE)
  }
  run     <- model_begin(file.path(basename(mpl_dir), "a-star"))
  if (is.null(run)) return(invisible(NULL))
  
  design          <- cfg$design
  consistent_only <- isTRUE(cfg$run$consistent_only)
  treatments      <- names(design$seq$treatments)
  samples         <- c("full", if (consistent_only) "consistent")
  
  e     <- design$seq$endowment
  xmin  <- design$seq$xmin
  p_win <- design$seq$coin_prob
  grid  <- 0:e
  r_tol <- get("R_TOL_DEFAULT", inherits = TRUE)
  
  # ---- Classification thresholds ----
  bn      <- design$a_flags$betting_normative
  tau_vec <- as.numeric(design$a_flags$tau)
  stopifnot(!is.null(bn), all(treatments %in% names(bn)),
            length(tau_vec) >= 1L, all(is.finite(tau_vec)),
            all(tau_vec > 0), all(tau_vec < 1))
  
  diag_rows   <- list()
  count_rows  <- list()
  
  for (tr in treatments) {
    m <- design$seq$treatments[[tr]]
    
    # r draws of the full-sample fit of this treatment
    hl <- readRDS(model_path("mpl", paste0("r-draws-", tr, ".rds")))
    stopifnot(is.list(hl), all(c("pid", "r_draws") %in% names(hl)))
    pid_full <- as.character(hl$pid)
    r_draws  <- hl$r_draws
    stopifnot(is.matrix(r_draws), ncol(r_draws) == length(pid_full))
    
    # a* for every draw and participant (discrete maximization over the grid)
    a_star_full <- matrix(NA_integer_, nrow = nrow(r_draws), ncol = ncol(r_draws))
    for (i in seq_len(ncol(r_draws))) {
      a_star_full[, i] <- vapply(
        r_draws[, i], FUN = a_star_discrete, FUN.VALUE = integer(1),
        m = m, e = e, p_win = p_win, xmin = xmin, r_tol = r_tol, grid = grid
      )
    }
    
    # Holt-Laury-consistent participants of this treatment
    scored <- data.table::fread(model_path("mpl", paste0("scored-", tr, ".csv")),
                                encoding = "UTF-8", colClasses = c(pid = "character"))
    consistent_pids <- scored[inconsistent == 0L, pid]
    stopifnot(length(consistent_pids) > 0L, all(consistent_pids %in% pid_full))
    
    for (sample in samples) {
      
      # The consistent sample is the full-sample draws restricted to the
      # consistent participants (valid because r_i is not pooled)
      keep       <- if (sample == "full") rep(TRUE, length(pid_full)) else pid_full %in% consistent_pids
      pid_levels <- pid_full[keep]
      a_star     <- a_star_full[, keep, drop = FALSE]
      iters      <- nrow(a_star)
      N          <- ncol(a_star)
      
      saveRDS(list(pid = pid_levels, a_star_draws = a_star),
              file.path(run, paste0("a-star-draws-", tr, "-", sample, ".rds")))
      
      # ---- Participant classification ----
      P0  <- colMeans(a_star == 0L, na.rm = TRUE)
      tbl <- data.table(pid = pid_levels, P0 = P0, Pplus = 1 - P0,
                        betting_normative = isTRUE(bn[[tr]]))
      
      for (tau in tau_vec) {
        tau_nm <- gsub("\\.", "", sprintf("%.2f", tau))
        if (isTRUE(bn[[tr]])) {
          col  <- paste0("exclude_nonbetters_tau", tau_nm)
          rule <- "exclude non-betters (P0 >= tau)"
          tbl[, (col) := (P0 >= tau)]
        } else {
          col  <- paste0("flag_betters_tau", tau_nm)
          rule <- "flag betters (Pplus >= tau)"
          tbl[, (col) := (Pplus >= tau)]
        }
        count_rows[[paste(tr, sample, tau)]] <- data.table(
          treatment = tr, sample = sample, tau = tau, rule = rule,
          n_participants = nrow(tbl), n_flagged = sum(tbl[[col]])
        )
      }
      
      setorder(tbl, pid)
      fwrite(tbl, file.path(run, paste0("a-star-flags-", tr, "-", sample, ".csv")))
      
      # ---- Consistency checks (each value must equal its threshold of zero) ----
      checks <- c(
        participants_not_unique = length(pid_levels) - length(unique(pid_levels)),
        missing_a_star          = sum(is.na(a_star)),
        a_star_out_of_range     = sum(a_star < 0L | a_star > e, na.rm = TRUE)
      )
      diag_rows[[paste(tr, sample)]] <- data.table(
        tag = paste(tr, sample, sep = "_"), treatment = tr, sample = sample,
        check = names(checks), value = as.numeric(checks),
        threshold = 0, pass = checks == 0
      )
    }
  }
  
  fwrite(rbindlist(count_rows), file.path(run, "a-star-flag-summary.csv"))
  write_diagnostics(run, rbindlist(diag_rows))
  
  model_end(run, cfg, settings = list(
    endowment = e, xmin = xmin, coin_prob = p_win,
    multipliers = design$seq$treatments, samples = samples,
    tau = tau_vec, betting_normative = bn
  ))
}
