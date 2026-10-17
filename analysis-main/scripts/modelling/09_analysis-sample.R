# ============================================================
# 09_analysis-sample.R
#
# PURPOSE
#   Defines, in one place, which participants and trials enter each
#   analysis. The script builds one trial-level table and, for every
#   named sample, a logical membership column. Models obtain their data
#   with analysis_trials(treatment, sample) and never define a sample
#   themselves; a sample that is not defined for a treatment stops the
#   caller with an error that lists the defined samples.
#
# SAMPLES (participants x trials)
#   participants  full          all participants of the treatment
#                 confirmatory  participants kept by the a* classification,
#                               in the Holt-Laury sample given by
#                               hl_sample(cfg); only for treatments with
#                               cfg$design$a_flags$betting_normative = TRUE
#   trials        all       all trials of the participant
#                 bets      trials with a bet (stake > 0)
#                 bets-min  trials with a bet (stake > 0) of participants
#                           with at least cfg$design$rq2$min_bets bets
#   Defined samples: rq1-full, rq1-confirmatory, rq2-full,
#   rq2-confirmatory, rq3-full, rq3-confirmatory (RQ1 and RQ3: all trials;
#   RQ2: bets-min), rq4-full (bets; no restriction on participants), and
#   synthesis-full (all trials; the model derives its bet and side
#   equations from the bet indicator and the side). A new sample is added to
#   the specification table in analysis_sample() and the folder is rebuilt.
#   The exploratory analyses EX1 and EX2 use participant-level quantities
#   built from the posterior draws of the RQ models, so their participants
#   are those of the RQ fits they build on and they have no sample of their own.
#
# FUNCTIONS
#   analysis_sample(cfg)           builds the table and the definitions
#   analysis_trials(tr, sample)    trials of a sample for one treatment
#   sample_definition(tr, sample)  its definition and participant and trial counts
#
# INPUT
#   path_der/master-sequences.csv
#   model "mpl": a-star/a-star-draws-<tr>-full.rds   (a* draws)
#               a-star/a-star-flags-<tr>-<sample>.csv (confirmatory set)
#
# OUTPUT (subfolder sample/ of the selected MPL model folder)
#   analysis-trials.csv    one row per trial: pid, treat, trial_id, seq_no,
#                          seq, block, pos, side, stake, screen_ms, bet,
#                          block_c, a_star (posterior mean), n_bets, and one logical
#                          column per sample (hyphens replaced by underscores)
#   sample-definitions.csv sample, treatment, participant and trial rule,
#                          number of participants and trials
#   diagnostics.csv        integrity checks
#   inputs.csv, run-info.txt
#
# NOTES
#   - a* of a participant is the same in every sample (posterior mean of the
#     draws of the full-sample fit)
#   - The subfolder is built as sample.tmp and renamed when complete; an
#     existing folder is never overwritten (delete it to rerun, for example
#     after changing min_bets or the a* step)
# ============================================================

analysis_sample <- function(cfg) {
  
  # sample subfolder inside the selected MPL model folder
  mpl_dir <- model_path("mpl")
  run     <- model_begin(file.path(basename(mpl_dir), "sample"))
  if (is.null(run)) return(invisible(NULL))
  
  design     <- cfg$design
  treatments <- names(design$seq$treatments)
  bn         <- design$a_flags$betting_normative
  tau        <- as.numeric(design$a_flags$tau[1])
  min_bets   <- as.integer(design$rq2$min_bets)
  n_seq      <- as.integer(design$seq$seq_n)
  stopifnot(length(min_bets) == 1L, min_bets >= 1L, length(n_seq) == 1L)
  
  # ---- Sample specification ----
  spec <- data.table(
    sample       = c("rq1-full", "rq1-confirmatory",
                     "rq2-full", "rq2-confirmatory",
                     "rq3-full", "rq3-confirmatory",
                     "rq4-full", "synthesis-full"),
    participants = c("full", "confirmatory",
                     "full", "confirmatory",
                     "full", "confirmatory",
                     "full", "full"),
    trials       = c("all", "all",
                     "bets-min", "bets-min",
                     "all", "all",
                     "bets", "all")
  )
  
  # ---- Trials ----
  infile <- file.path(path_der, "master-sequences.csv")
  stopifnot(file.exists(infile))
  keep_cols <- c("pid", "treat", "trial_id", "seq_no", "seq", "block", "pos",
                 "side", "stake", "screen_ms")
  trials <- fread(infile, encoding = "UTF-8")
  stopifnot(all(keep_cols %in% names(trials)))
  trials <- trials[, ..keep_cols]
  trials[, pid   := as.character(pid)]
  trials[, treat := as.character(treat)]
  trials[, seq   := as.character(seq)]
  trials[, stake := as.numeric(stake)]
  trials[, screen_ms := as.numeric(screen_ms)]
  
  tbl_list <- list()
  defs     <- list()
  
  for (tr in treatments) {
    
    d <- trials[treat == tr]
    
    # EU-optimal stake: posterior mean of the a* draws of the full-sample fit
    ast <- readRDS(model_path("mpl", file.path("a-star",
                                               paste0("a-star-draws-", tr, "-full.rds"))))
    stopifnot(is.list(ast), all(c("pid", "a_star_draws") %in% names(ast)))
    a_map <- data.table(pid    = as.character(ast$pid),
                        a_star = as.numeric(colMeans(ast$a_star_draws)))
    d <- merge(d, a_map, by = "pid", all.x = TRUE, sort = FALSE)
    
    d[, bet     := as.integer(!is.na(stake) & stake > 0)]
    d[, n_bets  := sum(bet), by = pid]
    d[, block_c := as.numeric(block) - 2.5]
    
    conf_pids <- if (isTRUE(bn[[tr]])) {
      a_star_pid_set(tr, hl_sample(cfg), tau, "keep")
    } else NULL
    
    for (i in seq_len(nrow(spec))) {
      s <- spec[i]
      if (s$participants == "confirmatory" && is.null(conf_pids)) next
      
      in_part  <- if (s$participants == "full") rep(TRUE, nrow(d)) else d$pid %in% conf_pids
      in_trial <- switch(s$trials,
                         "all"      = rep(TRUE, nrow(d)),
                         "bets"     = d$bet == 1L,
                         "bets-min" = d$bet == 1L & d$n_bets >= min_bets)
      member   <- in_part & in_trial
      set(d, j = gsub("-", "_", s$sample), value = member)
      
      defs[[paste(tr, s$sample)]] <- data.table(
        sample = s$sample, treatment = tr,
        participants = s$participants, trials = s$trials,
        min_bets = if (s$trials == "bets-min") min_bets else NA_integer_,
        n_participants = uniqueN(d$pid[member]), n_trials = sum(member)
      )
    }
    tbl_list[[tr]] <- d
  }
  
  out  <- rbindlist(tbl_list, fill = TRUE)
  defs <- rbindlist(defs)
  
  # a treatment without a sample has FALSE in its membership column
  for (col in gsub("-", "_", spec$sample)) {
    set(out, i = which(is.na(out[[col]])), j = col, value = FALSE)
  }
  
  setcolorder(out, c("pid", "treat", "trial_id", "seq_no", "seq", "block", "pos",
                     "side", "stake", "screen_ms", "bet", "block_c", "a_star", "n_bets"))
  
  fwrite(out,  file.path(run, "analysis-trials.csv"))
  fwrite(defs, file.path(run, "sample-definitions.csv"))
  
  # ---- Integrity checks (each value must equal its threshold) ----
  per_pid   <- out[, .N, by = pid]
  col_names <- gsub("-", "_", spec$sample)
  conf_cols <- col_names[spec$participants == "confirmatory"]
  bet_cols  <- col_names[spec$trials %in% c("bets", "bets-min")]
  checks <- rbindlist(list(
    data.table(check = "participants_without_all_sequences",
               value = sum(per_pid$N != n_seq)),
    data.table(check = "duplicated_trials",
               value = sum(duplicated(out[, .(pid, trial_id)]))),
    data.table(check = "missing_a_star", value = sum(is.na(out$a_star))),
    data.table(check = "empty_samples", value = sum(defs$n_trials == 0L)),
    data.table(check = "confirmatory_not_in_full", value = sum(vapply(conf_cols, function(cc)
      sum(out[[cc]] & !out[[sub("confirmatory", "full", cc)]]), numeric(1)))),
    data.table(check = "bets_sample_with_non_bets", value = sum(vapply(bet_cols, function(cc)
      sum(out[[cc]] & out$bet == 0L), numeric(1))))
  ))
  checks[, `:=`(tag = "sample", treatment = NA_character_, sample = NA_character_,
                threshold = 0, pass = value == 0)]
  setcolorder(checks, c("tag", "treatment", "sample", "check", "value", "threshold", "pass"))
  write_diagnostics(run, checks)
  
  model_end(run, cfg, settings = list(min_bets = min_bets, tau = tau,
                                      hl_sample = hl_sample(cfg),
                                      samples = spec$sample))
}

# Definition of a sample for one treatment: a one-row table with the
# participant and trial counts. Stops if the sample is not defined for the
# treatment, listing the defined samples.
sample_definition <- function(tr, sample) {
  stopifnot(length(tr) == 1L, length(sample) == 1L)
  
  defs <- data.table::fread(model_path("mpl", file.path("sample", "sample-definitions.csv")),
                            encoding = "UTF-8", colClasses = "character")
  # (the logical vector is built outside data.table's scope, where the
  # column "sample" would mask the argument of the same name)
  is_hit <- defs$sample == sample & defs$treatment == tr
  hit    <- defs[which(is_hit)]
  if (nrow(hit) == 0L) {
    stop("Sample '", sample, "' is not defined for treatment '", tr, "'. Defined: ",
         paste(unique(paste0(defs$sample, " (", defs$treatment, ")")), collapse = ", "),
         call. = FALSE)
  }
  hit[, `:=`(n_participants = as.integer(n_participants), n_trials = as.integer(n_trials))]
  hit
}

# Trials of a sample for one treatment (all columns of analysis-trials.csv).
analysis_trials <- function(tr, sample) {
  sample_definition(tr, sample)          # stops if the sample is not defined
  
  tbl <- data.table::fread(model_path("mpl", file.path("sample", "analysis-trials.csv")),
                           encoding = "UTF-8",
                           colClasses = c(pid = "character", treat = "character",
                                          seq = "character"))
  col    <- gsub("-", "_", sample)
  in_set <- tbl$treat == tr & tbl[[col]]
  tbl[which(in_set)]
}
