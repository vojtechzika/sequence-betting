# ============================================================
# 91_synthesis_sequences.R
#
# PURPOSE
#   For every sequence common to all four estimands, computes each
#   headline estimand as a PAIRED posterior difference from its own
#   grand mean: d_s^(t) = mu_s^(t) - grand^(t), read off the SAME
#   iteration t of the SAME fit (no independence assumption, no normal
#   approximation -- the same pairing already used by grand_label/
#   H_delta/T_delta in 13_/23_/33_/43_).
#
#   Reports, per sequence and per estimand: posterior median of d_s,
#   equal-tailed credible interval (level = `ci`), and sig = TRUE iff
#   the CI excludes 0 (plain zero-exclusion test, no ROPE/tolerance band
#   -- deliberately not the delta-band convention used by grand_label/
#   H_delta/T_delta in 13_/23_/33_/43_).
#
#   Estimands: h = heads share (RQ4, post$mu_h / post$hbar, paired in
#   Stan GQ); c = welfare loss (RQ3, post$mu_c, rowMeans as grand mean --
#   no Stan-level grand mean); b = betting rate (RQ1, post$mu_b /
#   post$grand_mean_bet, paired in Stan GQ); a = stake deviation (RQ2,
#   derived mu_s via the alpha+u+b -> delta_bar/sd_star mapping from
#   23_rq2_tables.R::extract_draws(), rowMeans as grand mean). RQ2/RQ3's
#   row-mean grand means are still paired automatically, since they're
#   derived from the same draws matrix as the sequence-level estimate.
#   RQ4 has a single "full" tag -- no full/confirmatory split, no alt
#   model (see 41_/43_ preambles).
#
# INPUT
#   RDS fit files under path_mod, produced upstream by 11_/13_/21_/23_/
#   31_/33_/41_/43_: rq{1,2,3}_fit_sequences_<tr>_<tag>[<sfx>].rds (plus
#   seq_levels / pid_levels / prepared companions for rq2), and
#   rq4_fit_sequences_<tr>_full.rds. Plus path_out/rq{1,2,3}_diagnostics.csv,
#   read to pick the primary vs. alternative-model suffix `<sfx>`. This
#   is the first script in the synthesis chain -- no CSV dependency on
#   another synthesis script. A missing fit for one estimand is not
#   fatal: that estimand is dropped and the rest are still reported.
#
# OUTPUT
#   path_out/synthesis_sequence_diffs_<tr>_<tag>.csv
#   path_out/synthesis_grand_means_<tr>_<tag>.csv
#     columns: estimand, suffix_used, grand_median, grand_lo, grand_hi,
#              grand_alt_mean_of_medians
#   Both are required input to 94_synthesis_mirror_pairs.R and
#   95_synthesis_mechanism_checks.R, which need the default
#   `sequences = NULL` (all 64) run to find all 32 mirror pairs.
#
# NOTES
#   - Only m25 has a confirmatory subset in this design (m19 is not
#     flagged betting_normative, see design_cfg()$a_flags$betting_normative),
#     so this runs a single tr/tag rather than looping over
#     cfg$run$treatment like rq1_tables(cfg)/96_synthesis_regressions.R
#     do. `tr` (default "m25"), `tag` (default "confirmatory"),
#     `sequences` (default NULL = every sequence common to all loaded
#     estimands), `ci` (default 0.95) and `quiet` (default TRUE) are all
#     overridable arguments.
#   - RQ4 ignores `tag` entirely (always fit on the full sample). For
#     RQ1, confirmatory == full in this data (zero normative
#     non-betters), so `b` is numerically identical either way -- a
#     mismatch there is not a tag issue. RQ2/RQ3 can genuinely move
#     between tags (HL-consistency / min-bets exclusions). Passing
#     tr="m19" with tag="confirmatory" finds no confirmatory fits;
#     missing fits are skipped with a warning, not an error.
#   - A warning fires (regardless of `quiet`) if grand_median and
#     grand_alt_mean_of_medians disagree by >0.01 for any estimand --
#     check synthesis_grand_means_<tr>_<tag>.csv first if you see one;
#     it usually means the two are computed from a different definition
#     of "grand mean" and/or a different selected-model suffix than
#     expected.
#   - `quiet = TRUE` (default) suppresses two purely informational
#     warnings: a sequence present in some but not all loaded estimands,
#     and a whole estimand whose fit files couldn't be found. Genuine
#     errors (no fits loadable at all, a requested sequence missing from
#     a loaded estimand, no sequence common to every loaded estimand)
#     always stop(). should_skip() gates the two CSV writes.
#   - Return value: invisible(list(sequences=, grand_means=)).
# ============================================================
synthesis_sequences <- function(cfg, tr = "m25", tag = "confirmatory",
                                sequences = NULL, ci = 0.95, quiet = TRUE) {
  
  # ---- Helper: selected-model suffix (mirrors get_selected_suffix in
  # 13_/23_/33_; "rq" prefix picks the right *_diagnostics.csv and, for
  # rq3, the three-way suffix map). Parameter named `tg`, not `tag`, to
  # avoid data.table NSE self-comparison inside `diag[tag == tg]`. ----
  get_selected_suffix <- function(rq, tr, tg) {
    f_diag <- file.path(path_out, paste0(rq, "_diagnostics.csv"))
    if (!file.exists(f_diag)) return("")
    diag <- fread(f_diag)
    row  <- diag[treatment == tr & tag == tg]
    if (nrow(row) == 0L) return("")
    
    if (rq == "rq3") {
      suffix_map <- list(primary = "", gamma_only = "_gamma", alternative = "_alt")
      sfx <- suffix_map[[row$selected_model]]
      return(if (is.null(sfx)) "" else sfx)
    }
    if (row$selected_model == "alternative") "_alt" else ""
  }
  
  # ---- Helper: paired difference from a grand-mean draw vector ----
  diff_from_grand <- function(mu_draws, grand_draws, seq_levels, ci) {
    stopifnot(is.matrix(mu_draws), length(grand_draws) == nrow(mu_draws))
    
    d_mat <- mu_draws - grand_draws          # recycled down rows -> paired by draw t
    
    lo_p <- (1 - ci) / 2
    hi_p <- 1 - lo_p
    d_lo <- apply(d_mat, 2, quantile, probs = lo_p)
    d_hi <- apply(d_mat, 2, quantile, probs = hi_p)
    sig  <- (d_lo > 0) | (d_hi < 0)   # CI excludes 0, no ROPE/tolerance band
    
    data.table(
      sequence = seq_levels,
      d_median = apply(d_mat, 2, median),
      d_lo     = d_lo,
      d_hi     = d_hi,
      sig      = sig
    )
  }
  
  # ---- Helper: RQ1 betting rate (b) -- mu_b / grand_mean_bet, paired
  # in Stan GQ ----
  load_rq1_draws <- function(tr, tg) {
    sfx   <- get_selected_suffix("rq1", tr, tg)
    f_fit <- file.path(path_mod, paste0("rq1_fit_sequences_", tr, "_", tg, sfx, ".rds"))
    f_seq <- file.path(path_mod, paste0("rq1_seq_levels_",    tr, "_", tg, sfx, ".rds"))
    if (!file.exists(f_fit) || !file.exists(f_seq)) return(NULL)
    
    seq_levels <- as.character(readRDS(f_seq))
    post <- rstan::extract(readRDS(f_fit))
    stopifnot(!is.null(post$mu_b), !is.null(post$grand_mean_bet))
    
    list(mu_draws    = post$mu_b,
         grand_draws = as.numeric(post$grand_mean_bet),
         seq_levels  = seq_levels,
         suffix      = sfx)
  }
  
  # ---- Helper: RQ2 stake deviation (a) -- derived mu_s (no Stan-level
  # grand mean), reusing the exact alpha/u/b -> delta_bar/sd_star
  # mapping from 23_rq2_tables.R::extract_draws() ----
  load_rq2_draws <- function(cfg, tr, tg) {
    sfx    <- get_selected_suffix("rq2", tr, tg)
    f_fit  <- file.path(path_mod, paste0("rq2_fit_sequences_", tr, "_", tg, sfx, ".rds"))
    f_pid  <- file.path(path_mod, paste0("rq2_pid_levels_",    tr, "_", tg, sfx, ".rds"))
    f_seq  <- file.path(path_mod, paste0("rq2_seq_levels_",    tr, "_", tg, sfx, ".rds"))
    f_prep <- file.path(path_mod, paste0("rq2_prepared_",      tr, "_", tg, sfx, ".rds"))
    if (!all(file.exists(c(f_fit, f_pid, f_seq, f_prep)))) return(NULL)
    
    e          <- as.numeric(cfg$design$seq$endowment)
    pid_levels <- as.character(readRDS(f_pid))
    seq_levels <- as.character(readRDS(f_seq))
    prep       <- readRDS(f_prep)
    prep[, pid := as.character(pid)]
    
    pid_map <- prep[, .(delta_bar = delta_bar[1], sd_star = sd_star[1]), by = pid]
    setkey(pid_map, pid)
    pid_map <- pid_map[.(pid_levels)]
    stopifnot(!anyNA(pid_map$pid))
    
    post <- rstan::extract(readRDS(f_fit))
    stopifnot(!is.null(post$alpha), !is.null(post$u), !is.null(post$b))
    
    alpha <- as.numeric(post$alpha)
    u     <- post$u
    b     <- post$b
    K     <- length(alpha)
    S     <- length(seq_levels)
    
    stopifnot(nrow(u) == K, ncol(u) == length(pid_levels), nrow(b) == K, ncol(b) == S)
    
    eta_base   <- sweep(u, 1, alpha, "+")
    mu_s_draws <- matrix(NA_real_, nrow = K, ncol = S)
    for (s in seq_len(S)) {
      eta_mat         <- eta_base + b[, s]
      mapped          <- sweep(eta_mat, 2, as.numeric(pid_map$sd_star),   "*")
      mapped          <- sweep(mapped,  2, as.numeric(pid_map$delta_bar), "+")
      mu_s_draws[, s] <- rowMeans(mapped) / e
    }
    
    list(mu_draws    = mu_s_draws,
         grand_draws = rowMeans(mu_s_draws),   # grand mean over ALL 64 sequences
         seq_levels  = seq_levels,
         suffix      = sfx)
  }
  
  # ---- Helper: RQ3 welfare loss (c) -- mu_c, no Stan-level grand mean ----
  load_rq3_draws <- function(tr, tg) {
    sfx   <- get_selected_suffix("rq3", tr, tg)
    f_fit <- file.path(path_mod, paste0("rq3_fit_sequences_", tr, "_", tg, sfx, ".rds"))
    f_seq <- file.path(path_mod, paste0("rq3_seq_levels_",    tr, "_", tg, sfx, ".rds"))
    if (!file.exists(f_fit) || !file.exists(f_seq)) return(NULL)
    
    seq_levels <- as.character(readRDS(f_seq))
    post <- rstan::extract(readRDS(f_fit))
    stopifnot(!is.null(post$mu_c))
    
    mu_s_draws <- post$mu_c
    list(mu_draws    = mu_s_draws,
         grand_draws = rowMeans(mu_s_draws),
         seq_levels  = seq_levels,
         suffix      = sfx)
  }
  
  # ---- Helper: RQ4 heads share (h) -- mu_h / hbar, paired in Stan GQ.
  # Single "full" tag, no alt model (see 41_/43_ preambles) ----
  load_rq4_draws <- function(tr) {
    f_fit <- file.path(path_mod, paste0("rq4_fit_sequences_", tr, "_full.rds"))
    f_seq <- file.path(path_mod, paste0("rq4_seq_levels_",    tr, "_full.rds"))
    if (!file.exists(f_fit) || !file.exists(f_seq)) return(NULL)
    
    seq_levels <- as.character(readRDS(f_seq))
    post <- rstan::extract(readRDS(f_fit))
    stopifnot(!is.null(post$mu_h), !is.null(post$hbar))
    
    list(mu_draws    = post$mu_h,
         grand_draws = as.numeric(post$hbar),
         seq_levels  = seq_levels,
         suffix      = "")   # RQ4 has no alt model, always "full" tag
  }
  
  # ---- Main body ----
  specs <- list(
    h = load_rq4_draws(tr),
    c = load_rq3_draws(tr, tag),
    b = load_rq1_draws(tr, tag),
    a = load_rq2_draws(cfg, tr, tag)
  )
  
  loaded <- specs[!vapply(specs, is.null, logical(1))]
  if (length(loaded) == 0L) {
    stop("synthesis_sequences: no fits could be loaded for tr=", tr, ", tag=", tag, ".")
  }
  
  if (is.null(sequences)) {
    # Report every sequence common to all four loaded estimands (typically
    # all 64) -- intersect rather than union so the merged table never
    # silently drops a column for a sequence one estimand happens to lack.
    seq_sets <- lapply(loaded, `[[`, "seq_levels")
    common   <- Reduce(intersect, seq_sets)
    all_seen <- Reduce(union, seq_sets)
    dropped  <- setdiff(all_seen, common)
    if (length(dropped) > 0L && !quiet) {
      warning("synthesis_sequences (tr=", tr, ", tag=", tag,
              "): sequence(s) present in some but not all loaded estimands, excluded: ",
              paste(dropped, collapse = ", "))
    }
    if (length(common) == 0L) {
      stop("synthesis_sequences: no sequence is common to all loaded estimands ",
           "for tr=", tr, ", tag=", tag, ".")
    }
    # Preserve the ordering of the first loaded estimand's seq_levels (all
    # four should be the same 64-sequence set/order in this design; this
    # just guards against a stray reordering upstream).
    sequences <- seq_sets[[1]][seq_sets[[1]] %in% common]
  }
  
  out        <- NULL
  grand_rows <- list()
  
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    if (is.null(sp)) {
      if (!quiet) {
        warning("synthesis_sequences: missing fit for estimand '", nm,
                "' (tr=", tr, ", tag=", tag, "). Skipping.")
      }
      next
    }
    missing_seq <- setdiff(sequences, sp$seq_levels)
    if (length(missing_seq) > 0L) {
      stop("synthesis_sequences: sequence(s) not found for '", nm, "': ",
           paste(missing_seq, collapse = ", "))
    }
    idx <- match(sequences, sp$seq_levels)
    
    tbl <- diff_from_grand(sp$mu_draws[, idx, drop = FALSE], sp$grand_draws,
                           sequences, ci)
    setnames(tbl, c("d_median", "d_lo", "d_hi", "sig"),
             paste0(nm, c("_d_median", "_d_lo", "_d_hi", "_sig")))
    
    out <- if (is.null(out)) tbl else merge(out, tbl, by = "sequence", sort = FALSE)
    
    # grand_median: posterior median of the paired grand-mean draws (what
    # the diff table above actually uses). grand_alt_mean_of_medians: mean,
    # across all 64 sequences, of each sequence's OWN posterior median --
    # a different (non-Bayesian) summary a downstream table might use
    # instead. Disagreement beyond rounding flags a mismatch (see NOTES).
    grand_rows[[nm]] <- data.table(
      estimand    = nm,
      suffix_used = sp$suffix,
      grand_median               = median(sp$grand_draws),
      grand_lo                   = as.numeric(quantile(sp$grand_draws, (1 - ci) / 2)),
      grand_hi                   = as.numeric(quantile(sp$grand_draws, 1 - (1 - ci) / 2)),
      grand_alt_mean_of_medians  = mean(apply(sp$mu_draws, 2, median))
    )
  }
  
  stopifnot(!is.null(out))
  out <- out[match(sequences, sequence)]   # restore requested display order
  grand_tbl <- rbindlist(grand_rows)
  
  gap <- abs(grand_tbl$grand_median - grand_tbl$grand_alt_mean_of_medians)
  bad <- grand_tbl$estimand[gap > 0.01]
  if (length(bad) > 0L) {
    warning("synthesis_sequences (tr=", tr, ", tag=", tag, "): grand_median ",
            "and grand_alt_mean_of_medians disagree by >0.01 for estimand(s): ",
            paste(bad, collapse = ", "),
            ". See synthesis_grand_means_", tr, "_", tag, ".csv -- this usually ",
            "means the two are computed from a different definition of 'grand ",
            "mean' (paired population draws vs. average of per-sequence medians) ",
            "and/or a different selected model suffix than expected.")
  }
  
  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  
  f_out <- file.path(path_out, paste0("synthesis_sequence_diffs_", tr, "_", tag, ".csv"))
  if (!should_skip(f_out, cfg, "output", paste0("Synthesis sequence diffs (", tr, "/", tag, ")"))) {
    fwrite(out, f_out)
    msg("Saved: ", f_out)
  }
  
  f_grand <- file.path(path_out, paste0("synthesis_grand_means_", tr, "_", tag, ".csv"))
  if (!should_skip(f_grand, cfg, "output", paste0("Synthesis grand means (", tr, "/", tag, ")"))) {
    fwrite(grand_tbl, f_grand)
    msg("Saved: ", f_grand)
  }
  
  invisible(list(sequences = out, grand_means = grand_tbl))
}