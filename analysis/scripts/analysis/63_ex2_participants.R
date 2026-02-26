# ============================================================
# scripts/analysis/63_ex2_participants.R
#   EX2: Participant-level outcomes table (posterior summaries)
#
# STRICT MODE (as requested):
# - NO fallbacks. We ONLY use participant-level generated quantities from Stan:
#     RQ1: mu_b_i
#     RQ2: mu_a_i
#     RQ3: mu_c_i
#     RQ4: mu_h_i
#
# - No global exclusion: output is on UNION(pid1,pid2,pid3,pid4).
# - RQ2 undefined rule (per prereg note): if n_bets < rq2_min_bets -> set mu_a_* = NA
#
# Output (per treatment):
#   data/clean/<ds>/output/ex2_participants_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

ex2_participants <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$plan))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  # RQ2 undefined rule threshold (keep consistent with EX2)
  rq2_min_bets <- cfg$run$ex2_rq2_min_bets
  if (is.null(rq2_min_bets)) rq2_min_bets <- 3L
  rq2_min_bets <- as.integer(rq2_min_bets)
  stopifnot(length(rq2_min_bets) == 1L, rq2_min_bets >= 1L)
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # master for bet counts
  f_master <- file.path(path_clean_ds(ds), "master_sequences.csv")
  stopifnot(file.exists(f_master))
  master <- fread(f_master, encoding = "UTF-8")
  
  reqm <- c("pid", "treat", "stake")
  missm <- setdiff(reqm, names(master))
  if (length(missm) > 0) stop("master_sequences.csv missing: ", paste(missm, collapse = ", "))
  
  master[, pid := as.character(pid)]
  master[, treat := as.character(treat)]
  master[, stake := as.numeric(stake)]
  master[is.na(stake), stake := 0]
  
  summ <- function(x) {
    qs <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    c(
      median = stats::median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      q025   = qs[1],
      q975   = qs[2]
    )
  }
  
  mat_to_tbl <- function(pid, mat, prefix) {
    stopifnot(is.matrix(mat), ncol(mat) == length(pid))
    tmp <- t(apply(mat, 2, summ))
    out <- data.table(pid = as.character(pid))
    out[, paste0(prefix, "_median") := tmp[, "median"]]
    out[, paste0(prefix, "_mean")   := tmp[, "mean"]]
    out[, paste0(prefix, "_q025")   := tmp[, "q025"]]
    out[, paste0(prefix, "_q975")   := tmp[, "q975"]]
    out
  }
  
  # strict getter for mu_*_i from a stanfit extract
  get_mu_i <- function(post, name, pid, label_for_error) {
    if (is.null(post[[name]])) {
      stop(
        "EX2 participants: missing required generated quantity '", name, "' in ", label_for_error, ".\n",
        "Available extracted names (first 60): ",
        paste(head(names(post), 60), collapse = ", ")
      )
    }
    mat <- post[[name]]
    if (!is.matrix(mat)) {
      stop("EX2 participants: '", name, "' in ", label_for_error, " is not a matrix (iters x N).")
    }
    if (ncol(mat) != length(pid)) {
      stop(
        "EX2 participants: '", name, "' column count (", ncol(mat),
        ") != length(pid_levels) (", length(pid), ") in ", label_for_error, "."
      )
    }
    mat
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    # bet counts for RQ2 restriction (stake>0)
    betN <- master[treat == tr & is.finite(stake) & stake > 0,
                   .(n_bets = .N), by = pid]
    setkey(betN, pid)
    
    # ---------- load fits + pid levels ----------
    f1 <- file.path(mod_dir, paste0("rq1_fit_sequences_", tr, ".rds"))
    p1 <- file.path(mod_dir, paste0("rq1_pid_levels_", tr, ".rds"))
    f2 <- file.path(mod_dir, paste0("rq2_fit_sequences_", tr, ".rds"))
    p2 <- file.path(mod_dir, paste0("rq2_pid_levels_", tr, ".rds"))
    f3 <- file.path(mod_dir, paste0("rq3_fit_sequences_", tr, ".rds"))
    p3 <- file.path(mod_dir, paste0("rq3_pid_levels_", tr, ".rds"))
    f4 <- file.path(mod_dir, paste0("rq4_fit_sequences_", tr, ".rds"))
    p4 <- file.path(mod_dir, paste0("rq4_pid_levels_", tr, ".rds"))
    
    stopifnot(file.exists(f1), file.exists(p1),
              file.exists(f2), file.exists(p2),
              file.exists(f3), file.exists(p3),
              file.exists(f4), file.exists(p4))
    
    fit1 <- readRDS(f1); pid1 <- as.character(readRDS(p1))
    fit2 <- readRDS(f2); pid2 <- as.character(readRDS(p2))
    fit3 <- readRDS(f3); pid3 <- as.character(readRDS(p3))
    fit4 <- readRDS(f4); pid4 <- as.character(readRDS(p4))
    
    post1 <- rstan::extract(fit1)
    post2 <- rstan::extract(fit2)
    post3 <- rstan::extract(fit3)
    post4 <- rstan::extract(fit4)
    
    # ---------- STRICT: require mu_*_i ----------
    mu_b_i <- get_mu_i(post1, "mu_b_i", pid1, paste0("RQ1 (tr=", tr, ")"))
    mu_a_i <- get_mu_i(post2, "mu_a_i", pid2, paste0("RQ2 (tr=", tr, ")"))
    mu_c_i <- get_mu_i(post3, "mu_c_i", pid3, paste0("RQ3 (tr=", tr, ")"))
    mu_h_i <- get_mu_i(post4, "mu_h_i", pid4, paste0("RQ4 (tr=", tr, ")"))
    
    rq1_tbl <- mat_to_tbl(pid1, mu_b_i, "mu_b")
    rq2_tbl <- mat_to_tbl(pid2, mu_a_i, "mu_a")
    rq3_tbl <- mat_to_tbl(pid3, mu_c_i, "mu_c")
    rq4_tbl <- mat_to_tbl(pid4, mu_h_i, "mu_h")
    
    # Apply RQ2 min-bets restriction (set NA, do not drop)
    rq2_tbl <- merge(rq2_tbl, betN, by = "pid", all.x = TRUE)
    rq2_tbl[is.na(n_bets), n_bets := 0L]
    rq2_tbl[n_bets < rq2_min_bets,
            c("mu_a_median","mu_a_mean","mu_a_q025","mu_a_q975") :=
              .(NA_real_, NA_real_, NA_real_, NA_real_)]
    
    # ---------- union of participants across models ----------
    pid_all <- sort(unique(c(pid1, pid2, pid3, pid4)))
    out <- data.table(pid = pid_all)
    
    out <- merge(out, rq1_tbl, by = "pid", all.x = TRUE)
    out <- merge(out, rq2_tbl[, c("pid","n_bets","mu_a_median","mu_a_mean","mu_a_q025","mu_a_q975")],
                 by = "pid", all.x = TRUE)
    out <- merge(out, rq3_tbl, by = "pid", all.x = TRUE)
    out <- merge(out, rq4_tbl, by = "pid", all.x = TRUE)
    
    out[, `:=`(dataset = ds, treatment = tr, rq2_min_bets = rq2_min_bets)]
    setcolorder(out, c(
      "dataset","treatment","pid","rq2_min_bets","n_bets",
      "mu_b_median","mu_b_mean","mu_b_q025","mu_b_q975",
      "mu_a_median","mu_a_mean","mu_a_q025","mu_a_q975",
      "mu_c_median","mu_c_mean","mu_c_q025","mu_c_q975",
      "mu_h_median","mu_h_mean","mu_h_q025","mu_h_q975"
    ))
    setorder(out, pid)
    
    f_out <- file.path(out_dir, paste0("ex2_participants_", tr, ".csv"))
    fwrite(out, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- out
  }
  
  invisible(outputs)
}

# Example:
# ex2_participants(cfg)