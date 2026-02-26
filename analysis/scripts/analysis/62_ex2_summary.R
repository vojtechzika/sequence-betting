# ============================================================
# scripts/analysis/62_ex2_summary.R
#   EX2: Summarize Stan results (pooled + outcome-specific slopes)
#
# Inputs (per dataset x treatment):
#   data/clean/<ds>/models/ex2_fit_<tr>.rds
#
# Output:
#   data/clean/<ds>/output/ex2_summary_<tr>.csv
# ============================================================

library(data.table)
library(rstan)

ex2_summary <- function(cfg) {
  
  stopifnot(is.list(cfg), !is.null(cfg$run), !is.null(cfg$plan))
  stopifnot(!is.null(cfg$run$dataset))
  stopifnot(!is.null(cfg$plan$by))
  
  ds <- as.character(cfg$run$dataset)
  tr_vec <- unique(as.character(cfg$plan$by))
  stopifnot(length(tr_vec) > 0L, all(nzchar(tr_vec)))
  
  mod_dir <- path_mod_ds(ds)
  out_dir <- path_out_ds(ds)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # robust summary with fixed q025/q975 names
  summ <- function(x) {
    qs <- stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    c(
      median = stats::median(x, na.rm = TRUE),
      mean   = mean(x, na.rm = TRUE),
      q025   = qs[1],
      q975   = qs[2],
      p_gt0  = mean(x > 0, na.rm = TRUE)
    )
  }
  
  # like summ(), but p_gt0 is meaningless (e.g., tau, sigma)
  summ_nop <- function(x) {
    s <- summ(x)
    s["p_gt0"] <- NA_real_
    s
  }
  
  outputs <- list()
  
  for (tr in tr_vec) {
    
    f_fit <- file.path(mod_dir, paste0("ex2_fit_", tr, ".rds"))
    if (!file.exists(f_fit)) {
      warning("EX2 summary: missing fit for treatment='", tr, "': ", f_fit)
      next
    }
    
    fit_obj <- readRDS(f_fit)
    
    # allow either raw stanfit, or list wrapper
    fit <- fit_obj
    if (is.list(fit_obj) && !inherits(fit_obj, "stanfit")) {
      if (!is.null(fit_obj$fit) && inherits(fit_obj$fit, "stanfit")) fit <- fit_obj$fit
    }
    if (!inherits(fit, "stanfit")) {
      stop("EX2 summary: object in ", f_fit, " is neither stanfit nor list(fit=stanfit,...).")
    }
    
    post <- rstan::extract(fit)
    pn <- names(post)
    
    # ---- infer K from outcome-specific slope matrices beta_*_k ----
    K <- NA_integer_
    if ("beta_opt_k" %in% pn && is.matrix(post$beta_opt_k)) K <- ncol(post$beta_opt_k)
    if (is.na(K) && "beta_rt_k" %in% pn && is.matrix(post$beta_rt_k)) K <- ncol(post$beta_rt_k)
    if (is.na(K) && "beta_r_k" %in% pn  && is.matrix(post$beta_r_k))  K <- ncol(post$beta_r_k)
    if (is.na(K)) {
      # Nothing outcome-specific available; still allow pooled-only summary
      K <- 0L
    }
    
    # prereg outcome labels (fixed order). If K != 4, fall back to k1..kK.
    outcomes <- if (K == 4L) c("b","a","c","h") else paste0("k", seq_len(K))
    
    rows <- list()
    
    # ---- pooled slopes (required) ----
    if ("beta_opt_bar" %in% pn) {
      rows[["pooled_beta_opt_bar"]] <- data.table(
        block = "pooled", outcome = NA_character_, term = "beta_opt_bar",
        t(summ(post$beta_opt_bar))
      )
    } else {
      warning("EX2 summary: missing beta_opt_bar in fit (tr='", tr, "').")
    }
    
    if ("beta_rt_bar" %in% pn) {
      rows[["pooled_beta_rt_bar"]] <- data.table(
        block = "pooled", outcome = NA_character_, term = "beta_rt_bar",
        t(summ(post$beta_rt_bar))
      )
    } else {
      warning("EX2 summary: missing beta_rt_bar in fit (tr='", tr, "').")
    }
    
    if ("beta_r_bar" %in% pn) {
      rows[["pooled_beta_r_bar"]] <- data.table(
        block = "pooled", outcome = NA_character_, term = "beta_r_bar",
        t(summ(post$beta_r_bar))
      )
    }
    
    # ---- pooling scales (tau_*) ----
    if ("tau_opt" %in% pn) rows[["pooled_tau_opt"]] <- data.table(block="pooled", outcome=NA_character_, term="tau_opt", t(summ_nop(post$tau_opt)))
    if ("tau_rt"  %in% pn) rows[["pooled_tau_rt"]]  <- data.table(block="pooled", outcome=NA_character_, term="tau_rt",  t(summ_nop(post$tau_rt)))
    if ("tau_r"   %in% pn) rows[["pooled_tau_r"]]   <- data.table(block="pooled", outcome=NA_character_, term="tau_r",   t(summ_nop(post$tau_r)))
    
    # ---- outcome-specific slopes (beta_*_k) ----
    add_matK <- function(param, term_name) {
      if (!(param %in% pn)) return(invisible(NULL))
      mat <- post[[param]]
      if (!is.matrix(mat)) return(invisible(NULL))
      if (ncol(mat) != K) return(invisible(NULL))
      for (k in seq_len(K)) {
        rows[[paste0(param, "_", k)]] <<- data.table(
          block = "outcome",
          outcome = outcomes[k],
          term = term_name,
          t(summ(mat[, k]))
        )
      }
      invisible(NULL)
    }
    
    if (K > 0) {
      add_matK("beta_opt_k", "beta_opt")
      add_matK("beta_rt_k",  "beta_rt")
      add_matK("beta_r_k",   "beta_r")
    }
    
    # ---- outcome-specific residual scales (sigma_k) if present ----
    # Some versions name this sigma_k; others sigma (vector[K]).
    if (K > 0) {
      if ("sigma_k" %in% pn && is.matrix(post$sigma_k) && ncol(post$sigma_k) == K) {
        for (k in seq_len(K)) {
          rows[[paste0("sigma_k_", k)]] <- data.table(
            block="outcome", outcome=outcomes[k], term="sigma",
            t(summ_nop(post$sigma_k[, k]))
          )
        }
      } else if ("sigma" %in% pn) {
        sig <- post$sigma
        if (is.matrix(sig) && ncol(sig) == K) {
          for (k in seq_len(K)) {
            rows[[paste0("sigma_", k)]] <- data.table(
              block="outcome", outcome=outcomes[k], term="sigma",
              t(summ_nop(sig[, k]))
            )
          }
        }
      }
    }
    
    tbl <- rbindlist(rows, fill = TRUE)
    tbl[, `:=`(dataset = ds, treatment = tr)]
    setcolorder(tbl, c("dataset","treatment","block","outcome","term","median","mean","q025","q975","p_gt0"))
    
    f_out <- file.path(out_dir, paste0("ex2_summary_", tr, ".csv"))
    fwrite(tbl, f_out)
    msg("Saved: ", f_out)
    
    outputs[[tr]] <- tbl
  }
  
  invisible(outputs)
}

# Example:
# ex2_summary(cfg)