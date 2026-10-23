# ============================================================
# 102_prior-sensitivity.R
#
# PURPOSE
#   Compares every RQ1-RQ4 model fitted under the main priors with its
#   refit under the broadened priors (prior_broad in settings.R, written
#   to the folder <model folder>-broad), and evaluates the robustness
#   criteria below. Estimands, thresholds and labels are those of the
#   reporting scripts 11, 21, 31 and 41 and must be kept identical to them
#   (labels are recomputed here from the draws, at the main rho).
#
# ROBUSTNESS CRITERIA (fixed before the broad fits were inspected;
#   not preregistered)
#   1. Shift: the posterior mean of every estimand changes by less than
#      0.1 posterior SD of the main fit.
#   2. Inference: no 95% interval changes its position relative to the
#      benchmark (entirely below, overlapping, or entirely above the
#      threshold(s) that define the label: 1 - rho for RQ1, -rho and rho for
#      RQ2, rho for RQ3, hbar -/+ rho for RQ4, with hbar at its posterior mean).
#   3. Labels: no label moves by more than one step on the ordered label
#      scale (RQ2 and RQ4 are signed, so a change of direction counts as the
#      sum of the steps to neutral).
#   The criteria are evaluated for sequences (estimand mu_s). The same
#   statistics are reported for participants (estimand mu_i) as
#   description. A comparison is valid only if the broad fit passed the
#   convergence checks of its diagnostics.csv (column converged_broad).
#
# ESTIMANDS (names in the fits)
#   RQ1 mu_b, mu_b_i | RQ2 mu_a, mu_a_i | RQ3 mu_c, mu_c_i | RQ4 mu_h, mu_h_i, hbar
#
# INPUT
#   data/models/<folder>/ and data/models/<folder>-broad/ for the folders
#   listed in the table of models below: fit-<tr>-<sample>.rds,
#   pid-levels-<tr>-<sample>.rds, seq-levels-<tr>-<sample>.rds,
#   diagnostics.csv (broad folder). Pairs without a broad folder are skipped.
#   selected-models.csv (marks the model used in the paper).
#
# OUTPUT (data/output/csv/)
#   prior-sensitivity-summary.csv      one row per model, sample and level:
#                                      the three criteria and the verdict
#   prior-sensitivity-units.csv        one row per sequence or participant:
#                                      means, SD, interval position, label
#                                      under both priors
#   prior-sensitivity-parameters.csv   scalar parameters under both priors
#
# NOTES
#   - Run after the models and the broad passes (cfg$run$prior_sensitivity)
#     are complete; output files are overwritten at every run
#   - The label rules duplicate those of the table scripts; if those change,
#     change the functions score_* below
# ============================================================

prior_sensitivity <- function(cfg) {
  
  # Nothing to compare unless broad fits were requested (cfg$run$prior_sensitivity)
  if (length(cfg$run$prior_sensitivity) == 0L) {
    msg("Prior sensitivity: no broad fits requested (cfg$run$prior_sensitivity is empty); skipped.")
    return(invisible(NULL))
  }
  
  design <- cfg$design
  
  # ---- Criteria ----
  thr_shift <- 0.1      # criterion 1: shift in posterior SD of the main fit
  thr_step  <- 1L       # criterion 3: largest admissible move on the label scale
  
  # ---- Models to compare (folder names; the broad twin has suffix -broad) ----
  models <- data.table(
    rq     = c("rq1", "rq2", "rq2", "rq2", "rq3", "rq3", "rq3", "rq4"),
    folder = c("rq1-bernoulli", "rq2-gaussian", "rq2-one-inflated",
               "rq2-one-inflated-seq", "rq3-hurdle-gamma", "rq3-gamma",
               "rq3-gaussian", "rq4-bernoulli"))
  sel <- read_selected_models()
  sel_folders <- if (is.null(sel)) character(0) else sel$folder
  models[, selected := folder %in% sel_folders]
  
  # ---- Estimand names and the main rho of each research question ----
  est <- list(rq1 = c(s = "mu_b", i = "mu_b_i"), rq2 = c(s = "mu_a", i = "mu_a_i"),
              rq3 = c(s = "mu_c", i = "mu_c_i"), rq4 = c(s = "mu_h", i = "mu_h_i"))
  rho <- vapply(c("rq1", "rq2", "rq3", "rq4"),
                function(r) as.numeric(design[[r]]$rho)[1], numeric(1))
  
  # ---- Exceedance probabilities and benchmark thresholds ----
  # pa, pb: probabilities that define the label (pb only for RQ2 and RQ4);
  # lo, hi: thresholds against which the 95% interval is placed
  exceed <- function(rq, mu, post) {
    r <- rho[[rq]]
    if (rq == "rq1") list(pa = colMeans(mu < 1 - r), pb = NA_real_, lo = 1 - r, hi = 1 - r)
    else if (rq == "rq2") list(pa = colMeans(mu < -r), pb = colMeans(mu > r), lo = -r, hi = r)
    else if (rq == "rq3") list(pa = colMeans(mu > r), pb = NA_real_, lo = r, hi = r)
    else {
      hbar <- as.numeric(post$hbar)
      list(pa = colMeans(mu > hbar + r), pb = colMeans(mu < hbar - r),
           lo = mean(hbar) - r, hi = mean(hbar) + r)
    }
  }
  
  # ---- Labels on an ordered integer scale (0 = neutral) ----
  cuts  <- list(sequence = c(0.50, 0.80, 0.95), participant = c(0.75, 0.90, 0.95))
  tier  <- function(p, cu) as.integer(p >= cu[1]) + as.integer(p >= cu[2]) + as.integer(p >= cu[3])
  score <- function(rq, level, pa, pb) {
    cu <- cuts[[level]]
    if (rq %in% c("rq1", "rq3")) return(tier(pa, cu))
    if (rq == "rq2") {
      under <- pa >= cu[1] & pa >= pb
      over  <- pb >= cu[1] & pb >  pa
      return(ifelse(under, -tier(pa, cu), ifelse(over, tier(pb, cu), 0L)))
    }
    # RQ4: strength by the larger of the two probabilities, floor at cu[1]
    m <- pmax(pa, pb)
    t <- 1L + as.integer(m >= cu[2]) + as.integer(m >= cu[3])
    ifelse(m < cu[1], 0L, ifelse(pa >= pb, t, -t))
  }
  tier_names <- list(sequence = c("weak", "moderate", "strong"),
                     participant = c("leaning", "likely", "solid"))
  label_of <- function(rq, level, s) {
    nm <- tier_names[[level]]
    a  <- abs(s)
    out <- rep("neutral", length(s))
    k <- a > 0
    out[k] <- if (rq %in% c("rq1", "rq3")) nm[a[k]]
              else if (rq == "rq2") paste0(ifelse(s[k] < 0, "under_", "over_"), nm[a[k]])
              else paste0(nm[a[k]], "_", ifelse(s[k] > 0, "head", "tail"))
    out
  }
  position <- function(q025, q975, lo, hi) ifelse(q975 < lo, -1L, ifelse(q025 > hi, 1L, 0L))
  
  # ---- One level (sequences or participants) of one fit pair ----
  compare_level <- function(rq, level, ids, mu_m, mu_b, post_m, post_b) {
    em <- exceed(rq, mu_m, post_m)
    eb <- exceed(rq, mu_b, post_b)
    sd_m <- apply(mu_m, 2, sd)
    q <- function(m) list(mean = colMeans(m), q025 = apply(m, 2, quantile, probs = 0.025),
                          q975 = apply(m, 2, quantile, probs = 0.975))
    qm <- q(mu_m); qb <- q(mu_b)
    sm <- score(rq, level, em$pa, em$pb)
    sb <- score(rq, level, eb$pa, eb$pb)
    data.table(
      level = level, id = ids,
      mean_main = qm$mean, mean_broad = qb$mean, sd_main = sd_m,
      shift_sd = (qb$mean - qm$mean) / sd_m,
      pos_main = position(qm$q025, qm$q975, em$lo, em$hi),
      pos_broad = position(qb$q025, qb$q975, eb$lo, eb$hi),
      pa_main = em$pa, pa_broad = eb$pa, pb_main = em$pb, pb_broad = eb$pb,
      label_main = label_of(rq, level, sm), label_broad = label_of(rq, level, sb),
      step = abs(sb - sm))
  }
  
  units_all <- list(); summ_all <- list(); par_all <- list()
  
  for (k in seq_len(nrow(models))) {
    rq <- models$rq[k]; folder <- models$folder[k]
    dir_m <- model_path(folder, must_exist = FALSE)
    dir_b <- model_path(paste0(folder, "-broad"), must_exist = FALSE)
    if (!dir.exists(dir_m) || !dir.exists(dir_b)) {
      msg("Prior sensitivity: skipped ", folder, " (main or broad folder missing)")
      next
    }
    diag_b <- fread(file.path(dir_b, "diagnostics.csv"))
    files  <- list.files(dir_m, pattern = "^fit-[^-]+-(full|confirmatory)\\.rds$")
    
    for (f in files) {
      tr     <- sub("^fit-([^-]+)-(full|confirmatory)\\.rds$", "\\1", f)
      smp    <- sub("^fit-([^-]+)-(full|confirmatory)\\.rds$", "\\2", f)
      suffix <- paste0("-", tr, "-", smp, ".rds")
      if (!file.exists(file.path(dir_b, f))) {
        msg("Prior sensitivity: no broad fit for ", folder, " ", tr, "/", smp); next
      }
      
      post_m <- rstan::extract(readRDS(file.path(dir_m, f)))
      post_b <- rstan::extract(readRDS(file.path(dir_b, f)))
      lev <- function(dir, what) as.character(readRDS(file.path(dir, paste0(what, suffix))))
      stopifnot(identical(lev(dir_m, "pid-levels"), lev(dir_b, "pid-levels")),
                identical(lev(dir_m, "seq-levels"), lev(dir_b, "seq-levels")))
      
      u <- rbind(
        compare_level(rq, "sequence", lev(dir_m, "seq-levels"),
                      post_m[[est[[rq]][["s"]]]], post_b[[est[[rq]][["s"]]]], post_m, post_b),
        compare_level(rq, "participant", lev(dir_m, "pid-levels"),
                      post_m[[est[[rq]][["i"]]]], post_b[[est[[rq]][["i"]]]], post_m, post_b))
      u[, `:=`(model = folder, rq = rq, treatment = tr, sample = smp)]
      units_all[[length(units_all) + 1L]] <- u
      
      # converged_broad: all sampler and data checks of the broad fit passed
      # (R-hat, ESS, divergences, treedepth, E-BFMI, sample size); posterior
      # predictive checks (check starting with "ppc") describe model fit, not
      # convergence, and are excluded
      # (tr_f and smp_f avoid the data.table columns of the same name)
      tr_f <- tr; smp_f <- smp
      conv <- diag_b[treatment == tr_f & sample == smp_f & !grepl("^ppc", check),
                     all(as.logical(pass))]
      s <- u[, .(n_units = .N,
                 max_abs_shift_sd = max(abs(shift_sd)),
                 n_position_changed = sum(pos_main != pos_broad),
                 n_label_changed = sum(label_main != label_broad),
                 max_step = max(step)), by = level]
      s[, `:=`(crit_shift = max_abs_shift_sd < thr_shift,
               crit_position = n_position_changed == 0L,
               crit_labels = max_step <= thr_step)]
      s[, `:=`(robust = crit_shift & crit_position & crit_labels,
               primary = level == "sequence", converged_broad = conv,
               model = folder, rq = rq, treatment = tr, sample = smp,
               selected = models$selected[k])]
      summ_all[[length(summ_all) + 1L]] <- s
      
      # scalar parameters under both priors
      scal <- function(p) names(p)[vapply(p, function(x) length(dim(x)) <= 1L, logical(1))]
      nms  <- setdiff(intersect(scal(post_m), scal(post_b)), "lp__")
      pr <- rbindlist(lapply(nms, function(nm) {
        xm <- as.numeric(post_m[[nm]]); xb <- as.numeric(post_b[[nm]])
        data.table(parameter = nm, mean_main = mean(xm), sd_main = sd(xm),
                   q025_main = unname(quantile(xm, 0.025)), q975_main = unname(quantile(xm, 0.975)),
                   mean_broad = mean(xb), sd_broad = sd(xb),
                   q025_broad = unname(quantile(xb, 0.025)), q975_broad = unname(quantile(xb, 0.975)),
                   shift_sd = (mean(xb) - mean(xm)) / sd(xm))
      }))
      pr[, `:=`(model = folder, rq = rq, treatment = tr, sample = smp)]
      par_all[[length(par_all) + 1L]] <- pr
      
      rm(post_m, post_b); invisible(gc())
    }
  }
  
  if (length(summ_all) == 0L) stop("Prior sensitivity: no model pair with a broad folder found.", call. = FALSE)
  
  out <- list(`prior-sensitivity-summary.csv` = rbindlist(summ_all),
              `prior-sensitivity-units.csv` = rbindlist(units_all),
              `prior-sensitivity-parameters.csv` = rbindlist(par_all))
  for (nm in names(out)) {
    fwrite(out[[nm]], file.path(path_csv, nm))
    msg("Saved: ", file.path(path_csv, nm))
  }
  invisible(TRUE)
}
