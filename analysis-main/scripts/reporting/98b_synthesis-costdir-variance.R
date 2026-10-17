# ============================================================
# 98b_synthesis-costdir-variance.R
#
# PURPOSE
#   Synthesis, who versus what: the posterior of the participant-level and
#   the sequence-level random-intercept standard deviations of the models
#   for betting (RQ1), welfare loss (RQ3) and side choice (RQ4), their
#   ratio and the participant share sigma_u^2 / (sigma_u^2 + sigma_s^2).
#   The scales differ across models (logit, log), so only the within-model
#   ratio and share are comparable across outcomes. The table tests the
#   closing claim that the participant drives cost whereas the sequence
#   drives direction. Cells report posterior medians and 95% credible
#   intervals.
#
# SAMPLE
#   Treatments with a confirmatory sample (cfg$design$a_flags$betting_normative):
#   RQ1 and RQ3 from the confirmatory fits; other treatments: full fits. RQ4
#   from the full fit.
#
# INPUT (models selected under the keys rq1, rq3, rq4)
#   fit-<tr>-<sample>.rds. Standard deviations by model folder:
#     rq1-bernoulli      sigma_u, sigma_s
#     rq3-hurdle-gamma   su0, sb0 (hurdle) and sup, sbp (positive loss)
#     rq3-gamma          sup, sbp
#     rq3-gaussian       sigma_u, sigma_s
#     rq4-bernoulli      sigma_u, sigma_s
#
# OUTPUT
#   data/output/csv/synthesis-<tr>-costdir-variance.csv
#   data/output/tables/synthesis-<tr>-costdir-variance.tex
#
# NOTES
#   - RQ2 is not included. The stake model has several participant and
#     sequence components (stake level and all-in propensity), which have no
#     common scale with the other outcomes
#   - The K draws (cfg$synthesis$costdir) are subsampled independently from
#     each fit; each ratio and share uses the paired draws of one model
#   - Seeded by cfg$run$seed
#   - Output files are overwritten at every run
# ============================================================

synthesis_costdir_variance <- function(cfg) {
  
  tr_vec <- unique(as.character(cfg$run$treatment))
  set.seed(as.integer(cfg$run$seed))
  tl <- synthesis_draw_tools(cfg)
  K  <- tl$K
  
  # participant and sequence standard deviations of each selected model
  vc_spec <- function(folder) switch(
    folder,
    "rq1-bernoulli"    = list(c("Betting", "logit", "sigma_u", "sigma_s")),
    "rq3-hurdle-gamma" = list(c("Welfare loss", "zero-loss hurdle", "su0", "sb0"),
                              c("Welfare loss", "positive-loss mean", "sup", "sbp")),
    "rq3-gamma"        = list(c("Welfare loss", "positive-loss mean", "sup", "sbp")),
    "rq3-gaussian"     = list(c("Welfare loss", "Gaussian", "sigma_u", "sigma_s")),
    "rq4-bernoulli"    = list(c("Side choice", "logit", "sigma_u", "sigma_s")),
    stop("No variance components defined for model folder '", folder, "'.", call. = FALSE))
  
  out <- list()
  
  for (tr in tr_vec) {
    
    conf_ok <- isTRUE(cfg$design$a_flags$betting_normative[[tr]])
    tg_main <- if (conf_ok) "confirmatory" else "full"
    
    vc_rows <- list()
    for (rq in list(c("rq1", tg_main), c("rq3", tg_main), c("rq4", "full"))) {
      folder <- basename(model_path(rq[1L]))
      for (sp in vc_spec(folder)) {
        L <- tryCatch({
          f <- model_path(rq[1L], paste0("fit-", tr, "-", rq[2L], ".rds"))
          rstan::extract(readRDS(f), pars = sp[3:4])
        }, error = function(e) NULL)
        if (is.null(L)) {
          warning("Synthesis variance components: ", rq[1L], " fit missing for ", tr, "/", rq[2L],
                  "; skipped.", call. = FALSE)
          next
        }
        su <- as.vector(L[[sp[3L]]]); ss <- as.vector(L[[sp[4L]]])
        idx <- tl$pick(length(su)); su <- su[idx]; ss <- ss[idx]
        d <- list(sd_participant = su, sd_sequence = ss, ratio = su / ss,
                  share_participant = su^2 / (su^2 + ss^2))
        vc_rows[[length(vc_rows) + 1L]] <- rbindlist(lapply(names(d), function(q)
          tl$sumrow(d[[q]], treatment = tr, outcome = sp[1L], component = sp[2L],
                    model = folder, sample = rq[2L], quantity = q)))
      }
    }
    vc <- rbindlist(vc_rows)
    if (!nrow(vc)) next
    vc[, n_draws := K]
    f1 <- file.path(path_csv, paste0("synthesis-", tr, "-costdir-variance.csv"))
    fwrite(vc, f1);  msg("Saved: ", f1)
    
    rowkeys <- unique(vc[, .(outcome, component, sample)])
    body <- unlist(lapply(seq_len(nrow(rowkeys)), function(i) {
      k  <- rowkeys[i]
      cs <- lapply(c("sd_participant", "sd_sequence", "ratio", "share_participant"), function(q)
        tl$cell(vc[outcome == k$outcome & component == k$component & quantity == q], stars = FALSE))
      lab <- paste0(k$outcome, ", ", k$component)
      c(paste(lab, "&", paste(sapply(cs, `[`, 1), collapse = " & "), "\\\\"),
        paste(" &", paste(sapply(cs, `[`, 2), collapse = " & "),
              if (i < nrow(rowkeys)) "\\\\[3pt]" else "\\\\"))
    }))
    tl$write_tex(
      file.path(path_tab, paste0("synthesis-", tr, "-costdir-variance.tex")),
      caption = "Participant versus sequence heterogeneity, by outcome",
      label   = paste0("Tab:synthesis-varcomp", if (tr != tr_vec[1L]) paste0("-", tr)),
      colspec = "l*{4}{c}",
      header  = " & $\\sigma_{\\text{participant}}$ & $\\sigma_{\\text{sequence}}$ & Ratio & Participant share \\\\",
      body    = body,
      footer  = paste("Posterior draws & \\multicolumn{4}{c}{", K, "} \\\\"),
      note    = paste0("Standard deviations of the participant- and sequence-level random intercepts in the ",
                       "selected model for each outcome (RQ1 betting, RQ3 welfare loss, RQ4 side choice), ",
                       "on each model's linear-predictor scale. Because scales differ across models, only the ",
                       "within-model ratio $\\sigma_{\\text{participant}}/\\sigma_{\\text{sequence}}$ and the ",
                       "participant share $\\sigma^2_{\\text{participant}}/(\\sigma^2_{\\text{participant}}+",
                       "\\sigma^2_{\\text{sequence}})$ are comparable across rows. Cells report posterior medians ",
                       "and 95\\% credible intervals. Side choice uses the full sample; the other outcomes the ",
                       if (conf_ok) "confirmatory" else "full", " sample."))
    
    out[[tr]] <- vc
  }
  
  invisible(out)
}
