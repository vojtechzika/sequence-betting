# ============================================================
# 53_ex1-regression-tables.R
#
# PURPOSE
#   Summarises the preregistered EX1 analyses that follow the index.
#   EX1-II: posterior summary of the regression of the participant index chi_i
#   on standardised optimism, log response time and risk aversion (pooled
#   over the fits of 54_ex1-stan-regression.R): median, 95% credible
#   interval and P(parameter > 0).
#   EX1-III: association between the pure-sequence score chi_pure
#   (52_ex1-pure-sequences.R) and the full index, as the posterior
#   distribution of cor(chi_pure, chi_i) across the draws of chi_i, with
#   median, 95% interval and P(cor > 0), and the posterior median of chi_i by
#   pure-sequence type.
#
# INPUT
#   model "ex1-regression": pooled-<tr>.rds, diagnostics.csv
#   model "ex1-chi":        chi-draws-<tr>-eps<main>.rds
#   data/output/csv/ex1-<tr>-pure-participants.csv   (52_ex1-pure-sequences.R)
#
# OUTPUT (data/output/csv/)
#   ex1-<tr>-regression-summary.csv   one row per parameter
#   ex1-<tr>-pure-vs-index.csv        EX1-III correlation and chi_i by type
#
# NOTES
#   - Draws of chi_i are those at the main tolerance (first value of
#     cfg$design$ex1$eps). If chi_i has no variation across participants in
#     a draw the correlation is undefined for that draw; the number of
#     defined draws is reported
#   - Pearson correlation, over participants with a classified pure-sequence
#     score
#   - Output files are overwritten at every run
# ============================================================

ex1_regression_tables <- function(cfg) {

  tr_vec <- unique(as.character(cfg$run$treatment))
  eps_main <- as.numeric(cfg$design$ex1$eps)[1]
  eps_nm   <- gsub("\\.", "", sprintf("%.2f", eps_main))

  summ <- function(x) c(median = median(x), q025 = unname(quantile(x, 0.025)),
                        q975 = unname(quantile(x, 0.975)), p_pos = mean(x > 0))

  for (tr in tr_vec) {

    # ---- EX1-II ----
    reg <- readRDS(model_path("ex1-regression", paste0("pooled-", tr, ".rds")))
    pm  <- as.matrix(reg$pooled)
    tab <- rbindlist(lapply(colnames(pm), function(p) {
      as.list(c(parameter = p, summ(pm[, p])))
    }))
    tab[, c("median", "q025", "q975", "p_pos") := lapply(.SD, as.numeric),
        .SDcols = c("median", "q025", "q975", "p_pos")]
    tab[, `:=`(treatment = tr, n_participants = reg$N, n_fits = reg$Trep)]
    f1 <- file.path(path_csv, paste0("ex1-", tr, "-regression-summary.csv"))
    fwrite(tab, f1);  msg("Saved: ", f1)

    # ---- EX1-III ----
    chi  <- readRDS(model_path("ex1-chi", paste0("chi-draws-", tr, "-eps", eps_nm, ".rds")))
    pure <- fread(file.path(path_csv, paste0("ex1-", tr, "-pure-participants.csv")),
                  encoding = "UTF-8", colClasses = c(pid = "character"))
    pure <- pure[match(as.character(chi$pid_levels), pid)]
    ok   <- !is.na(pure$chi_pure)

    cors <- suppressWarnings(
      as.numeric(cor(t(chi$chi_i[, ok, drop = FALSE]), pure$chi_pure[ok])))
    cors <- cors[is.finite(cors)]
    row_cor <- if (length(cors)) as.list(summ(cors)) else
      list(median = NA_real_, q025 = NA_real_, q975 = NA_real_, p_pos = NA_real_)
    cor_tbl <- data.table(treatment = tr, quantity = "cor(chi_pure, chi_i)",
                          group = "all classified", n = sum(ok), n_draws_defined = length(cors))
    cor_tbl <- cbind(cor_tbl, as.data.table(row_cor))

    by_type <- rbindlist(lapply(c(1L, 0L, -1L), function(v) {
      idx <- which(pure$chi_pure == v)
      if (!length(idx)) return(NULL)
      draws <- rowMeans(chi$chi_i[, idx, drop = FALSE])
      cbind(data.table(treatment = tr, quantity = "mean chi_i",
                       group = c(`1` = "hothand_consistent", `0` = "inconsistent",
                                 `-1` = "gambler_consistent")[[as.character(v)]],
                       n = length(idx), n_draws_defined = length(draws)),
            as.data.table(as.list(summ(draws))))
    }))
    f2 <- file.path(path_csv, paste0("ex1-", tr, "-pure-vs-index.csv"))
    fwrite(rbind(cor_tbl, by_type), f2);  msg("Saved: ", f2)
  }

  invisible(TRUE)
}
