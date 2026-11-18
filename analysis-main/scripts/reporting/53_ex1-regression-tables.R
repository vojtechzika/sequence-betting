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
#   EX1-I: classification of the sequences and of the participants by the
#   posterior probability that the directional index lies beyond the
#   threshold delta: H(.) = P(chi_s > delta), T(.) = P(chi_s < -delta) for
#   sequences and HH(.) = P(chi_i > delta), G(.) = P(chi_i < -delta) for
#   participants, with the preregistered labels (sequences: strong,
#   moderate, weak at 0.95, 0.80, 0.50; participants: solid, likely,
#   leaning at 0.95, 0.90, 0.75; neutral otherwise).
#
# INPUT
#   model "ex1-regression": pooled-<tr>.rds, diagnostics.csv
#   model "ex1-chi":        chi-draws-<tr>-eps<main>.rds
#   data/output/csv/ex1-<tr>-pure-participants.csv   (52_ex1-pure-sequences.R)
#   cfg$design$ex1$delta (first value = main threshold)
#
# OUTPUT (data/output/csv/)
#   ex1-<tr>-regression-summary.csv   one row per parameter
#   ex1-<tr>-pure-vs-index.csv        EX1-III correlation and chi_i by type
#   ex1-<tr>-sequences.csv            chi_s median and 95% interval, number of
#                                     bets, H(.), T(.), P(|chi_s| > delta) and
#                                     the class, for each delta
#   ex1-<tr>-participants.csv         the same for chi_i (HH(.), G(.))
#   ex1-<tr>-classification-summary.csv  number of sequences and participants
#                                     in each class, for each delta
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

    # ---- EX1-I: classification by the posterior probability beyond delta ----
    deltas <- as.numeric(cfg$design$ex1$delta)
    stopifnot(length(deltas) >= 1L, all(is.finite(deltas)), all(deltas > 0))
    dnm <- function(d) gsub("\\.", "", sprintf("%.2f", d))

    # Class from the probabilities of the positive and the negative direction
    class_dir <- function(p_pos, p_neg, cuts, nm, suf_pos, suf_neg) {
      lev <- function(p) ifelse(p >= cuts[1], 1L, ifelse(p >= cuts[2], 2L,
                          ifelse(p >= cuts[3], 3L, 0L)))
      lp <- lev(p_pos);  ln <- lev(p_neg)
      pos <- lp > 0L & (ln == 0L | p_pos >= p_neg)
      neg <- ln > 0L & !pos
      out <- rep("neutral", length(p_pos))
      out[pos] <- paste(nm[lp[pos]], suf_pos)
      out[neg] <- paste(nm[ln[neg]], suf_neg)
      out
    }

    unit_table <- function(M, ids, n, id_name, nm_pos, nm_neg, cuts, nm, suf_pos, suf_neg) {
      qs  <- apply(M, 2L, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
      tab <- data.table(id = ids, n_bets = as.integer(n),
                        chi_median = qs[2, ], chi_q025 = qs[1, ], chi_q975 = qs[3, ])
      setnames(tab, "id", id_name)
      for (d in deltas) {
        p_pos <- colMeans(M > d);  p_neg <- colMeans(M < -d)
        tab[, (paste0(nm_pos, "_delta_", dnm(d))) := p_pos]
        tab[, (paste0(nm_neg, "_delta_", dnm(d))) := p_neg]
        tab[, (paste0("P_abs_delta_", dnm(d))) := colMeans(abs(M) > d)]
        tab[, (paste0("class_delta_", dnm(d))) :=
              class_dir(p_pos, p_neg, cuts, nm, suf_pos, suf_neg)]
      }
      tab
    }

    seq_tab <- unit_table(chi$chi_s, as.character(chi$seq_levels), chi$n_s, "sequence",
                          "H", "T", c(0.95, 0.80, 0.50), c("strong", "moderate", "weak"),
                          "H-similar", "T-similar")
    pid_tab <- unit_table(chi$chi_i, as.character(chi$pid_levels), chi$n_i, "pid",
                          "HH", "G", c(0.95, 0.90, 0.75), c("solid", "likely", "leaning"),
                          "hot-handish", "gamblerish")
    seq_tab[, treatment := tr];  pid_tab[, treatment := tr]
    f3 <- file.path(path_csv, paste0("ex1-", tr, "-sequences.csv"))
    fwrite(seq_tab, f3);  msg("Saved: ", f3)
    f4 <- file.path(path_csv, paste0("ex1-", tr, "-participants.csv"))
    fwrite(pid_tab, f4);  msg("Saved: ", f4)

    cls_sum <- rbindlist(lapply(deltas, function(d) {
      cc <- paste0("class_delta_", dnm(d))
      rbind(seq_tab[, .(level = "sequence", delta = d, class = get(cc))][, .N, by = .(level, delta, class)],
            pid_tab[, .(level = "participant", delta = d, class = get(cc))][, .N, by = .(level, delta, class)])
    }))
    cls_sum[, treatment := tr]
    f5 <- file.path(path_csv, paste0("ex1-", tr, "-classification-summary.csv"))
    fwrite(cls_sum[order(level, delta, -N)], f5);  msg("Saved: ", f5)
  }

  invisible(TRUE)
}
