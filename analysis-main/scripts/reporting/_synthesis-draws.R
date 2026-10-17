# ============================================================
# _synthesis-draws.R
#
# PURPOSE
#   Shared helpers of the cost-versus-direction scripts of the Synthesis
#   (98a_synthesis-costdir-sequences.R, 98b_synthesis-costdir-variance.R,
#   98c_synthesis-costdir-participants.R). The scripts combine posterior
#   draws of the selected RQ1-RQ4 models and of the joint synthesis model.
#   The helpers draw K posterior draws from a fit, summarise a vector of
#   draws, and write LaTeX tables.
#
# FUNCTION
#   synthesis_draw_tools(cfg)   returns a list with
#     K, min_bets               number of draws and bets threshold
#                               (cfg$synthesis$costdir)
#     load_rq(key, tr, sample, mu_s, mu_i)
#                               K draws of the sequence means (K x S) and the
#                               participant means (K x N) of a selected model
#     load_rq4(tr)              as load_rq for RQ4 (full sample), with the
#                               baseline hbar
#     load_syn(tr)              K draws of the personal side-choice slopes of
#                               the joint model (run and imbalance)
#     sumrow(v, ...)            median, 95% interval and P(v > 0) of draws v
#     fmt, cell, write_tex      LaTeX table helpers
#     wspear, wmean, wmed, dirichlet
#                               weighted statistics and Bayesian-bootstrap
#                               weights
#
# NOTES
#   - Each model is estimated separately, so draws are subsampled
#     independently from each fit and paired by row index; within a model the
#     draws stay joint across sequences and participants
#   - Sampling of draws uses the current random seed; the calling script sets
#     it from cfg$run$seed
# ============================================================

synthesis_draw_tools <- function(cfg) {
  
  K        <- as.integer(cfg$synthesis$costdir$n_draw)
  min_bets <- as.integer(cfg$synthesis$costdir$min_bets)
  
  pick  <- function(n) sample.int(n, K, replace = n < K)
  # rstan returns scalar parameters as one-dimensional arrays
  rows  <- function(x, idx) {
    if (length(dim(x)) <= 1L) as.vector(x)[idx]
    else if (length(dim(x)) == 2L) x[idx, , drop = FALSE]
    else x[idx, , , drop = FALSE]
  }
  named <- function(m, nm) { colnames(m) <- nm; m }
  
  qs <- function(x) c(median = median(x, na.rm = TRUE),
                      q025   = unname(quantile(x, .025, na.rm = TRUE)),
                      q975   = unname(quantile(x, .975, na.rm = TRUE)),
                      p_gt0  = mean(x > 0, na.rm = TRUE))
  sumrow <- function(v, ...) data.table(..., t(qs(v)))
  
  # ---- model loaders ----
  # Fit and levels of the model selected under `key`; NULL if absent
  load_fit <- function(key, tr, sample, pars) {
    f <- function(stem)
      tryCatch(model_path(key, paste0(stem, "-", tr, "-", sample, ".rds")),
               error = function(e) NA_character_)
    paths <- c(f("fit"), f("pid-levels"), f("seq-levels"))
    if (anyNA(paths) || !all(file.exists(paths))) return(NULL)
    list(post = rstan::extract(readRDS(paths[1L]), pars = pars),
         pid  = as.character(readRDS(paths[2L])),
         seq  = as.character(readRDS(paths[3L])),
         folder = basename(model_path(key)))
  }
  load_rq <- function(key, tr, sample, mu_s, mu_i) {
    L <- load_fit(key, tr, sample, c(mu_s, mu_i))
    if (is.null(L)) return(NULL)
    idx <- pick(nrow(L$post[[mu_s]]))
    list(seq = named(rows(L$post[[mu_s]], idx), L$seq),
         pid = named(rows(L$post[[mu_i]], idx), L$pid),
         folder = L$folder)
  }
  load_rq4 <- function(tr) {
    L <- load_fit("rq4", tr, "full", c("mu_h", "hbar", "mu_h_i"))
    if (is.null(L)) return(NULL)
    idx <- pick(length(L$post$hbar))
    list(seq  = named(rows(L$post$mu_h, idx), L$seq),
         hbar = rows(L$post$hbar, idx),
         pid  = named(rows(L$post$mu_h_i, idx), L$pid),
         folder = L$folder)
  }
  # Personal side-choice slopes of the joint model: fixed effect plus the
  # participant deviation (side run: beta_side[2] + u[5], side imbalance:
  # beta_side[3] + u[6])
  load_syn <- function(tr) {
    L <- load_fit("synthesis", tr, "full", c("beta_side", "u"))
    if (is.null(L)) return(NULL)
    idx <- pick(nrow(L$post$beta_side))
    u   <- rows(L$post$u, idx); bs <- rows(L$post$beta_side, idx)
    list(run = named(bs[, 2] + u[, 5, ], L$pid),
         imb = named(bs[, 3] + u[, 6, ], L$pid))
  }
  
  # ---- weighted statistics ----
  wcor <- function(x, y, w) {
    mx <- sum(w * x); my <- sum(w * y)
    sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
  }
  wspear    <- function(x, y, w) wcor(rank(x), rank(y), w)
  dirichlet <- function(n) { g <- rgamma(n, 1); g / sum(g) }
  wmean     <- function(x, w) sum(w * x) / sum(w)
  wmed      <- function(x, w) {
    o <- order(x); cw <- cumsum(w[o]) / sum(w); unname(x[o][which(cw >= 0.5)[1]])
  }
  
  # ---- LaTeX ----
  fmt  <- function(x) { x <- ifelse(abs(x) < 0.005, 0, x); sprintf("%.2f", x) }
  cell <- function(r, digits_fmt = fmt, stars = TRUE) {
    star <- if (stars && (r$q025 > 0 || r$q975 < 0)) "^{*}" else ""
    c(sprintf("$%s%s$", digits_fmt(r$median), star),
      sprintf("$[%s, %s]$", digits_fmt(r$q025), digits_fmt(r$q975)))
  }
  write_tex <- function(f, caption, label, colspec, header, body, footer, note) {
    writeLines(c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", caption, "}"),
                 paste0("\\label{", label, "}"), "\\begin{threeparttable}",
                 paste0("\\begin{tabular}{", colspec, "}"), "\\toprule", header, "\\midrule",
                 body, "\\midrule", footer, "\\bottomrule", "\\end{tabular}",
                 "\\begin{tablenotes}[flushleft]", "\\footnotesize",
                 paste0("\\item \\textit{Notes:} ", note), "\\end{tablenotes}",
                 "\\end{threeparttable}", "\\end{table}"), f)
    msg("Saved: ", f)
  }
  
  list(K = K, min_bets = min_bets, pick = pick, rows = rows, named = named,
       sumrow = sumrow, load_rq = load_rq, load_rq4 = load_rq4, load_syn = load_syn,
       wspear = wspear, wmean = wmean, wmed = wmed, dirichlet = dirichlet,
       fmt = fmt, cell = cell, write_tex = write_tex)
}
