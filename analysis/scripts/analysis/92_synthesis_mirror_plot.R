# ============================================================
# 92_synthesis_mirror_pairs.R
#
# PURPOSE
#   Scatter plot of ALL 64 sequences in side-choice x betting-rate
#   space: x = |heads-share diff from grand mean| (abs(h_d_median)),
#   y = betting-rate diff from grand mean (b_d_median, signed).
#
#   Every sequence is plotted and every mirror pair gets a connecting
#   line (see mirror_seq() below). Visual weight is reserved for pairs
#   with at least one CREDIBLE betting-rate shift (b_sig = TRUE, 95% CI
#   excludes 0): those points/lines use the normal colour and get a
#   sequence-name label. Pairs where NEITHER member is credible are
#   still drawn -- point and connecting line in grey, faded, thinner --
#   but unlabeled, so the full distribution stays visible without
#   turning the plot into 64 overlapping text labels. A pair with
#   exactly one credible member shows both endpoints labeled and
#   coloured individually (the credible one in normal colour, its
#   non-credible partner in grey), with a mid-weight connecting line --
#   distinct from both the "both credible" and "neither credible" tiers
#   so the eye isn't drawn to it as strongly as a fully-credible pair,
#   but it doesn't recede into the background either.
#
#   x is an ABSOLUTE VALUE, not the signed delta: within a mirror pair,
#   heads-share deltas are near-antisymmetric (flipping every H<->T just
#   relabels which side is "heads"), so the signed value would put
#   mirror partners on opposite sides of x=0 for no informative reason.
#   Taking |h_d_median| puts both members of a pair at (roughly) the
#   same x, so the connecting line's vertical extent directly shows
#   whether the pair's betting rate moves the same way (short line) or
#   oppositely (long line) despite an equally strong side pull. The axis
#   is labelled side-bias MAGNITUDE, not "heads share," since the sign
#   is discarded by construction.
#
#   A shaded region (very light, mostly transparent) marks where side
#   bias itself is not credible: the boundary is the empirical midpoint
#   between the largest abs(h_d_median) among h_sig = FALSE sequences
#   and the smallest among h_sig = TRUE sequences. This is only drawn as
#   a single line because abs(h_d_median) happens to separate h_sig
#   TRUE/FALSE cleanly in this data (checked at runtime, see stop()
#   below) -- it is an empirical property of this particular posterior,
#   not a fixed a priori cutoff, and would need re-checking after any
#   refit.
#
#   Two aesthetics carry three pieces of information: point/line COLOUR
#   marks b_sig (normal colour = TRUE, grey = FALSE); point SHAPE is a
#   combined direction x side-bias-credibility code -- filled circle
#   (Head-biased, h_sig TRUE), hollow circle (Head-biased, h_sig FALSE),
#   filled diamond (Tail-biased, h_sig TRUE), hollow diamond (Tail-biased,
#   h_sig FALSE). Plain posterior-median points only, no CI whiskers.
#   Point alpha and connecting-line weight/alpha are a further, purely
#   presentational channel (not a new legend entry) used only to mute the
#   unlabeled
#   "neither member credible" pairs so they read as background context
#   rather than competing with the labeled findings.
#
# INPUT
#   path_out/synthesis_sequence_diffs_<tr>_<tag>.csv, produced by
#   synthesis_sequences(cfg) in 93_synthesis_sequences.R. Must be the
#   ALL-SEQUENCES run (default, unrestricted `sequences` argument) -- a
#   restricted run missing one member of a pair still plots the
#   surviving member with no line, and warns (unless `quiet = TRUE`).
#
# OUTPUT
#   <path_fig or path_out>/betting_side_scatter_mirror_<tr>_<tag>.png
#
# NOTES
#   - `quiet = TRUE` (default) suppresses the one informational warning
#     this file can raise (a mirror sequence missing from the input
#     data entirely, excluding its pair from the connecting-line layer).
#     Genuine errors (missing input file/columns, empty input, or
#     abs(h_d_median) failing to separate h_sig TRUE/FALSE) always
#     stop(). should_skip() gates the PNG write; the plot object is
#     still built and returned either way.
#   - Return value: invisible(list(plot=, data=, pairs=, path=)).
# ============================================================
synthesis_mirror_plot <- function(cfg, tr = "m25", tag = "confirmatory",
                                  out_file = NULL, quiet = TRUE) {
  
  # ---- Helper: H<->T mirror of a sequence string, e.g. "HHHTTT" ->
  # "TTTHHH". For 6-character H/T sequences this can never map a
  # sequence to itself, so all 64 sequences partition into exactly 32
  # disjoint mirror pairs ----
  mirror_seq <- function(s) chartr("HT", "TH", s)
  
  # ---- Helper: build one row per mirror PAIR (not per sequence) with
  # both sequences' plotting coordinates, so geom_segment can draw a
  # single connecting line per pair without doubling up (A->B and B->A) ----
  build_mirror_pairs <- function(dt, quiet) {
    dt <- copy(dt)
    dt[, mirror := mirror_seq(sequence)]
    
    missing_mirror <- setdiff(dt$mirror, dt$sequence)
    if (length(missing_mirror) > 0L && !quiet) {
      warning("synthesis_mirror_plot: mirror sequence(s) not present in the ",
              "input data, so their pair(s) are dropped from the mirror-line ",
              "layer (points for the sequences that ARE present still plot): ",
              paste(unique(missing_mirror), collapse = ", "))
    }
    dt <- dt[mirror %in% sequence]
    
    # Canonical, order-independent key so each pair is processed once
    # regardless of which member appears first in the input.
    dt[, pair_key := paste(pmin(sequence, mirror), pmax(sequence, mirror))]
    
    pair_keys <- unique(dt$pair_key)
    rbindlist(lapply(pair_keys, function(pk) {
      two <- dt[pair_key == pk]
      if (nrow(two) != 2L) return(NULL)  # defensive: skip anything incomplete
      data.table(
        seq1 = two$sequence[1], seq2 = two$sequence[2],
        b1 = two$b_d_median[1], h1 = two$h_d_median[1],
        b2 = two$b_d_median[2], h2 = two$h_d_median[2]
      )
    }))
  }
  
  # ---- Main body ----
  f_in <- file.path(path_out, paste0("synthesis_sequence_diffs_", tr, "_", tag, ".csv"))
  if (!file.exists(f_in)) {
    stop("synthesis_mirror_plot: expected input not found: ", f_in,
         " -- run synthesis_sequences(cfg) first (93_synthesis_sequences.R).")
  }
  dt <- fread(f_in)
  
  required_cols <- c("sequence", "h_d_median", "b_d_median", "b_sig", "h_sig")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("synthesis_mirror_plot: input CSV is missing column(s): ",
         paste(missing_cols, collapse = ", "), ". Was it produced by the ",
         "current version of 93_synthesis_sequences.R?")
  }
  
  # Show every sequence -- credibility now controls colour/label/line
  # weight, not inclusion.
  dt[, mirror := mirror_seq(sequence)]
  sig_seqs <- dt[b_sig == TRUE, sequence]
  dt_show <- copy(dt)
  if (nrow(dt_show) == 0L) {
    stop("synthesis_mirror_plot: input has no rows -- nothing to plot.")
  }
  
  dt_show[, credible   := b_sig]
  dt_show[, show_label := credible | mirror %in% sig_seqs]  # pair has >=1 credible member
  dt_show[, pt_alpha   := fifelse(show_label, 1, 0.55)]
  dt_show[, direction  := fifelse(h_d_median > 0, "Head-biased", "Tail-biased")]
  dt_show[, direction  := factor(direction, levels = c("Head-biased", "Tail-biased"))]
  
  # Combined shape: direction (circle vs. diamond) x side-bias credibility
  # (filled vs. hollow), so both are read off the symbol itself rather than
  # a separate fill channel -- filled circle / hollow circle / filled
  # diamond / hollow diamond.
  dt_show[, shape_cat := fifelse(
    direction == "Head-biased" & h_sig,  "Head-biased, side bias credible",
    fifelse(direction == "Head-biased" & !h_sig, "Head-biased, side bias not credible",
            fifelse(direction == "Tail-biased" & h_sig,  "Tail-biased, side bias credible",
                    "Tail-biased, side bias not credible")))]
  dt_show[, shape_cat := factor(shape_cat, levels = c(
    "Head-biased, side bias credible", "Head-biased, side bias not credible",
    "Tail-biased, side bias credible", "Tail-biased, side bias not credible"
  ))]
  
  # ggplot draws solid pch 15-20 and open pch 0-14 at inconsistent visual
  # sizes for the same `size` value -- solid glyphs render smaller,
  # diamonds (18 vs. 5) especially so. These per-shape multipliers are an
  # empirical correction so filled and hollow symbols read as the same
  # size on the page; re-check visually after rendering and adjust if
  # your device/DPI shows a different ratio.
  size_mult <- c(
    `Head-biased, side bias credible`     = 1.25,  # pch 16, filled circle
    `Head-biased, side bias not credible` = 1.00,  # pch 1,  hollow circle
    `Tail-biased, side bias credible`     = 1.55,  # pch 18, filled diamond
    `Tail-biased, side bias not credible` = 1.00   # pch 5,  hollow diamond
  )
  dt_show[, pt_size := 1.6 * size_mult[as.character(shape_cat)]]
  
  # Empirical side-bias significance threshold. Verify abs(h_d_median)
  # actually separates h_sig TRUE/FALSE before trusting a single cutoff
  # line -- if it doesn't, a line would misclassify some sequences and
  # this stops rather than drawing something misleading.
  x_sig    <- abs(dt_show[h_sig == TRUE,  h_d_median])
  x_nonsig <- abs(dt_show[h_sig == FALSE, h_d_median])
  if (length(x_sig) == 0L || length(x_nonsig) == 0L || min(x_sig) <= max(x_nonsig)) {
    stop("synthesis_mirror_plot: abs(h_d_median) does not cleanly separate ",
         "h_sig TRUE/FALSE for tr=", tr, ", tag=", tag, " -- a single side-bias ",
         "threshold line would misclassify some sequences. Inspect h_sig vs. ",
         "abs(h_d_median) directly instead of using this plot's shaded region.")
  }
  h_thresh <- mean(c(max(x_nonsig), min(x_sig)))
  
  pairs <- build_mirror_pairs(dt_show, quiet = quiet)
  
  p <- ggplot(dt_show, aes(x = abs(h_d_median), y = b_d_median)) +
    # Fully transparent -- the side-bias n.s. region is now marked only by
    # the dashed vline + annotation below, no background tint.
    annotate("rect", xmin = -Inf, xmax = h_thresh, ymin = -Inf, ymax = Inf,
             fill = "grey50", alpha = 0, colour = NA) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey40") +
    geom_vline(xintercept = h_thresh, linetype = "dashed", colour = "grey50") +
    annotate("text", x = h_thresh, y = Inf, label = "side-bias n.s.  ←",
             hjust = 1.05, vjust = 1.5, size = 3, colour = "grey50")
  
  if (nrow(pairs) > 0L) {
    # Three visual tiers instead of a single credible/not-credible split:
    # pairs with NEITHER member credible are unlabeled background context,
    # so they get the thinnest, faintest line; pairs with exactly ONE
    # credible member sit in the middle; pairs with BOTH credible are the
    # normal, most prominent line. Drawn as separate literal-aesthetic
    # layers (not a single mapped `colour`/`linewidth` scale) so this
    # doesn't collide with the point layer's credible-based colour scale.
    pairs[, any_credible  := seq1 %in% dt_show[credible == TRUE, sequence] |
            seq2 %in% dt_show[credible == TRUE, sequence]]
    pairs[, both_credible := seq1 %in% dt_show[credible == TRUE, sequence] &
            seq2 %in% dt_show[credible == TRUE, sequence]]
    
    pairs_neither <- pairs[any_credible == FALSE]
    pairs_one     <- pairs[any_credible == TRUE & both_credible == FALSE]
    pairs_both    <- pairs[both_credible == TRUE]
    
    if (nrow(pairs_neither) > 0L) {
      p <- p + geom_segment(
        data = pairs_neither, inherit.aes = FALSE,
        aes(x = abs(h1), y = b1, xend = abs(h2), yend = b2),
        colour = "grey85", linewidth = 0.18, alpha = 0.5, linetype = "dashed"
      )
    }
    if (nrow(pairs_one) > 0L) {
      p <- p + geom_segment(
        data = pairs_one, inherit.aes = FALSE,
        aes(x = abs(h1), y = b1, xend = abs(h2), yend = b2),
        colour = "grey55", linewidth = 0.3, alpha = 0.8, linetype = "solid"
      )
    }
    if (nrow(pairs_both) > 0L) {
      p <- p + geom_segment(
        data = pairs_both, inherit.aes = FALSE,
        aes(x = abs(h1), y = b1, xend = abs(h2), yend = b2),
        colour = "steelblue4", linewidth = 0.45, alpha = 0.9, linetype = "solid"
      )
    }
  }
  
  p <- p +
    geom_point(aes(colour = credible, shape = shape_cat, alpha = I(pt_alpha),
                   size = I(pt_size)),
               stroke = 0.6)
  
  # Labels only for sequences in a pair with >=1 credible member -- this
  # is what keeps a 64-sequence plot readable, not the point styling.
  dt_label <- dt_show[show_label == TRUE]
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = dt_label, aes(label = sequence, colour = credible), size = 2.8,
      max.overlaps = Inf, seed = 1, show.legend = FALSE
    )
  } else {
    p <- p + geom_text(
      data = dt_label, aes(label = sequence, colour = credible), size = 2.8,
      vjust = -0.8, hjust = 0.3, show.legend = FALSE
    )
  }
  
  p <- p +
    scale_colour_manual(
      values = c(`TRUE` = "steelblue4", `FALSE` = "grey65"),
      labels = c(`TRUE` = "Betting rate credible (95% CI excl. 0)", `FALSE` = "Not credible"),
      name = NULL
    ) +
    # Shape alone carries two pieces of information: circle vs. diamond for
    # direction (Head-biased vs. Tail-biased), filled vs. hollow for
    # side-bias credibility (h_sig TRUE vs. FALSE). No separate fill aes.
    scale_shape_manual(
      values = c(
        `Head-biased, side bias credible`     = 16,  # filled circle
        `Head-biased, side bias not credible` = 1,   # hollow circle
        `Tail-biased, side bias credible`     = 18,  # filled diamond
        `Tail-biased, side bias not credible` = 5    # hollow diamond
      ),
      name = NULL
    ) +
    labs(
      x = "\n Side-bias Δ",
      y = "Betting-rate Δ \n"
    ) +
    # theme_classic(): white background, x/y axis lines, no gridlines at
    # all -- avoids the theme_minimal() grid-vs-LaTeX weight issue outright
    # rather than just thinning it.
    theme_classic(base_size = 10) +
    theme(legend.position = "none", legend.box = "vertical")
  
  fig_dir <- if (exists("path_fig")) path_fig else path_out
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(out_file)) {
    out_file <- file.path(fig_dir, paste0("betting_side_scatter_mirror_", tr, "_", tag, ".png"))
  }
  
  if (!should_skip(out_file, cfg, "output", paste0("Synthesis mirror-pair scatter (", tr, "/", tag, ")"))) {
    # ggsave itself stays silent here because width/height are supplied
    # explicitly -- it only prints a "Saving NxM in image" message when
    # left to guess the size from the current device.
    ggsave(out_file, p, width = 10, height = 5, dpi = 300)
    msg("Saved: ", out_file)
  }
  
  invisible(list(plot = p, data = dt_show, pairs = pairs, path = out_file))
}