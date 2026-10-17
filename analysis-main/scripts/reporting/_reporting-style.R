# ============================================================
# _reporting-style.R
#
# PURPOSE
#   Defines how the treatments are labelled and coloured in tables and
#   figures of the reporting stage, so that all scripts use the same
#   labels and colours. The multipliers are read from the design
#   (cfg$design$seq$treatments); the treatment type (FN: false negative,
#   FP: false positive) and the colours are presentation choices made here.
#
# FUNCTION
#   treatment_style(cfg)   data.table with one row per treatment in
#                          cfg$run$treatment: treat (code), type, label
#                          (e.g. "FN (m = 2.5)"), colour, colour_light
#
# NOTES
#   - colour_light is used where a second shade is needed (for example men
#     in the figures by sex)
#   - Stops if a treatment of the run has no type or colour defined here
# ============================================================

treatment_style <- function(cfg) {
  
  tr   <- unique(as.character(cfg$run$treatment))
  mult <- unlist(cfg$design$seq$treatments)
  
  type         <- c(m25 = "FN",      m19 = "FP")
  colour       <- c(m25 = "#2166AC", m19 = "#C0392B")
  colour_light <- c(m25 = "#74ADD1", m19 = "#E8927C")
  
  stopifnot(all(tr %in% names(mult)), all(tr %in% names(type)))
  
  data.table(
    treat        = tr,
    type         = unname(type[tr]),
    label        = sprintf("%s (m = %s)", type[tr], format(mult[tr])),
    colour       = unname(colour[tr]),
    colour_light = unname(colour_light[tr])
  )
}
