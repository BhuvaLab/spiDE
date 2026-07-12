# Shared, colour-blind-friendly palette + ggplot2 theme for the spiDE vignettes.
# Sourced from each vignette via
#   source(system.file("spide_palette.R", package = "spiDE"))
# so a given series (method, combiner, cell type, layout) always maps to the
# same colour across every figure. Colours are drawn from ColorBrewer's
# qualitative Dark2 / Set2 palettes (both colour-blind safe).

spide_pal <- list(
  # mixed-effects modes
  method = c(fixed = "#D95F02", intercept = "#1B9E77", slope = "#7570B3"),
  # within-gene p-value combiners
  combine = c(Brown = "#7570B3", Cauchy = "#D95F02"),
  # cell types (Set2)
  celltype = c(A = "#66C2A5", B = "#FC8D62", C = "#8DA0CB", D = "#E78AC3"),
  # spatial niche layouts (Dark2)
  layout = c(gradient = "#1B9E77", clustered = "#D95F02",
             random = "#7570B3", multiniche = "#E7298A"),
  accent = "#1B9E77",   # single-series highlight
  ref = "grey55"        # reference lines (diagonal / nominal level)
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))
  ggplot2::theme_update(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 12)
  )
}

# convenience scale constructors ------------------------------------------------
scale_colour_spide <- function(which, name = NULL, ...) {
  ggplot2::scale_colour_manual(values = spide_pal[[which]], name = name, ...)
}
scale_fill_spide <- function(which, name = NULL, ...) {
  ggplot2::scale_fill_manual(values = spide_pal[[which]], name = name, ...)
}
