# Shared, colour-blind-friendly palette + ggplot2 theme for the spiDE vignettes.
# Sourced from each vignette via
#   source(system.file("spide_palette.R", package = "spiDE"))
# so a given series (method, combiner, cell type, layout) always maps to the
# same colour across every figure.
#
# EACH QUANTITY GETS ITS OWN HUE FAMILY. Two different quantities must not be
# drawn from the same set of hues, or the same colour would mean two different
# things across (or within) figures -- `method` and `category`, for instance,
# co-occur in the spatial-robustness breakdown. The families below are
# deliberately disjoint:
#
#   method   teal / orange / purple      (ColorBrewer Dark2)
#   category steel-blue / gold / crimson
#   layout   green / wine / navy / olive
#   celltype pink / brown / amber / blue
#
# Every family below was checked with the dataviz palette validator (six
# checks: lightness band, chroma floor, CVD separation, normal-vision floor,
# contrast) and PASSES on a light surface. Do NOT hand-edit a hex here without
# re-running that validator -- the previous Set2 `celltype` palette looked fine
# but failed hard: #E78AC3 vs #8DA0CB separated by only dE 1.5 under
# protanopia and dE 14.1 even in normal vision (below the 15 floor), and the
# old `category` green/brown pair failed CVD separation at dE 5.3 (deutan).
#
# Note the hue space for colour-blind-safe categorical palettes is finite:
# `condition` and `combine` are 2-level quantities that appear only in figures
# where no other categorical colour is present, so they reuse the Dark2 hues
# without ambiguity.

spide_pal <- list(
  # mixed-effects modes plus the two-stage estimator. The first three are
  # Dark2; the fourth is NOT Dark2's own pink #E7298A -- that hex fails the
  # dataviz validator against this family outright (worst all-pairs deutan
  # dE 1.7 against the intercept teal: a deuteranope cannot tell the twostage
  # line from the intercept line at all). #C51B7D is the best fourth hue found
  # under the validator with the trio held fixed: every check passes except
  # that its worst pair (vs #1B9E77, deutan) sits at dE 7.0, in the 6-8 floor
  # band that is legal only WITH a secondary encoding. The twostage series
  # therefore always carries a triangle point marker (`method_shape` below) in
  # line figures; grouped-bar figures get their secondary channel from the
  # dodge gaps. Without an entry here twostage plots as NA grey wherever
  # scale_colour_spide() is used, which is how it first appeared in the
  # simulation vignette.
  method = c(fixed = "#D95F02", intercept = "#1B9E77", slope = "#7570B3",
             twostage = "#C51B7D"),
  # point markers for the method family: the secondary identity channel the
  # floor-band twostage/intercept pair requires (see above). Map BOTH
  # colour and shape to `method` and the legends merge into one.
  method_shape = c(fixed = 16, intercept = 16, slope = 16, twostage = 17),
  # within-gene p-value combiners (own figures only; no clash in context)
  combine = c(Brown = "#7570B3", Cauchy = "#D95F02"),
  # one- vs two-sided input to a combiner. Brown keeps its `combine` colour;
  # the two Cauchy variants are split off it so the pair reads as one family.
  sided = c("Brown (1-sided)" = "#7570B3", "Cauchy (1-sided)" = "#D95F02",
            "Cauchy (2-sided)" = "#A02C3C"),
  # cell types -- pink / brown / amber / blue
  celltype = c(A = "#D8315B", B = "#FFD166", C = "#1E1B18", D = "#06D6A0"),
  # spatial niche layouts -- green / wine / navy / olive
  layout = c(gradient = "#117733", clustered = "#882255",
             random = "#2B6699", multiniche = "#8C6D1F"),
  # gene categories in the simulation -- steel-blue / gold / crimson
  category = c("cell-type marker" = "#3B6BA5", housekeeping = "#C77D00",
               background = "#A02C3C"),

  # planted differential-expression status
  de = c(DE = "#D95F02", "non-DE" = "grey70"),
  # experimental condition
  condition = c(Responder = "#1B9E77", "Non-responder" = "#7570B3"),
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

# convenience scale constructors -----------------------------------------------
scale_colour_spide <- function(which, name = NULL, ...) {
  ggplot2::scale_colour_manual(values = spide_pal[[which]], name = name, ...)
}
scale_fill_spide <- function(which, name = NULL, ...) {
  ggplot2::scale_fill_manual(values = spide_pal[[which]], name = name, ...)
}
scale_shape_spide <- function(which = "method_shape", name = NULL, ...) {
  ggplot2::scale_shape_manual(values = spide_pal[[which]], name = name, ...)
}
