# The spiDE plotting theme and colours (0.99.34; design/specs/2026-10-01-plots.md).
# One colour per role, the same in every plot, validated for colour blindness
# across roles (all pairs, the dataviz validator):
#   direction (up / down, and the diverging t scale): amber / violet
#   the two conditions (first level / second): teal / rufous, from the coolors
#     palette 001219-005f73-0a9396-94d2bd-e9d8a6-ee9b00-ca6702-bb3e03-ae2012-9b2226
#   magnitude (density, R^2, counts): scico "lapaz" (Crameri 2018), light to dark
#   a nuisance factor (slide, batch): Okabe-Ito
# The theme is vissE::bhuvad_theme() -- black panel border, no grid, boxed
# strips, italic legend titles -- with black axis ticks.

.spideCols <- c(up = "#B8791A", down = "#5B4B9A", first = "#0A9396", second = "#AE2012",
                mid = "#F7F7F7", faint = "grey92", light = "grey80", grey = "grey45")
# scico::scico(11, palette = "lapaz"), reversed to run light to dark; written
# out so that scico is not a dependency
.lapaz <- c("#FEF2F2", "#F1D5C4", "#CAB79D", "#A0A695", "#7B9A9E", "#5B8BA2", "#4177A1",
            "#315E98", "#27468B", "#212A78", "#190C64")
.okabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' The spiDE plot theme
#'
#' The theme of every spiDE plot: a black panel border, no grid lines, black
#' axis ticks, boxed facet strips and italic legend titles, with text scaled
#' by \code{rl}. Add it to your own ggplots so that they match spiDE's, and
#' change any element with a further \code{+ theme(...)}.
#'
#' @param rl a positive number, the text size relative to ggplot2's default.
#' @return a ggplot2 theme.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_spiDE()
#' @seealso [spiDEColours()]
#' @import ggplot2
#' @export
theme_spiDE <- function(rl = 1.1) {
  if (!is.numeric(rl) || length(rl) != 1L || !is.finite(rl) || rl <= 0) {
    stop("'rl' must be a single positive number", call. = FALSE)
  }
  theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.title = element_text(size = rel(rl) * 1.1),
          axis.text = element_text(size = rel(rl), colour = "black"),
          plot.title = element_text(size = rel(rl) * 1.2),
          strip.background = element_rect(fill = NA, colour = "black"),
          strip.text = element_text(size = rel(rl)),
          legend.text = element_text(size = rel(rl)),
          legend.title = element_text(size = rel(rl), face = "italic"))
}

#' The spiDE colours
#'
#' spiDE gives each role one colour and uses it in every plot: amber and violet
#' for a slope that rises or falls with niche density (and for the diverging
#' t scale), teal and rufous for the first and second level of the condition,
#' the scico "lapaz" ramp for magnitudes, and the Okabe-Ito colours for a
#' nuisance factor such as slide. All pairs stay distinct under the common
#' forms of colour blindness. Use them to make your own figures match, or pass
#' others to a ggplot2 scale to override a spiDE plot's.
#'
#' @param role \code{NULL} (the named fixed colours) or one of
#'   \code{"direction"}, \code{"condition"}, \code{"magnitude"},
#'   \code{"diverging"} and \code{"nuisance"}.
#' @param n the number of colours of a ramp, or of nuisance levels.
#' @return a character vector of colours.
#' @examples
#' spiDEColours("condition")
#' spiDEColours("magnitude", 5)
#' @seealso [theme_spiDE()]
#' @export
spiDEColours <- function(role = NULL, n = 9L) {
  if (is.null(role)) return(.spideCols)
  role <- match.arg(role, c("direction", "condition", "magnitude", "diverging", "nuisance"))
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("'n' must be a positive integer", call. = FALSE)
  }
  n <- as.integer(n)
  switch(role,
         direction = .spideCols[c("up", "down")],
         condition = .spideCols[c("first", "second")],
         magnitude = .magnitudeRamp(n),
         diverging = grDevices::colorRampPalette(.spideCols[c("down", "mid", "up")], space = "Lab")(n),
         nuisance = if (n <= length(.okabeIto)) .okabeIto[seq_len(n)] else scales::hue_pal()(n))
}

# The magnitude ramp, light to dark. skip drops that fraction of its light end:
# on a white panel the lowest values must stay visible (the tissue map).
.magnitudeRamp <- function(n = 256L, skip = 0) {
  r <- grDevices::colorRampPalette(.lapaz, space = "Lab")(1000L)
  r <- r[seq(floor(skip * 999) + 1, 1000)]
  r[round(seq(1, length(r), length.out = n))]
}

.scaleDirection <- function(labels, aesthetics = "colour", extra = NULL, name = NULL) {
  brk <- c("up", "down", names(extra))
  scale_discrete_manual(aesthetics, values = stats::setNames(unname(c(.spideCols[c("up", "down")], extra)), brk),
                        breaks = brk, labels = c(unname(labels[c("up", "down")]), names(extra)),
                        name = name, drop = FALSE)
}

.scaleCondition <- function(levels, aesthetics = "colour", name = NULL) {
  vals <- if (length(levels) == 2L) unname(.spideCols[c("first", "second")]) else rep("black", length(levels))
  scale_discrete_manual(aesthetics, values = stats::setNames(vals, levels), name = name)
}

.scaleT <- function(limit, name = "t") {
  scale_fill_gradient2(low = .spideCols[["down"]], mid = .spideCols[["mid"]], high = .spideCols[["up"]],
                       midpoint = 0, limits = c(-limit, limit), oob = scales::squish, name = name)
}

# plotmath labels: gene symbols in italics; quotes and backslashes escaped
.pmString <- function(x) gsub("([\"\\\\])", "\\\\\\1", x)
.itGene <- function(g) sprintf('italic("%s")', .pmString(g))
.itTriplet <- function(gene, index, niche) {
  sprintf('italic("%s")~~"%s | %s"', .pmString(gene), .pmString(index), .pmString(niche))
}
.parseLabels <- function(x) parse(text = x, keep.source = FALSE)
.fmtq <- function(q) {
  ifelse(!is.finite(q), "NA",
         ifelse(q < 1e-3, formatC(q, digits = 0, format = "e"), formatC(q, digits = 2, format = "g")))
}
