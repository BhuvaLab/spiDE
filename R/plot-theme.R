# The spiDE plotting theme and colours (0.99.34; palette "V", option chosen by the PI).
# One colour per role, the same in every plot, checked under simulated
# deuteranopia, protanopia and tritanopia:
#   direction (up / down, and the diverging t scale): gold / violet
#   the two conditions (first level / second): dark violet / pink
#   magnitude (density, R^2, counts): viridis "magma", light to dark, on the
#     light page: rev(viridisLite::magma(11, begin = 0.08, end = 0.98))
#   spatial (tissue coordinates): magma(begin = 0.25, end = 0.98), dark to
#     bright, for the dark panel of plotNicheMap()
#   a nuisance factor (slide, batch): Okabe-Ito
# The theme is vissE::bhuvad_theme() -- black panel border, no grid, boxed
# strips, italic legend titles -- with black axis ticks. A spatial plot darkens
# only its panel (the sticker's ink); legends, titles and axes keep the page's.

.spideCols <- c(up = "#C99A16", down = "#5B2C9E", first = "#3E1A78", second = "#E05A9E",
                mid = "#F7F7F7", faint = "grey92", light = "grey80", grey = "grey45",
                ink = "#140B2E", dim = "#3B3158")
# Near white, then viridisLite::magma(11, begin = 0.08, end = 0.98) reversed to
# run light to dark (cream, orange, pink, violet, near black); written out so
# that viridisLite is not a dependency. Near white so that low values look empty.
.magma <- c("#FDFAF4", "#FCF4B6", "#FECA8D", "#FE9F6D", "#F8745C", "#E44F64", "#C23B75", "#9C2E7F",
            "#782281", "#55137D", "#2E1162", "#0F0B2C")
# viridisLite::magma(11, begin = 0.25, end = 0.98), dark to bright
.magmaSpatial <- c("#51127C", "#6D1D81", "#8B2981", "#A9337D", "#C73D73", "#E34E65", "#F46A5C",
                   "#FC8C63", "#FEAF78", "#FED294", "#FCF4B6")
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
#' spiDE gives each role one colour and uses it in every plot: gold and violet
#' for a slope that rises or falls with niche density (and for the diverging
#' t scale), dark violet and pink for the first and second level of the
#' condition, the viridis "magma" ramp for magnitudes (light to dark on the
#' page), a brighter part of the same ramp (dark to bright) for plots in
#' tissue coordinates, drawn on a dark panel, and the Okabe-Ito colours for a
#' nuisance factor such as slide. The colours were checked under simulated
#' deuteranopia, protanopia and tritanopia. Use them to make your own figures
#' match, or pass
#' others to a ggplot2 scale to override a spiDE plot's.
#'
#' @param role \code{NULL} (the named fixed colours) or one of
#'   \code{"direction"}, \code{"condition"}, \code{"magnitude"} (the
#'   light-panel ramp, light to dark), \code{"spatial"} (the dark-panel ramp
#'   of tissue plots, dark to bright), \code{"diverging"} and
#'   \code{"nuisance"}.
#' @param n the number of colours of a ramp, or of nuisance levels.
#' @return a character vector of colours.
#' @examples
#' spiDEColours("condition")
#' spiDEColours("magnitude", 5)
#' spiDEColours("spatial", 5)
#' @seealso [theme_spiDE()]
#' @export
spiDEColours <- function(role = NULL, n = 9L) {
  if (is.null(role)) return(.spideCols)
  role <- match.arg(role, c("direction", "condition", "magnitude", "spatial", "diverging", "nuisance"))
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("'n' must be a positive integer", call. = FALSE)
  }
  n <- as.integer(n)
  switch(role,
         direction = .spideCols[c("up", "down")],
         condition = .spideCols[c("first", "second")],
         magnitude = .magnitudeRamp(n),
         spatial = .magnitudeRamp(n, spatial = TRUE),
         diverging = grDevices::colorRampPalette(.spideCols[c("down", "mid", "up")], space = "Lab")(n),
         nuisance = if (n <= length(.okabeIto)) .okabeIto[seq_len(n)] else scales::hue_pal()(n))
}

# The magnitude ramp, light to dark (magma). spatial = TRUE gives the
# dark-panel ramp of tissue plots, dark to bright. skip drops that fraction of
# the first end of the ramp.
.magnitudeRamp <- function(n = 256L, skip = 0, spatial = FALSE) {
  base <- if (spatial) .magmaSpatial else .magma
  if (skip == 0 && n == length(base)) return(base)
  r <- grDevices::colorRampPalette(base, space = "Lab")(1000L)
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
