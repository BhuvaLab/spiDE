# Drawing pieces for the spiDE hex sticker (data-raw/sticker/make_sticker.R): a pointy-top hexagon (the Bioconductor
# sticker shape), comic halftone, offset-print ("glitch") lettering, a sagging orb web, speed streaks and the index cell.
# The look is inspired by the Spider-Verse films' comic style; no character, logo or emblem of them is used.
suppressPackageStartupMessages(library(ggplot2))
S3 <- sqrt(3)
pal <- list(navy = "#14284B", blue = "#1F4E99", sky = "#6FA8DC", teal = "#2A9D8F", verm = "#E76F51",
            purple = "#6A4C9C", pale = "#EAF1FB", grey = "#8FA3BF", white = "#FFFFFF")

hex_poly <- function(R = 1, cx = 0, cy = 0) {
  a <- (90 + 60 * (0:5)) * pi / 180
  data.frame(x = cx + R * cos(a), y = cy + R * sin(a))
}
in_hex <- function(x, y, R = 1) abs(x) <= S3 / 2 * R & abs(y) <= R - abs(x) / S3
circle <- function(cx, cy, r, n = 200) {
  a <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}
rand_in_hex <- function(n, R = 0.9, seed = 1) {
  set.seed(seed)
  x <- runif(n * 3, -R, R); y <- runif(n * 3, -R, R)
  k <- in_hex(x, y, R)
  data.frame(x = x[k], y = y[k])[seq_len(n), ]
}

# the sticker canvas: filled hexagon, content added by the recipe, outline drawn last by `finish()`
canvas <- function(fill) {
  ggplot() +
    geom_polygon(data = hex_poly(1), aes(x, y), fill = fill, colour = NA) +
    coord_fixed(xlim = c(-0.9, 0.9), ylim = c(-1.0, 1.0), expand = FALSE) +
    theme_void()
}
finish <- function(p, border, lw = 2.2) {
  p + geom_polygon(data = hex_poly(1), aes(x, y), fill = NA, colour = border, linewidth = lw)
}
wordmark <- function(p, y = -0.66, size = 17, col1 = "white", col2 = pal$verm) {
  p +
    annotate("text", x = 0.012, y = y, label = "spi", hjust = 1, size = size, fontface = "bold", colour = col1) +
    annotate("text", x = 0.012, y = y, label = "DE", hjust = 0, size = size, fontface = "bold", colour = col2)
}

# ---- Spider-Verse-style pieces: comic halftone, offset-print glitch, cell webs --------------------------
verse <- list(ink = "#0E0E14", red = "#E8253A", cyan = "#16D9F5", magenta = "#FF3FA4", yellow = "#FFE34D",
              paper = "#F4F1EC")

# comic halftone: a dot grid whose dot size grows towards (fx, fy); clipped to the hexagon
halftone <- function(step = 0.045, fx = 0.6, fy = -0.6, smax = 2.2, colour = "white", alpha = 0.18) {
  g <- expand.grid(x = seq(-1, 1, by = step), y = seq(-1, 1, by = step * S3 / 2))
  g$x <- g$x + ifelse(round(g$y / (step * S3 / 2)) %% 2 == 0, 0, step / 2)
  g <- g[in_hex(g$x, g$y, 0.97), ]
  d <- sqrt((g$x - fx)^2 + (g$y - fy)^2)
  g$s <- smax * pmax(0, 1 - d / 1.9)^1.4
  g <- g[g$s > 0.08, ]
  geom_point(data = g, aes(x, y, size = s), colour = colour, alpha = alpha, show.legend = FALSE)
}

# misregistered print: cyan and red copies offset behind the main layer
glitch_text <- function(p, label, x, y, size, colour, hjust = 0.5, fontface = "bold.italic", off = 0.014,
                        c1 = verse$cyan, c2 = verse$red, angle = 0) {
  p +
    annotate("text", x = x - off, y = y + off * 0.6, label = label, size = size, colour = c1, hjust = hjust,
             fontface = fontface, angle = angle) +
    annotate("text", x = x + off, y = y - off * 0.6, label = label, size = size, colour = c2, hjust = hjust,
             fontface = fontface, angle = angle) +
    annotate("text", x = x, y = y, label = label, size = size, colour = colour, hjust = hjust,
             fontface = fontface, angle = angle)
}

# k-nearest-neighbour edges of a cell layout (the neighbourhood network)
knn_edges <- function(pts, k = 3) {
  D <- as.matrix(dist(pts[, c("x", "y")]))
  diag(D) <- Inf
  e <- do.call(rbind, lapply(seq_len(nrow(pts)), function(i) {
    j <- order(D[i, ])[seq_len(k)]
    data.frame(i = pmin(i, j), j = pmax(i, j))
  }))
  e <- unique(e)
  data.frame(x = pts$x[e$i], y = pts$y[e$i], xend = pts$x[e$j], yend = pts$y[e$j],
             d = (sqrt(pts$x[e$i]^2 + pts$y[e$i]^2) + sqrt(pts$x[e$j]^2 + pts$y[e$j]^2)) / 2)
}

# an orb web whose knots are cells: spokes from the index cell, rings joining neighbouring spokes
cell_web <- function(cx = 0, cy = 0, rings = c(0.15, 0.27, 0.39, 0.51, 0.63), spokes = 11, jitter = 0.018,
                     seed = 3) {
  set.seed(seed)
  a0 <- seq(0, 2 * pi, length.out = spokes + 1)[-1] + 0.2
  nodes <- expand.grid(s = seq_len(spokes), k = seq_along(rings))
  nodes$a <- a0[nodes$s] + rnorm(nrow(nodes), 0, 0.05)
  nodes$r <- rings[nodes$k] * (1 + rnorm(nrow(nodes), 0, 0.05))
  nodes$x <- cx + nodes$r * cos(nodes$a) + rnorm(nrow(nodes), 0, jitter)
  nodes$y <- cy + nodes$r * sin(nodes$a) + rnorm(nrow(nodes), 0, jitter)
  key <- function(s, k) (k - 1) * spokes + s
  sp <- do.call(rbind, lapply(seq_len(spokes), function(s) {
    idx <- key(s, seq_along(rings))
    data.frame(x = c(cx, nodes$x[idx[-length(idx)]]), y = c(cy, nodes$y[idx[-length(idx)]]),
               xend = nodes$x[idx], yend = nodes$y[idx], k = seq_along(rings))
  }))
  rg <- do.call(rbind, lapply(seq_along(rings), function(k) {
    i <- key(seq_len(spokes), k); j <- key(c(2:spokes, 1), k)
    data.frame(x = nodes$x[i], y = nodes$y[i], xend = nodes$x[j], yend = nodes$y[j], k = k)
  }))
  list(nodes = nodes, edges = rbind(sp, rg), rings = rings)
}

# ---- neon pieces after the two 2026 mockups (spiDE_sticker_v1/v2 on RDM) -------------------------------
neon <- list(violet = "#2A0E52", violet2 = "#4A1786", magenta = "#FF1F9E", pink = "#FF5FC1", blue = "#2D3FE0",
             royal = "#1E2BB8", cyan = "#2FD2FF", ink = "#140A2E", paper = "#FBFAFD")

# a glowing stroke: the same path drawn wide and faint, then narrow and solid
glow_path <- function(df, colour, lw = 1.2, grp = NULL, halo = c(5, 3.2, 2), halo_alpha = c(0.07, 0.14, 0.28)) {
  aes0 <- if (is.null(grp)) aes(x, y) else aes(x, y, group = .data[[grp]])
  c(lapply(seq_along(halo), function(i) geom_path(data = df, mapping = aes0, colour = colour,
                                                   linewidth = lw * halo[i], alpha = halo_alpha[i], lineend = "round")),
    list(geom_path(data = df, mapping = aes0, colour = colour, linewidth = lw, lineend = "round")))
}
glow_seg <- function(df, colour, lw = 1.2, halo = c(5, 3.2, 2), halo_alpha = c(0.07, 0.14, 0.28)) {
  c(lapply(seq_along(halo), function(i) geom_segment(data = df, aes(x, y, xend = xend, yend = yend), colour = colour,
                                                      linewidth = lw * halo[i], alpha = halo_alpha[i], lineend = "round")),
    list(geom_segment(data = df, aes(x, y, xend = xend, yend = yend), colour = colour, linewidth = lw, lineend = "round")))
}

# sagging web threads: between neighbouring spokes at each ring radius, bowed towards the centre
web_threads <- function(cx, cy, angles, rings, sag = 0.12, n = 30) {
  out <- list(); id <- 0
  for (r in rings) for (i in seq_along(angles)) {
    a1 <- angles[i]; a2 <- if (i < length(angles)) angles[i + 1] else angles[1] + 2 * pi
    t <- seq(0, 1, length.out = n)
    a <- a1 + (a2 - a1) * t
    rr <- r * (1 - sag * sin(pi * t))
    id <- id + 1
    out[[id]] <- data.frame(x = cx + rr * cos(a), y = cy + rr * sin(a), id = id, r = r)
  }
  do.call(rbind, out)
}

# diagonal speed streaks (comic motion lines), kept inside the hexagon and away from a central disc
streaks <- function(n = 26, angle = -32, cx = 0, cy = 0, keep_out = 0.42, seed = 2, cols = NULL) {
  set.seed(seed)
  th <- angle * pi / 180
  x0 <- runif(n * 4, -0.95, 0.95); y0 <- runif(n * 4, -0.95, 0.95); L <- runif(n * 4, 0.12, 0.42)
  d <- data.frame(x = x0, y = y0, xend = x0 + L * cos(th), yend = y0 + L * sin(th), lw = runif(n * 4, 0.6, 3.2))
  mid <- sqrt(((d$x + d$xend) / 2 - cx)^2 + ((d$y + d$yend) / 2 - cy)^2)
  d <- d[in_hex(d$x, d$y, 0.9) & in_hex(d$xend, d$yend, 0.9) & mid > keep_out, ]
  d <- d[seq_len(min(n, nrow(d))), ]
  d$col <- if (is.null(cols)) neon$blue else sample(cols, nrow(d), replace = TRUE)
  d
}

# the index cell of the mockups: a many-sided cell, a glowing hexagonal core filled with halftone, spokes between
index_cell <- function(cx, cy, r = 0.21, sides = 11, core = 0.105, rot = 0.15) {
  a <- seq(0, 2 * pi, length.out = sides + 1)[-1] + rot
  outer <- data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
  h <- hex_poly(core, cx, cy)
  ah <- (90 + 60 * (0:5)) * pi / 180
  spokes <- data.frame(x = cx + core * cos(ah), y = cy + core * sin(ah),
                       xend = cx + 0.92 * r * cos(ah), yend = cy + 0.92 * r * sin(ah))
  g <- expand.grid(x = seq(-core, core, by = 0.016), y = seq(-core, core, by = 0.016))
  g <- g[in_hex(g$x, g$y, core * 0.9), ]
  g$x <- g$x + cx; g$y <- g$y + cy
  list(outer = outer, core = h, spokes = spokes, dots = g)
}
