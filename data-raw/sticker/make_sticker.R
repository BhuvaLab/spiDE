# The spiDE hex sticker: one index cell (the halftone-core polygon) at the hub of a web whose knots are its neighbouring
# cells; six neighbours on glowing threads fade from magenta (near) to blue (far), the niche influence spiDE tests.
# Writes spiDE.png (900 x 1039 px, about 520 dpi at the standard 43.99 x 50.8 mm hexagon, transparent corners) and
# spiDE.svg (43.99 x 50.8 mm, text as outlines so it renders without the font), each under BiocStickers' 500 KB limit,
# into the folder given, by default this one:
#   Rscript data-raw/sticker/make_sticker.R [out_dir]
# Needs ggplot2 and ragg.
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
out <- if (length(args)) args[1] else here
source(file.path(here, "sticker-helper.R"))

cx <- 0
cy <- -0.22
deg <- (0:8) * 40 + 20
ang <- deg * pi / 180
rings <- c(0.31, 0.42, 0.53, 0.63)
inside <- function(x, y) in_hex(x, y, 0.92) & y < 0.42

# the web, faint behind: threads (whole arcs only), spokes, and a small knot (a cell) where threads meet a spoke
th <- web_threads(cx, cy, ang, rings, sag = 0.13)
th <- th[ave(inside(th$x, th$y), th$id, FUN = all) == 1, ]
ends <- do.call(rbind, lapply(split(th, th$id), function(a) a[c(1, nrow(a)), c("x", "y", "r")]))
ends$s <- (round(atan2(ends$y - cy, ends$x - cx) * 180 / pi - 20) %% 360) / 40 + 1   # spoke index of each end
kn <- unique(data.frame(s = round(ends$s), r = ends$r))
kn$x <- cx + kn$r * cos(ang[kn$s])
kn$y <- cy + kn$r * sin(ang[kn$s])
# each spoke runs from the index cell to its outermost knot, never past the last thread
reach <- tapply(kn$r, kn$s, max)
sp <- data.frame(s = seq_along(ang), r = ifelse(is.na(reach[as.character(seq_along(ang))]), rings[1],
                                                  reach[as.character(seq_along(ang))]))
sp$x <- cx + 0.21 * cos(ang[sp$s])
sp$y <- cy + 0.21 * sin(ang[sp$s])
sp$xend <- cx + sp$r * cos(ang[sp$s])
sp$yend <- cy + sp$r * sin(ang[sp$s])

# the neighbour cells: on web knots, nearer = stronger influence
tier <- data.frame(ring = 2:4, col = c(neon$magenta, "#C66BFF", "#4C5BFF"),
                   size = c(10.5, 7.8, 5.8), lw = c(2.7, 1.7, 1.1))
nb <- data.frame(deg = c(140, 300, 20, 180, 60, 220), ring = c(2, 2, 3, 3, 4, 4))
nb <- merge(nb, tier, by = "ring")
nb$x <- cx + rings[nb$ring] * cos(nb$deg * pi / 180)
nb$y <- cy + rings[nb$ring] * sin(nb$deg * pi / 180)
nb$x0 <- cx + 0.21 * cos(nb$deg * pi / 180)
nb$y0 <- cy + 0.21 * sin(nb$deg * pi / 180)

st <- streaks(30, angle = -30, cx = cx, cy = cy, keep_out = 0.27, seed = 6,
              cols = c(neon$royal, neon$magenta, "#7A2BD0", neon$blue))
st <- st[st$y < 0.42 & st$yend < 0.42, ]
# no streak passes under a neighbour cell
near <- sapply(seq_len(nrow(st)), function(i) {
  t <- seq(0, 1, length.out = 12)
  px <- st$x[i] + t * (st$xend[i] - st$x[i]); py <- st$y[i] + t * (st$yend[i] - st$y[i])
  min(outer(px, nb$x, "-")^2 + outer(py, nb$y, "-")^2) < 0.09^2
})
st <- st[!near, ]
ic <- index_cell(cx, cy)

p <- canvas(neon$violet) +
  halftone(fx = -0.7, fy = -0.7, colour = neon$pink, alpha = 0.22, smax = 2.2) +
  scale_size_identity() +
  geom_segment(data = st, aes(x, y, xend = xend, yend = yend, linewidth = lw), colour = st$col, alpha = 0.55,
               lineend = "round") +
  scale_linewidth_identity() +
  glow_path(th, neon$royal, lw = 0.8, grp = "id", halo_alpha = c(0.03, 0.07, 0.14)) +
  glow_seg(sp, neon$royal, lw = 0.9, halo_alpha = c(0.03, 0.07, 0.14)) +
  geom_point(data = kn, aes(x, y), colour = "#6C78FF", size = 1.6, alpha = 0.9)
# each neighbour's thread to the index cell, glowing in its own colour and weight
for (i in seq_len(nrow(nb))) {
  p <- p + glow_seg(data.frame(x = nb$x0[i], y = nb$y0[i], xend = nb$x[i], yend = nb$y[i]), nb$col[i], lw = nb$lw[i])
}
p <- p +
  geom_point(data = nb, aes(x, y, size = size + 4), colour = nb$col, alpha = 0.22) +
  geom_point(data = nb, aes(x, y, size = size), shape = 21, fill = nb$col, colour = neon$cyan, stroke = 2) +
  geom_polygon(data = ic$outer, aes(x, y), fill = "#24145C", colour = neon$cyan, linewidth = 1.8) +
  geom_segment(data = ic$spokes, aes(x, y, xend = xend, yend = yend), colour = neon$magenta, linewidth = 0.8) +
  geom_polygon(data = ic$core, aes(x, y), fill = neon$magenta, colour = NA, alpha = 0.35) +
  geom_point(data = ic$dots, aes(x, y), colour = neon$pink, size = 0.7) +
  glow_path(rbind(ic$core, ic$core[1, ]), "#FFD3EE", lw = 1.1)
p <- glitch_text(p, "spiDE", -0.03, 0.55, 24, neon$magenta, off = 0.013, c1 = neon$cyan, c2 = neon$blue, angle = 2)
sticker <- finish(p, neon$magenta, 2.6)


# the canvas is exactly the hexagon (its border stroke included), with the 2:sqrt(3) sticker aspect
r_out <- 1.015
sticker <- suppressMessages(sticker + coord_fixed(xlim = c(-1, 1) * r_out * S3 / 2, ylim = c(-1, 1) * r_out,
                                                  expand = FALSE))
w_in <- 5                        # the design was drawn on a 5 inch wide canvas; point, line and text sizes follow it
h_in <- w_in * 2 / S3
ragg::agg_png(file.path(out, "spiDE.png"), width = 900, height = 1039, units = "px", res = 1039 / h_in,
              background = "transparent")
print(sticker)
invisible(dev.off())
grDevices::svg(file.path(out, "spiDE.svg"), width = w_in, height = h_in, bg = "transparent")
print(sticker)
invisible(dev.off())
# print size 43.99 x 50.8 mm (the viewBox stays, so the drawing scales as one piece), then a light minify: hex colours
# for Cairo's rgb(% ...), its default attributes dropped, path coordinates to 0.1 pt (about 0.01 mm printed) and every
# other number to 0.01
f <- file.path(out, "spiDE.svg")
s <- paste(readLines(f), collapse = "\n")
s <- sub('width="[^"]*" height="[^"]*"', 'width="43.99mm" height="50.8mm"', s)
m <- gregexpr("rgb\\(([0-9.]+)%, ([0-9.]+)%, ([0-9.]+)%\\)", s, perl = TRUE)
regmatches(s, m) <- lapply(regmatches(s, m), function(x) vapply(x, function(v) {
  p <- as.numeric(regmatches(v, gregexpr("[0-9.]+", v))[[1]])
  sprintf("#%02X%02X%02X", round(p * 2.55)[1], round(p * 2.55)[2], round(p * 2.55)[3])
}, character(1)))
s <- gsub(' fill-rule="nonzero"| stroke-miterlimit="10"', "", s)
m <- gregexpr(' d="[^"]*"', s)
regmatches(s, m) <- lapply(regmatches(s, m), function(x) gsub("([0-9]+\\.[0-9])[0-9]+", "\\1", x))
s <- gsub("([0-9]+\\.[0-9]{2})[0-9]+", "\\1", s)
writeLines(s, f)
message("wrote ", file.path(out, c("spiDE.png", "spiDE.svg")), collapse = " ")
