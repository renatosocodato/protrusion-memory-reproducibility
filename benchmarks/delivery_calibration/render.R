args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply ring and decision JSON files and an output PDF")
ring <- jsonlite::read_json(args[1])$cases
ai <- jsonlite::read_json(args[2])$domains
sage <- "#c8e2d5"
lilac <- "#e0d6ed"
rose <- "#f2cedb"
green <- "#2d6857"
violet <- "#644972"
pink <- "#9a4e66"
ink <- "#262935"
labels <- c("Actin-mediated", "Direct only", "Transferred effect", "Incomplete block", "Transfer fails")
y <- rev(seq_along(labels))

grDevices::cairo_pdf(args[3], width = 10, height = 6.7, family = "Arial")
layout(matrix(1:4, 2, byrow = TRUE), widths = c(1, 1), heights = c(.95, 1.05))
par(mgp = c(2.1, .55, 0), las = 1, bty = "n", family = "Arial", col.axis = ink,
    col.lab = ink, cex.axis = .83)
title_panel <- function(letter, text) {
  mtext(paste(letter, text, sep = "    "), side = 3, adj = 0, line = .7, font = 2, col = ink)
}

par(mar = c(2.4, .8, 2.4, .7))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
title_panel("A", "Calibrate the same spatial mode")
starts <- c(.03, .36, .69)
for (i in seq_along(starts)) {
  x <- starts[i]
  rect(x, .38, x + .28, .78, col = c(sage, lilac, sage)[i], border = NA)
  text(x + .14, .64, c("Angular\nstandards", "Actin-disabled\nmechanics", "Actin-intact\nreserved test")[i],
       cex = .86, col = ink, font = 2)
  text(x + .14, .45, c("8 preparations", "10 preparations", "12 preparations")[i], cex = .71, col = ink)
}
arrows(starts[1:2] + .282, .58, starts[2:3] - .006, .58, length = .045, col = violet)
text(.5, .16, "Fixed upstream input · mode 2 · 10–30 s", cex = .94, col = ink)

par(mar = c(3.3, 8.05, 2.4, .75))
plot(NULL, xlim = c(-.16, .62), ylim = c(.55, 5.45), yaxt = "n", xaxt = "n",
     xlab = expression("Signed mode-2 deformation ("*mu*"m)"), ylab = "")
title_panel("B", "Mechanics can mimic an actin effect")
rect(-.08, .44, .08, 5.56, col = sage, border = NA)
abline(v = 0, lty = 2, col = green)
axis(1, at = seq(-.1, .6, by = .1))
axis(2, at = y, labels = labels, tick = FALSE)
for (i in seq_along(ring)) {
  case <- ring[[i]]$representative
  wide <- unlist(case$shape_interval_with_leakage)
  col <- if (i == 5) pink else violet
  if (length(wide)) segments(wide[1], y[i], wide[2], y[i], col = if (i == 5) rose else lilac, lwd = 6)
  small <- unlist(case$shape_odd_interval)
  segments(small[1], y[i], small[2], y[i], col = col, lwd = 1.5)
  points(case$shape_odd_mean, y[i], pch = if (length(wide)) 19 else 4, col = col, cex = .85)
}

par(mar = c(3.3, 8.05, 2.4, .75))
plot(NULL, xlim = c(0, 60), ylim = c(-.7, 5.62), yaxt = "n", xlab = "Independent synthetic cohorts", ylab = "")
title_panel("C", "Recovery depends on transfer")
axis(2, at = y, labels = labels, tick = FALSE)
for (i in seq_along(ring)) {
  counts <- ring[[i]]$counts
  get <- function(key) if (is.null(counts[[key]])) 0 else counts[[key]]
  called <- c(get("deformation_beyond_leakage"), get("mechanical_control_unqualified"))
  end <- cumsum(c(0, called, 60 - sum(called)))
  rect(head(end, -1), y[i] - .3, tail(end, -1), y[i] + .3,
       col = c(if (i == 5) rose else sage, "#a8a1b4", lilac), border = NA)
}
legend(-2, .60, xpd = NA, cex = .67, bty = "n", ncol = 2,
       fill = c(sage, rose, lilac, "#a8a1b4"), border = NA,
       legend = c("Recovered mediation", "False attribution", "No excess deformation", "Control unqualified"))

par(mar = c(3.3, 4.3, 2.4, .95))
plot(NULL, xlim = c(0, 1.05), ylim = c(.25, 1.06), xlab = "Mechanical controls acquired per case",
     ylab = "Fraction of cases decided")
title_panel("D", "Acquire controls on ambiguous cases")
abline(h = 1, lty = 3, col = "#a8a1b4")
points_for <- function(domain, style, color) {
  for (name in c("always_calibrate", "analytic", "learned")) {
    item <- ai[[domain]][[name]]
    points(item$measurements_per_case, item$coverage, pch = style, bg = color, col = ink, cex = 1.2)
  }
}
points_for("transfer", 21, sage)
points_for("nontransfer", 24, rose)
text(.67, .958, "Learned", pos = 1, cex = .75, offset = .4)
text(.71, 1, "Analytic", pos = 2, cex = .75, offset = .45)
text(1, 1, "Universal", pos = 1, cex = .75, offset = .55)
legend("left", legend = c("Transfer qualified", "Transfer fails in twins"), pt.bg = c(sage, rose),
       pch = c(21, 24), bty = "n", cex = .76)
invisible(dev.off())
