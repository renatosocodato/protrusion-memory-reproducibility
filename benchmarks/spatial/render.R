args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply geometry JSON, recovery JSON and output PDF")
geometry <- jsonlite::read_json(args[1])
recovery <- jsonlite::read_json(args[2])
frames <- geometry$frames
time <- vapply(frames, function(x) x$time_s, 0) / 60
xy <- t(vapply(frames, function(x) unlist(x$centroid_xy_um), numeric(2)))
xy <- sweep(xy, 2, xy[1, ])
modes <- c("2", "3", "4")
centered <- lapply(modes, function(m) {
    values <- t(vapply(frames, function(x) {
        if (is.null(x$centered_radius_modes_um)) return(c(NA_real_, NA_real_))
        unlist(x$centered_radius_modes_um[[m]])
    }, numeric(2)))
    reference <- colMeans(values[1:6, ], na.rm = TRUE)
    sqrt(rowSums(sweep(values, 2, reference)^2))
})
pick <- function(r, model) {
    r$clock == "G01" && r$actin_diffusivity == .015 && r$model == model
}
shared <- Filter(function(r) pick(r, "common_gain"), recovery$results)
free <- Filter(function(r) pick(r, "mode_specific_gains"), recovery$results)
phase <- vapply(shared, function(r) r$phase, 0)
fraction <- function(rows, truth) {
    vapply(rows, function(r) r$confusion[[truth]][[truth]], 0) / recovery$test_paths_per_class
}
dark <- c("#2d6857", "#644972", "#9a4e66")
light <- c("#bfdfcd", "#dcd0eb", "#f2cdd9")
grDevices::cairo_pdf(args[3], width = 11.2, height = 3.7, family = "Arial")
par(mfrow = c(1, 3), mar = c(4.2, 4.35, 2.35, 0.75), family = "Arial",
    las = 1, bty = "l", mgp = c(2.6, .8, 0), cex.lab = .94)
plot(xy, type = "l", lwd = 2.1, col = dark[1], asp = 1,
     xlab = expression(Delta*x~"centroid ("*mu*m*")"),
     ylab = expression(Delta*y~"centroid ("*mu*m*")"))
points(xy[c(1,nrow(xy)), , drop = FALSE], pch = 21, bg = light[1:2], col = dark[1], cex = 1.35)
title("Centroid displacement", adj = 0, cex.main = 1.01)
legend("topleft", c("First frame","Last frame"), pch = 21, pt.bg = light[1:2], bty = "n", cex = .78)
matplot(time, do.call(cbind, centered), type = "l", lty = c(1, 2, 3),
        col = dark, lwd = 2.1, xlab = "Time (min)",
        ylab = expression("Change in mode coefficient ("*mu*m*")"))
title("Centroid-frame deformation", adj = 0, cex.main = 1.01)
legend("topleft", paste("Mode", modes), lty = c(1,2,3), col = dark, bty = "n", cex = .78)
plot(phase, 100*fraction(shared, "mediated"), type = "b", pch = 21, bg = light[1], col = dark[1],
     lwd = 2, ylim = c(0, 73), xlab = "Pulse phase (fraction of frame)", ylab = "Correct assignments (%)")
lines(phase, 100*fraction(shared, "direct"), type = "b", pch = 22, bg = light[2],
      col = dark[2], lwd = 2)
lines(phase, 100*fraction(free, "direct"), lwd = 2, col = dark[3], lty = 3)
title("Pathway recovery at 65.589 s", adj = 0, cex.main = 1.01)
legend("topleft", c("Mediated, shared gain", "Direct, shared gain", "Either, separate gains"),
       pch = c(21,22,NA), pt.bg = light, lty = c(1,1,3), col = dark, bty = "n", cex = .77)
invisible(dev.off())
