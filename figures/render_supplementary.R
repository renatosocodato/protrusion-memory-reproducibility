#!/usr/bin/env Rscript
# Draw the seven supplementary figures from frozen source outputs.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(here, "bj_style.R"))

read_data <- function(name) {
  read.csv(file.path(data_dir, name), stringsAsFactors = FALSE)
}
read_result <- function(name) {
  jsonlite::fromJSON(file.path(data_dir, name), simplifyVector = FALSE)
}

supp_export <- function(plots, name, ncol, height, footer = NULL,
                        footer_height = 9, widths = rep(1, ncol)) {
  nrow <- length(plots) / ncol
  plate_plots <- if (is.null(footer)) plots else lapply(plots, `+`, no_legend)
  combined <- labelled_grid(plate_plots, ncol, widths)
  if (!is.null(footer)) {
    combined <- plot_grid(combined, footer, ncol = 1,
                          rel_heights = c(height - footer_height, footer_height))
  }
  export(combined, name, width = 180, height = height)
  panel_dir <- file.path(output_dir, "panels")
  dir.create(panel_dir, showWarnings = FALSE)
  for (index in seq_along(plots)) {
    label <- letters[index]
    export(tag_panel(plots[[index]], label), paste0(name, label),
           width = 180 * widths[(index - 1) %% ncol + 1] / sum(widths),
           height = height/nrow, directory = panel_dir)
  }
}
panel <- function(p, x, y, title = NULL) {
  p + labs(x = x, y = y, title = title) + th()
}

# Supplementary Figure 1: paired timeseries from the source workbook, and
# dimensionless load trajectories on their original model clock.
amiri <- read_result("S01_amiri_traces.json")
series <- c(sensor_calyculin = "TS + calyculin", sensor_blebbistatin = "TS + blebbistatin",
            control_calyculin = "CTS + calyculin")
paired <- do.call(rbind, lapply(names(series), function(key) {
  x <- as.data.frame(do.call(rbind, lapply(amiri$conditions[[key]]$values, unlist)))
  names(x) <- c("time", "traction", "efficiency")
  x$segment <- cumsum(c(TRUE, diff(x$time) > 5))
  x$efficiency <- 100 * x$efficiency
  x$condition <- unname(series[[key]])
  x
}))
paired$condition <- factor(paired$condition, levels = unname(series))
load <- read_data("S01_load.csv")
load$drive <- factor(load$drive, levels = c("autonomous", "periodic", "stochastic"),
                     labels = c("Autonomous", "Periodic", "Stochastic"))

paired_plot <- function(metric, label, title) {
  color <- setNames(c(forest, plum, green_ink), unname(series))
  p <- ggplot(paired, aes(time, .data[[metric]], colour = condition, linetype = condition,
                         group = interaction(condition, segment))) +
    geom_line(linewidth = trace_width) + geom_point(size = .8) +
    scale_colour_manual(values = color) +
    scale_linetype_manual(values = comparison_lines) +
    guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
  panel(p, "Time (min)", label, title)
}
a <- paired_plot("traction", "Mean traction (Pa)", "Traction")
b <- paired_plot("efficiency", "FRET efficiency (%)", "Sensor response")
c <- panel(ggplot(load, aes(time, load, colour = drive, linetype = drive)) +
             geom_line(linewidth = trace_width) +
             scale_colour_manual(values = c(forest, plum, green_ink)) +
             scale_linetype_manual(values = comparison_lines) +
             guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)),
           "Assigned model time", "Model load", "Synthetic load")
line_legend <- function(plot, title = NULL) {
  footer_key(plot + guides(colour = guide_legend(nrow = 1),
                           linetype = guide_legend(nrow = 1)), title)
}
keys <- plot_grid(line_legend(a, "Traction and FRET conditions"),
                  line_legend(c, "Load drives"), nrow = 1, rel_widths = c(1.75, 1.25))
supp_export(list(a, b, c), "S01", ncol = 3, height = 66,
            footer = keys, footer_height = 13)

# Supplementary Figure 2: time trace, state projection and tangent estimates.
trace <- read_data("S02_state.csv")
trace_short <- subset(trace, time <= 60)
long <- stack(trace_short[c("q", "a", "x")])
long$time <- rep(trace_short$time, 3)
a <- panel(ggplot(long, aes(time, values, colour = ind, linetype = ind)) +
             geom_line(linewidth = trace_width) +
             scale_colour_manual(values = c(forest, plum, green_ink)) +
             scale_linetype_manual(values = comparison_lines) +
             coord_cartesian(ylim = c(.2, .85)),
           "Assigned model time", "Normalized proxy", "Autonomous trace")
b <- panel(ggplot(trace, aes(q, m)) + geom_path(colour = forest, alpha = .45, linewidth = .35),
           "Gate proxy, q", "Feedback proxy, m", "State projection")
chaos <- read_result("chaos_exponents.json")
cases <- c("low_control", "canonical", "feedback", "initial_b", "initial_c", "half_step")
case_name <- c("r = 10", "Lorenz", "r = 28, A", "r = 28, B", "r = 28, C", "r = 28, h/2")
exponents <- data.frame(case = factor(case_name, levels = rev(case_name)),
                        value = vapply(cases, function(key) chaos$cases[[key]]$spectrum[[1]], 0),
                        role = c("control", "canonical", rep("modified", 4)))
c <- panel(ggplot(exponents, aes(value, case, fill = role, shape = role)) +
             geom_segment(aes(x = 0, xend = value, yend = case),
                          colour = muted_ink, linewidth = .33) +
             geom_vline(xintercept = 0, colour = muted_ink, linewidth = reference_width) +
             geom_point(size = 2.2, colour = ink, stroke = mark_stroke) +
             scale_fill_manual(values = c(control = green_pale, canonical = pink_fill, modified = lilac_fill)) +
             scale_shape_manual(values = c(control = 21, canonical = 22, modified = 23)) +
             coord_cartesian(xlim = c(-.75, 1.15)),
           "Exponent / model time", NULL, "Tangent estimate") + no_legend
keys <- plot_grid(line_legend(a), NULL, NULL, nrow = 1)
supp_export(list(a, b, c), "S02", ncol = 3, height = 63, footer = keys)

# Supplementary Figure 3: one trace per model and regime. Observation markers
# appear only at the prescribed 60-second displacement sample times.
org <- read_data("S03_trajectory.csv")
org$model <- factor(org$model, levels = c("uniform", "finite_recruitment", "demixing_analogue"),
                    labels = c("Uniform", "Finite recruitment", "Demixing analogue"))
trajectory_panel <- function(concentration, metric, label) {
  data <- subset(org, m == concentration)
  p <- ggplot(data, aes(time, .data[[metric]], colour = model, linetype = model)) +
    annotate("rect", xmin = 20, xmax = 200, ymin = -Inf, ymax = Inf,
             fill = lilac_tint) +
    geom_line(linewidth = trace_width) +
    scale_colour_manual(values = c(forest, green_ink, plum)) +
    scale_linetype_manual(values = comparison_lines)
  if (metric == "displacement") {
    p <- p + geom_point(data = subset(data, time %% 60 == 0), size = 1.35)
  }
  panel(p, "Assigned time (s)", label, paste0("m = ", formatC(concentration, digits = 1, format = "f"))) +
    theme(legend.key.width = grid::unit(8, "mm")) +
    guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
}
series_panels <- function(concentration) {
  list(trajectory_panel(concentration, "contrast", "Effective contrast"),
       trajectory_panel(concentration, "actin", "Actin output"),
       trajectory_panel(concentration, "displacement", expression(Delta*x)))
}
row1 <- series_panels(.3)
row2 <- series_panels(1)
supp_export(c(row1, row2), "S03", ncol = 3, height = 117,
            footer = line_legend(row1[[1]]))

# Supplementary Figure 4: within-generating-class normalized confusion grids.
conf <- read_data("S04_heldout_predictions.csv")
classes <- c("uniform", "finite_recruitment", "demixing_analogue")
short_classes <- c("Uniform", "Recruitment", "Demixing")
matrix_panel <- function(concentration, mode, mode_label) {
  data <- subset(conf, mean_concentration == concentration & channel == mode)
  ct <- table(factor(data$true_class, levels = classes), factor(data$predicted_class, levels = classes))
  stopifnot(all(rowSums(ct) == 300))
  tbl <- as.data.frame(as.table(prop.table(ct, 1)))
  names(tbl) <- c("generated", "predicted", "fraction")
  tbl$generated <- factor(short_classes[match(tbl$generated, classes)], levels = rev(short_classes))
  tbl$predicted <- factor(short_classes[match(tbl$predicted, classes)], levels = short_classes)
  p <- ggplot(tbl, aes(predicted, generated, fill = fraction)) +
    geom_tile(linewidth = .35, colour = "white") +
    geom_text(aes(label = sprintf("%.2f", fraction)), colour = ink,
              family = figure_font, size = 2.9) +
    scale_fill_gradientn(colours = c(green_tint, green_pale, green_fill, green_mid),
                         limits = c(0, 1), guide = "none") +
    scale_x_discrete(expand = expansion(add = .5)) +
    scale_y_discrete(expand = expansion(add = .5)) +
    coord_fixed()
  panel(p, "Predicted class", "Generating class",
        paste0("m = ", formatC(concentration, digits = 1, format = "f"), ": ", mode_label)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
}
supp_export(list(matrix_panel(.3, "position", "Position"),
                 matrix_panel(.3, "position_plus_organization", "Position + contrast"),
                 matrix_panel(1, "position", "Position"),
                 matrix_panel(1, "position_plus_organization", "Position + contrast")),
            "S04", ncol = 2, height = 151)

# Preserve the original 2.5/25/50/75/97.5 percentile box geometry explicitly.
box_panel <- function(data, value, label, title, colors, threshold = 0,
                      show_outliers = FALSE, axis_labels = waiver()) {
  values <- split(data[[value]], data$model)
  summary <- do.call(rbind, lapply(names(values), function(k) {
    q <- quantile(values[[k]], c(.025, .25, .50, .75, .975), names = FALSE)
    data.frame(model = k, ymin = q[1], lower = q[2], middle = q[3], upper = q[4], ymax = q[5])
  }))
  summary$model <- factor(summary$model, levels = levels(data$model))
  p <- ggplot(summary, aes(model, fill = model)) +
    geom_hline(yintercept = threshold, colour = muted_ink, linewidth = reference_width,
               linetype = if (threshold == 0) "solid" else "dashed") +
    geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
                 stat = "identity", colour = ink, linewidth = .35, width = .62) +
    scale_fill_manual(values = colors) + scale_y_continuous(labels = axis_labels)
  if (length(levels(data$model)) > 2) {
    p <- p + scale_x_discrete(guide = guide_axis(n.dodge = 2))
  }
  if (show_outliers) {
    outlines <- do.call(rbind, lapply(names(values), function(k) {
      limits <- summary[summary$model == k, c("ymin", "ymax")]
      z <- values[[k]]
      z <- z[z < limits$ymin | z > limits$ymax]
      data.frame(model = rep(k, length(z)), value = z)
    }))
    outlines$model <- factor(outlines$model, levels = levels(data$model))
    p <- p + geom_point(data = outlines, aes(model, value), inherit.aes = FALSE,
                        size = .8, colour = ink, alpha = .75)
  }
  panel(p, NULL, label, title) + no_legend + theme(axis.text.x = element_text(size = 8))
}

s2 <- read_data("S05_second_challenge_parents.csv")
s3 <- read_data("S05_third_challenge_parents.csv")
names(s3)[names(s3) == "force_change"] <- "delta_force_measured"
s3$model <- "H"
s5 <- rbind(s2, s3[names(s2)])
s5$model <- factor(s5$model, levels = c("D", "R", "A", "H"),
                   labels = c("Demix", "Passive", "Active", "Lookalike"))
s5_fills <- c(Demix = lilac_fill, Passive = green_pale, Active = lilac_pale, Lookalike = green_mid)
stopifnot(all(table(s5$m, s5$model) == 200))
class_panel <- function(data, title) {
  rates <- rbind(
    aggregate(spectral_positive == "True" ~ model, data = data, mean) |> setNames(c("model", "fraction")) |> transform(rule = "Spectrum"),
    aggregate(joint_positive == "True" ~ model, data = data, mean) |> setNames(c("model", "fraction")) |> transform(rule = "Spectrum + force")
  )
  rates$rule <- factor(rates$rule, levels = c("Spectrum", "Spectrum + force"))
  dodge <- position_dodge(width = .78)
  p <- ggplot(rates, aes(model, fraction, fill = rule)) +
    geom_col(width = .72, position = dodge, colour = ink, linewidth = .22) +
    geom_text(aes(label = sprintf("%.2f", fraction),
                  y = fraction + ifelse(rule == "Spectrum", .055, .16)),
              position = dodge, size = 2.8, family = figure_font, colour = ink) +
    scale_fill_manual(values = c(green_fill, lilac_fill)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(limits = c(0, 1.24), breaks = c(0, .5, 1),
                       expand = expansion(mult = c(0, .02))) +
    guides(fill = guide_legend(nrow = 1))
  panel(p, NULL, "Demixing call fraction", title) + theme(axis.text.x = element_text(size = 8))
}
s5_row <- function(concentration) {
  data <- subset(s5, m == concentration)
  title <- paste0("m = ", formatC(concentration, digits = 1, format = "f"))
  list(box_panel(data, "S", "S (per assigned s)", title, s5_fills),
       box_panel(data, "delta_force_measured", "Force change", title, s5_fills),
       class_panel(data, title))
}
s5_plots <- c(s5_row(.3), s5_row(1))
s5_widths <- c(1, 1, 1.15)
keys <- plot_grid(NULL, NULL, footer_key(s5_plots[[3]]), nrow = 1, rel_widths = s5_widths)
supp_export(s5_plots, "S05", ncol = 3, height = 131, footer = keys, widths = s5_widths)

s6 <- read_data("S06_fourth_challenge_parents.csv")
s6$model <- factor(s6$model, levels = c("D", "H"), labels = c("Demixing", "Active lookalike"))
s6_fills <- c(Demixing = lilac_fill, "Active lookalike" = green_mid)
stopifnot(all(table(s6$m, s6$model) == 200))
s6_row <- function(concentration) {
  data <- subset(s6, m == concentration)
  title <- paste0("m = ", formatC(concentration, digits = 1, format = "f"))
  list(box_panel(data, "Q", "Nonlinear ratio Q", title, s6_fills, threshold = 1.2, show_outliers = TRUE),
       box_panel(data, "force_change", "Force change (dimensionless)", title, s6_fills,
                 show_outliers = TRUE, axis_labels = scales::label_number(accuracy = .0001)))
}
supp_export(c(s6_row(.3), s6_row(1)), "S06", ncol = 2, height = 121)

# Supplementary Figure 7: actual frozen replacement test at m=1.0.
conformation <- read_result("S07_conformation.json")
ex <- conformation[c("0.3", "1.0")]
regimes <- c("m = 0.3", "m = 1.0")
accuracy <- data.frame(m = rep(regimes, each = 2),
                       channel = factor(rep(c("Position", "Position + gate"), 2),
                                        levels = c("Position", "Position + gate")),
                       pct = 100 * unlist(lapply(ex, function(e) c(e$position_only$accuracy, e$position_plus_gate$accuracy))))
dodge <- position_dodge(width = .76)
a <- panel(ggplot(accuracy, aes(m, pct, fill = channel)) +
             geom_hline(yintercept = 50, colour = muted_ink, linetype = "dotted", linewidth = reference_width) +
             geom_col(position = dodge, width = .72, colour = ink, linewidth = .23) +
             geom_text(aes(label = sprintf("%g%%", pct), y = pct + 4),
                       position = dodge, size = 2.8, family = figure_font, colour = ink) +
             scale_fill_manual(values = c(green_fill, lilac_fill)) +
             scale_y_continuous(limits = c(0, 112), breaks = c(0, 25, 50, 75, 100)),
           NULL, "Held-out accuracy (%)", "Classification")
gain <- data.frame(m = factor(regimes, levels = rev(regimes)),
                   value = 100 * vapply(ex, function(e) e$paired_gain, 0),
                   low = 100 * vapply(ex, function(e) e$paired_gain_bootstrap95[[1]], 0),
                   high = 100 * vapply(ex, function(e) e$paired_gain_bootstrap95[[2]], 0))
b <- panel(ggplot(gain, aes(value, m)) +
             geom_vline(xintercept = 0, colour = muted_ink, linewidth = reference_width) +
             geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = .16,
                           colour = forest, linewidth = .56) +
             geom_point(shape = 21, fill = lilac_fill, colour = plum, size = 2.6, stroke = mark_stroke) +
             geom_text(aes(x = high + .52, label = sprintf("%g pp", value)),
                       hjust = 0, size = 2.8, family = figure_font, colour = ink) +
             coord_cartesian(xlim = c(-1, 24)),
           "Paired accuracy gain (pp)", NULL, "Parent bootstrap 95% interval")
keys <- plot_grid(footer_key(a), NULL, nrow = 1)
supp_export(list(a, b), "S07", ncol = 2, height = 74, footer = keys)
