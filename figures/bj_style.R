# Shared typography, inks and fills for the Biophysical Journal figures.
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(svglite)
})
cowplot::set_null_device("agg")

# Liberation Sans and Arimo share Arial's character metrics.
font_candidates <- c("Arial", "Liberation Sans", "Arimo")
available_fonts <- systemfonts::system_fonts()$family
figure_font <- font_candidates[font_candidates %in% available_fonts][1]
if (is.na(figure_font)) {
  stop("Install Arial, Liberation Sans, or Arimo before rendering.", call. = FALSE)
}
if (!capabilities("cairo")) {
  stop("R must have Cairo support to render the figure PDFs.", call. = FALSE)
}

render_args <- commandArgs(trailingOnly = TRUE)
if (length(render_args) > 1L) {
  stop("Pass at most one argument (the figure output directory).", call. = FALSE)
}
output_dir <- if (length(render_args)) render_args[1] else {
  file.path(dirname(here), "build", "figures")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
  stop(paste("Could not create figure output directory", output_dir), call. = FALSE)
}
data_dir <- Sys.getenv("REPRO_DATA_DIR", unset = file.path(here, "data"))
if (!dir.exists(data_dir)) {
  stop(paste("Figure data directory does not exist", data_dir), call. = FALSE)
}
message("Figure font: ", figure_font)
ink <- "#202D2B"
muted_ink <- "#546661"
neutral_fill <- "#EDF1F4"

forest <- "#214A3E"
green_ink <- "#427361"
green_mid <- "#79B69F"
green_fill <- "#A8D6BF"
green_pale <- "#D3EDDF"
green_tint <- "#EBF7F0"

plum <- "#694B79"
lilac_fill <- "#CBB3E6"
lilac_pale <- "#E3D3F1"
lilac_tint <- "#F5EFF9"
pink_fill <- "#F1C9D7"
pink_tint <- "#FAEAF0"

trace_width <- .74
reference_width <- .34
mark_stroke <- .32
comparison_lines <- c("solid", "dashed", "dotted")

th <- function() {
  theme_classic(base_size = 10, base_family = figure_font) + theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink, size = 9),
    axis.title = element_text(size = 10),
    axis.line = element_line(colour = ink, linewidth = .32),
    axis.ticks = element_line(colour = ink, linewidth = .31),
    axis.ticks.length = grid::unit(1.25, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.position = "bottom",
    legend.key.width = grid::unit(6, "mm"),
    legend.key.height = grid::unit(3.6, "mm"),
    legend.key.spacing.x = grid::unit(4, "pt"),
    legend.margin = margin(1, 0, 1, 0),
    legend.box.spacing = grid::unit(4, "pt"),
    plot.title = element_text(size = 10, face = "bold", colour = forest,
                              margin = margin(t = 7, b = 5)),
    plot.margin = margin(11, 8, 8, 8),
    plot.background = element_rect(fill = "white", colour = NA))
}

tag_panel <- function(panel, label) {
  tag <- grid::textGrob(label, x = grid::unit(5, "pt"),
                        y = grid::unit(1, "npc") - grid::unit(5, "pt"),
                        just = c("left", "top"),
                        gp = grid::gpar(fontsize = 12, col = ink, fontfamily = figure_font,
                                       fontface = "bold"))
  ggdraw(panel) + draw_grob(tag)
}

aligned_row <- function(panels, labels = NULL, widths = rep(1, length(panels))) {
  panels <- align_plots(plotlist = panels, align = "h", axis = "tb")
  if (!is.null(labels)) panels <- Map(tag_panel, panels, labels)
  plot_grid(plotlist = panels, nrow = 1, rel_widths = widths)
}

labelled_grid <- function(panels, ncol, widths = rep(1, ncol), labels = letters[seq_along(panels)]) {
  rows <- split(seq_along(panels), ceiling(seq_along(panels) / ncol))
  for (column in seq_len(ncol)) {
    indices <- seq(column, length(panels), by = ncol)
    panels[indices] <- align_plots(plotlist = panels[indices], align = "v", axis = "lr")
  }
  plot_grid(plotlist = lapply(rows, function(indices) {
    aligned_row(panels[indices], labels[indices], widths)
  }), ncol = 1)
}

no_legend <- theme(legend.position = "none")

footer_key <- function(plot, title = NULL) {
  key <- get_legend(plot, legend = "bottom")
  if (is.null(title)) return(ggdraw(key))
  title <- ggdraw() + draw_label(title, size = 8.5, fontfamily = figure_font,
                               fontface = "bold", colour = forest)
  plot_grid(title, ggdraw(key), ncol = 1, rel_heights = c(.65, 1))
}

export <- function(g, name, width = 180, height = 82, directory = output_dir) {
  stem <- file.path(directory, name)
  ggsave(paste0(stem, ".pdf"), g, width = width, height = height,
         units = "mm", device = cairo_pdf, family = figure_font, bg = "white")
  ggsave(paste0(stem, ".svg"), g, width = width, height = height,
         units = "mm", device = svglite::svglite, bg = "white")
  ggsave(paste0(stem, ".png"), g, width = width, height = height,
         units = "mm", device = ragg::agg_png, dpi = 600, bg = "white")
}
