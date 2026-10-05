#!/usr/bin/env Rscript

# High-resolution DAG for the Parker B. Francis sepsis-ARDS pollution project.

suppressPackageStartupMessages({
  library(grid)
  library(grDevices)
})

out_dir <- file.path("output", "figures", "pbf_study_dag")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- data.frame(
  id = c("conf", "pollution", "host", "microbe", "sepsis", "ards", "outcomes", "culture"),
  label = c(
    "Patient, social, seasonal,\ngeographic, and site factors",
    "Air pollution exposure\nFine particulate matter,\nnitrogen dioxide, ozone",
    "Airway injury and\nimpaired host defense",
    "Respiratory\npathogen ecology",
    "Sepsis",
    "Incident sepsis-associated\nARDS",
    "ARDS severity and\nclinical outcomes",
    "Respiratory\nculture detection"
  ),
  x = c(0.49, 0.13, 0.35, 0.35, 0.13, 0.59, 0.86, 0.59),
  y = c(0.86, 0.63, 0.70, 0.45, 0.25, 0.58, 0.58, 0.24),
  w = c(0.29, 0.23, 0.21, 0.21, 0.13, 0.23, 0.21, 0.20),
  h = c(0.12, 0.15, 0.12, 0.12, 0.09, 0.12, 0.12, 0.12),
  fill = c("#EEF2F6", "#FFF4D6", "#E7F0FA", "#E8F5EC", "#F7E8EA", "#FDECEC", "#F1EBF8", "#F3F3F3"),
  stroke = c("#2F3A45", "#2F3A45", "#2F3A45", "#2F3A45", "#2F3A45", "#9D2F2F", "#2F3A45", "#2F3A45"),
  stringsAsFactors = FALSE
)

get_node <- function(id) nodes[nodes$id == id, ][1, ]

point_on <- function(node, side) {
  if (side == "right") return(c(node$x + node$w / 2, node$y))
  if (side == "left") return(c(node$x - node$w / 2, node$y))
  if (side == "top") return(c(node$x, node$y + node$h / 2))
  if (side == "bottom") return(c(node$x, node$y - node$h / 2))
  stop("Unknown side: ", side)
}

draw_arrow <- function(from, to, from_side = "right", to_side = "left",
                       curvature = 0, col = "#3E4A56", lwd = 1.6) {
  a <- point_on(get_node(from), from_side)
  b <- point_on(get_node(to), to_side)
  grid.curve(
    x1 = unit(a[1], "npc"),
    y1 = unit(a[2], "npc"),
    x2 = unit(b[1], "npc"),
    y2 = unit(b[2], "npc"),
    curvature = curvature,
    angle = 90,
    square = FALSE,
    inflect = FALSE,
    arrow = arrow(type = "closed", length = unit(0.12, "inches")),
    gp = gpar(col = col, lwd = lwd, lineend = "round")
  )
}

draw_node <- function(node) {
  grid.roundrect(
    x = unit(node$x, "npc"),
    y = unit(node$y, "npc"),
    width = unit(node$w, "npc"),
    height = unit(node$h, "npc"),
    r = unit(0.045, "snpc"),
    gp = gpar(fill = node$fill, col = node$stroke, lwd = 1.25)
  )
  grid.text(
    node$label,
    x = unit(node$x, "npc"),
    y = unit(node$y, "npc"),
    gp = gpar(col = "#111111", fontsize = 12.5, fontfamily = "Helvetica", lineheight = 0.9)
  )
}

draw_dag <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  grid.text(
    "Confounding context",
    x = unit(0.25, "npc"),
    y = unit(0.965, "npc"),
    gp = gpar(col = "#5B636A", fontsize = 10.5, fontfamily = "Helvetica")
  )
  grid.text(
    "Primary clinical pathway",
    x = unit(0.70, "npc"),
    y = unit(0.78, "npc"),
    gp = gpar(col = "#5B636A", fontsize = 10.5, fontfamily = "Helvetica")
  )
  grid.text(
    "Clinical sampling process",
    x = unit(0.36, "npc"),
    y = unit(0.07, "npc"),
    gp = gpar(col = "#5B636A", fontsize = 10.5, fontfamily = "Helvetica")
  )

  # Primary biologic pathway
  draw_arrow("pollution", "host", lwd = 1.8)
  draw_arrow("pollution", "microbe", lwd = 1.8)
  draw_arrow("host", "ards", curvature = -0.05, lwd = 1.8)
  draw_arrow("microbe", "ards", curvature = 0.04, lwd = 1.8)
  draw_arrow("sepsis", "ards", curvature = -0.18, lwd = 1.8)
  draw_arrow("ards", "outcomes", lwd = 1.8)
  draw_arrow("microbe", "outcomes", curvature = -0.22, lwd = 1.35)

  # Culture observation process
  draw_arrow("microbe", "culture", curvature = -0.08, col = "#6B7280", lwd = 1.35)
  draw_arrow("sepsis", "culture", curvature = 0.06, col = "#6B7280", lwd = 1.35)
  draw_arrow("ards", "culture", from_side = "bottom", to_side = "top", col = "#6B7280", lwd = 1.35)

  # Shared confounding context
  draw_arrow("conf", "pollution", from_side = "bottom", to_side = "top", curvature = 0.16, col = "#7A828A", lwd = 1.1)
  draw_arrow("conf", "sepsis", from_side = "bottom", to_side = "top", curvature = -0.10, col = "#7A828A", lwd = 1.1)
  draw_arrow("conf", "ards", from_side = "bottom", to_side = "top", curvature = 0.02, col = "#7A828A", lwd = 1.1)
  draw_arrow("conf", "outcomes", from_side = "bottom", to_side = "top", curvature = -0.12, col = "#7A828A", lwd = 1.1)
  draw_arrow("conf", "culture", from_side = "bottom", to_side = "top", curvature = 0.12, col = "#7A828A", lwd = 1.1)

  invisible(apply(nodes, 1, function(row) draw_node(as.list(row))))
}

png_path <- file.path(out_dir, "pbf_air_pollution_sepsis_ards_dag.png")
pdf_path <- file.path(out_dir, "pbf_air_pollution_sepsis_ards_dag.pdf")
svg_path <- file.path(out_dir, "pbf_air_pollution_sepsis_ards_dag.svg")

png(png_path, width = 13.2, height = 6.3, units = "in", res = 600, bg = "white")
draw_dag()
dev.off()

pdf(pdf_path, width = 13.2, height = 6.3, bg = "white", useDingbats = FALSE)
draw_dag()
dev.off()

svg(svg_path, width = 13.2, height = 6.3, bg = "white")
draw_dag()
dev.off()

message("Wrote: ", normalizePath(png_path, winslash = "/"))
message("Wrote: ", normalizePath(pdf_path, winslash = "/"))
message("Wrote: ", normalizePath(svg_path, winslash = "/"))
