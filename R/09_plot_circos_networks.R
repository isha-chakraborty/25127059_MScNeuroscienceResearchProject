
# Circos representations of recurrent metal-associated candidates
#
# Generates:
#   Figure 9: all 13 recurrent copper-associated genes
#   Figure 10: the 20 most recurrent broad-iron genes
#
# Each coloured link joins occurrences of the same gene in different retained
# module–disease comparisons. Links show repeated gene–module representation;
# they do not indicate physical, regulatory or other molecular interactions.

required_packages <- c("circlize", "dplyr", "RColorBrewer", "readr", "stringr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
  library(RColorBrewer)
  library(readr)
  library(stringr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# 1. Project-relative paths and corrected recurrence inputs
# -----------------------------------------------------------------------------

find_project_root <- function() {
  candidate <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (length(list.files(candidate, pattern = "[.]Rproj$")) > 0L ||
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  normalizePath(getwd(), mustWork = TRUE)
}

project_root <- find_project_root()
figure_dir <- file.path(project_root, "outputs", "figures")
source_data_dir <- file.path(project_root, "outputs", "figure_source_data")
table_dir <- file.path(project_root, "outputs", "tables")
invisible(lapply(
  c(figure_dir, source_data_dir, table_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

copper_file <- file.path(
  table_dir, "nine_module_recurrent_copper_candidates.csv"
)
broad_iron_file <- file.path(
  table_dir, "nine_module_recurrent_broad_iron_candidates.csv"
)
occurrence_file <- file.path(
  project_root, "outputs", "derived_data", "nine_module_candidate_occurrence_long.csv"
)

if (!all(file.exists(c(copper_file, broad_iron_file, occurrence_file)))) {
  stop(
    "Corrected recurrence inputs were not found. Run ",
    "R/07_recurrence_analysis_9_modules.R first."
  )
}

copper_recurrence <- read_csv(copper_file, show_col_types = FALSE)
broad_iron_recurrence <- read_csv(broad_iron_file, show_col_types = FALSE)
occurrence <- read_csv(occurrence_file, show_col_types = FALSE) |>
  distinct(gene_symbol, module_id)

required_recurrence_columns <- c(
  "gene_symbol", "occurrences", "retained_modules"
)
if (!all(required_recurrence_columns %in% names(copper_recurrence)) ||
    !all(required_recurrence_columns %in% names(broad_iron_recurrence))) {
  stop("The recurrence files do not contain the expected Script 07 columns.")
}
if (!all(c("gene_symbol", "module_id") %in% names(occurrence))) {
  stop("The occurrence input does not contain the expected Script 07 columns.")
}

# -----------------------------------------------------------------------------
# 2. Thesis module order, labels and colours
# -----------------------------------------------------------------------------

module_order <- c(
  "AD_ERC_darkred",
  "AD_ERC_greenyellow",
  "AD_HIPPO_yellow",
  "FTLD1_red",
  "FTLD2_darkmagenta",
  "FTLD3_black",
  "FTLD3_turquoise",
  "MSA_violet",
  "PSP_pink"
)

module_display_labels <- c(
  AD_ERC_darkred = "AD ERC darkred",
  AD_ERC_greenyellow = "AD ERC greenyellow",
  AD_HIPPO_yellow = "AD HIPPO yellow",
  FTLD1_red = "FTLD1 red",
  FTLD2_darkmagenta = "FTLD2 darkmagenta",
  FTLD3_black = "FTLD3 black",
  FTLD3_turquoise = "FTLD3 turquoise",
  MSA_violet = "MSA violet",
  PSP_pink = "PSP pink"
)

module_base_colours <- c(
  AD_ERC_darkred = "#8B1A1A",
  AD_ERC_greenyellow = "#6BAF00",
  AD_HIPPO_yellow = "#E6B800",
  FTLD1_red = "#D73027",
  FTLD2_darkmagenta = "#8B006B",
  FTLD3_black = "#222222",
  FTLD3_turquoise = "#009E9A",
  MSA_violet = "#7651A8",
  PSP_pink = "#D95F9A"
)

lighten_colour <- function(colour, amount = 0.77) {
  rgb_values <- col2rgb(colour) / 255
  light_values <- rgb_values + (1 - rgb_values) * amount
  rgb(light_values[1, ], light_values[2, ], light_values[3, ])
}

module_fill_colours <- vapply(
  module_base_colours,
  lighten_colour,
  character(1)
)

unknown_modules <- setdiff(unique(occurrence$module_id), module_order)
if (length(unknown_modules) > 0L) {
  stop(
    "The occurrence file contains modules outside the corrected nine: ",
    paste(unknown_modules, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 3. Select and validate genes
# -----------------------------------------------------------------------------

rank_genes <- function(data, n = NULL) {
  ranked <- data |>
    mutate(
      gene_symbol = str_to_upper(str_trim(as.character(gene_symbol))),
      occurrences = as.integer(occurrences)
    ) |>
    arrange(desc(occurrences), gene_symbol)

  if (!is.null(n)) ranked <- slice_head(ranked, n = n)
  ranked
}

copper_selected <- rank_genes(copper_recurrence)
broad_iron_selected <- rank_genes(broad_iron_recurrence, n = 20L)

if (nrow(copper_selected) != 13L || any(copper_selected$occurrences < 2L)) {
  stop("Figure 9 requires all 13 recurrent copper-associated genes.")
}
if (nrow(broad_iron_selected) != 20L || any(broad_iron_selected$occurrences < 2L)) {
  stop("Figure 10 requires the 20 most recurrent broad-iron genes.")
}
if (anyDuplicated(copper_selected$gene_symbol) ||
    anyDuplicated(broad_iron_selected$gene_symbol)) {
  stop("Selected Circos gene lists must contain unique gene symbols.")
}

make_gene_module_long <- function(selected_genes) {
  long_data <- occurrence |>
    semi_join(selected_genes, by = "gene_symbol") |>
    filter(module_id %in% module_order) |>
    left_join(
      selected_genes |>
        select(gene_symbol, occurrences),
      by = "gene_symbol"
    ) |>
    mutate(module_id = factor(module_id, levels = module_order)) |>
    arrange(module_id, gene_symbol)

  observed_counts <- long_data |>
    count(gene_symbol, name = "observed_occurrences")
  count_check <- selected_genes |>
    select(gene_symbol, expected_occurrences = occurrences) |>
    left_join(observed_counts, by = "gene_symbol")

  if (any(count_check$expected_occurrences != count_check$observed_occurrences)) {
    print(filter(count_check, expected_occurrences != observed_occurrences))
    stop("Parsed gene–module occurrences do not match the recurrence counts.")
  }

  long_data
}

# -----------------------------------------------------------------------------
# 4. Draw a Circos plot while retaining all nine module sectors
# -----------------------------------------------------------------------------

make_circos_plot <- function(
    selected_genes,
    output_stem,
    module_cex,
    gene_cex,
    gap_size,
    link_alpha,
    gene_track_height,
    raster_size,
    pdf_size) {

  gene_module_long <- make_gene_module_long(selected_genes)

  # Every retained module is displayed, including a module with no occurrence
  # in the selected set (for example, FTLD2 darkmagenta in Figure 9).
  module_counts <- table(factor(gene_module_long$module_id, levels = module_order))
  equal_sector_width <- max(as.integer(module_counts), 1L) + 2
  sector_sizes <- tibble(
    module_id = module_order,
    xmin = 0,
    xmax = equal_sector_width
  )

  position_list <- lapply(module_order, function(current_module) {
    current <- gene_module_long |>
      filter(as.character(module_id) == current_module) |>
      arrange(gene_symbol)

    if (nrow(current) == 0L) return(current)
    current$position <- if (nrow(current) == 1L) {
      equal_sector_width / 2
    } else {
      seq(1.1, equal_sector_width - 1.1, length.out = nrow(current))
    }
    current
  })
  gene_positions <- bind_rows(position_list)

  # Gene colours follow ranked gene order so all output formats use the same
  # multicolour link palette.
  genes <- selected_genes$gene_symbol
  gene_colours <- setNames(
    colorRampPalette(brewer.pal(8, "Dark2"))(length(genes)),
    genes
  )

  draw_plot <- function() {
    par(mar = rep(0.35, 4), xpd = NA, family = "sans")
    circos.clear()
    circos.par(
      start.degree = 90,
      gap.degree = gap_size,
      track.margin = c(0.0005, 0.0005),
      cell.padding = c(0.001, 0.001, 0.001, 0.001),
      canvas.xlim = c(-1.13, 1.13),
      canvas.ylim = c(-1.13, 1.13),
      points.overflow.warning = FALSE
    )

    circos.initialize(
      factors = sector_sizes$module_id,
      xlim = cbind(sector_sizes$xmin, sector_sizes$xmax)
    )

    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = 0.105,
      bg.col = module_fill_colours[sector_sizes$module_id],
      bg.border = module_base_colours[sector_sizes$module_id],
      bg.lwd = 1.8,
      panel.fun = function(x, y) {
        current_module <- get.cell.meta.data("sector.index")
        current_xlim <- get.cell.meta.data("xlim")
        current_label <- module_display_labels[current_module]
        fitted_cex <- if (nchar(current_label) >= 18L) {
          module_cex * 0.78
        } else if (nchar(current_label) >= 15L) {
          module_cex * 0.86
        } else {
          module_cex
        }
        circos.text(
          mean(current_xlim),
          0.52,
          current_label,
          facing = "bending.inside",
          niceFacing = TRUE,
          cex = fitted_cex,
          font = 2,
          col = "#202020"
        )
      }
    )

    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = gene_track_height,
      bg.border = NA,
      bg.col = NA,
      panel.fun = function(x, y) {
        current_module <- get.cell.meta.data("sector.index")
        current_data <- gene_positions |>
          filter(as.character(module_id) == current_module)

        if (nrow(current_data) > 0L) {
          circos.text(
            x = current_data$position,
            y = 0.52,
            labels = current_data$gene_symbol,
            facing = "clockwise",
            niceFacing = TRUE,
            cex = gene_cex,
            font = 2,
            col = gene_colours[current_data$gene_symbol]
          )
        }
      }
    )

    # Pairwise links connect module occurrences of the same recurrent gene.
    for (current_gene in genes) {
      gene_data <- gene_positions |>
        filter(gene_symbol == current_gene)

      if (nrow(gene_data) >= 2L) {
        gene_pairs <- combn(seq_len(nrow(gene_data)), 2L)
        for (pair_number in seq_len(ncol(gene_pairs))) {
          first_position <- gene_data[gene_pairs[1, pair_number], ]
          second_position <- gene_data[gene_pairs[2, pair_number], ]
          circos.link(
            sector.index1 = as.character(first_position$module_id),
            point1 = first_position$position,
            sector.index2 = as.character(second_position$module_id),
            point2 = second_position$position,
            col = adjustcolor(gene_colours[current_gene], alpha.f = link_alpha),
            lwd = 1.25,
            border = NA
          )
        }
      }
    }

    circos.clear()
  }

  pdf_path <- file.path(figure_dir, paste0(output_stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(output_stem, ".png"))

  pdf(pdf_path, width = pdf_size, height = pdf_size, useDingbats = FALSE)
  draw_plot()
  dev.off()

  png(
    png_path,
    width = raster_size,
    height = raster_size,
    res = 400,
    bg = "white"
  )
  draw_plot()
  dev.off()

  source_data <- gene_positions |>
    transmute(
      gene_symbol,
      occurrences,
      recurrence_label = paste0(occurrences, "/9"),
      module_id = as.character(module_id),
      position,
      gene_colour = unname(gene_colours[gene_symbol]),
      module_colour = unname(module_base_colours[as.character(module_id)]),
      link_meaning = "Repeated occurrence of the same gene across modules"
    )

  write_csv(
    source_data,
    file.path(source_data_dir, paste0(output_stem, "_source_data.csv"))
  )

  invisible(list(pdf = pdf_path, png = png_path))
}

# -----------------------------------------------------------------------------
# 5. Export selections and create Figures 9 and 10
# -----------------------------------------------------------------------------

write_csv(
  copper_selected |>
    mutate(recurrence_label = paste0(occurrences, "/9")),
  file.path(table_dir, "Figure_9_all_13_recurrent_copper_genes.csv")
)
write_csv(
  broad_iron_selected |>
    mutate(recurrence_label = paste0(occurrences, "/9")),
  file.path(table_dir, "Figure_10_top_20_recurrent_broad_iron_genes.csv")
)

make_circos_plot(
  selected_genes = copper_selected,
  output_stem = "Figure_9_recurrent_copper_circos",
  module_cex = 1.12,
  gene_cex = 1.12,
  gap_size = 3.0,
  link_alpha = 0.54,
  gene_track_height = 0.19,
  raster_size = 6800,
  pdf_size = 16
)

make_circos_plot(
  selected_genes = broad_iron_selected,
  output_stem = "Figure_10_recurrent_broad_iron_circos",
  module_cex = 1.12,
  gene_cex = 1.02,
  gap_size = 3.0,
  link_alpha = 0.43,
  gene_track_height = 0.23,
  raster_size = 8400,
  pdf_size = 20
)

message("R/09_plot_circos_networks.R completed successfully.")
message("Figure 9 contains all 13 recurrent copper-associated genes.")
message("Figure 10 contains the 20 most recurrent broad-iron genes.")
message("Circos links represent repeated gene–module occurrence only.")
