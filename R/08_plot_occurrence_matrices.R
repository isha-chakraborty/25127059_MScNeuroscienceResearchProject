
# Occurrence matrices for recurrent metal-associated candidates
#
# Generates:
#   Figure 7: all 13 recurrent copper-associated genes
#   Figure 8: the 25 highest-recurrence broad-iron genes
#   Supplementary Figure S1: all 77 recurrent broad-iron genes
#
# The plots show gene occurrence across the corrected nine metal-enriched
# module–disease comparisons. An occupied cell indicates that a gene was assigned
# to that module; it does not represent a molecular interaction or enrichment
# strength.

required_packages <- c("dplyr", "ggplot2", "readr", "stringr", "tibble", "tidyr")
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
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 1. Project-relative paths
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

output_dirs <- list(
  figures = file.path(project_root, "outputs", "figures"),
  source_data = file.path(project_root, "outputs", "figure_source_data"),
  tables = file.path(project_root, "outputs", "tables")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

recurrence_file <- file.path(
  project_root, "outputs", "tables", "nine_module_all_recurrent_candidates.csv"
)
occurrence_file <- file.path(
  project_root, "outputs", "derived_data", "nine_module_candidate_occurrence_long.csv"
)

if (!file.exists(recurrence_file) || !file.exists(occurrence_file)) {
  stop(
    "Corrected recurrence inputs were not found. Run ",
    "R/07_recurrence_analysis_9_modules.R first."
  )
}

recurrence <- read_csv(recurrence_file, show_col_types = FALSE)
occurrence <- read_csv(occurrence_file, show_col_types = FALSE)

required_recurrence_columns <- c(
  "gene_symbol", "occurrences", "retained_modules", "copper",
  "broad_iron", "core_iron"
)
required_occurrence_columns <- c("gene_symbol", "module_id")

if (!all(required_recurrence_columns %in% names(recurrence))) {
  stop("The recurrence table does not have the expected Script 07 columns.")
}
if (!all(required_occurrence_columns %in% names(occurrence))) {
  stop("The occurrence table does not have the expected Script 07 columns.")
}

# -----------------------------------------------------------------------------
# 2. Module labels, order and original WGCNA colours
# -----------------------------------------------------------------------------

module_key <- tribble(
  ~module_order, ~module_id,             ~module_label,          ~disease_group,       ~module_colour,
  1L,            "AD_ERC_greenyellow",  "ERC\ngreenyellow",    "Alzheimer's disease", "#6BAF00",
  2L,            "AD_ERC_darkred",      "ERC\ndarkred",        "Alzheimer's disease", "#8B1A1A",
  3L,            "AD_HIPPO_yellow",     "HIPPO\nyellow",       "Alzheimer's disease", "#E6B800",
  4L,            "FTLD1_red",           "FTLD1\nred",          "Frontotemporal lobar degeneration", "#D73027",
  5L,            "FTLD2_darkmagenta",   "FTLD2\ndarkmagenta",  "Frontotemporal lobar degeneration", "#8B006B",
  6L,            "FTLD3_black",         "FTLD3\nblack",        "Frontotemporal lobar degeneration", "#222222",
  7L,            "FTLD3_turquoise",     "FTLD3\nturquoise",    "Frontotemporal lobar degeneration", "#009E9A",
  8L,            "MSA_violet",          "MSA\nviolet",         "Movement disorders", "#7651A8",
  9L,            "PSP_pink",            "PSP\npink",           "Movement disorders", "#D95F9A"
) |>
  mutate(
    module_id = factor(module_id, levels = module_id),
    module_label = factor(module_label, levels = module_label),
    disease_group = factor(
      disease_group,
      levels = c(
        "Alzheimer's disease",
        "Frontotemporal lobar degeneration",
        "Movement disorders"
      )
    )
  )

module_colours <- setNames(
  module_key$module_colour,
  as.character(module_key$module_label)
)

unknown_modules <- setdiff(unique(occurrence$module_id), as.character(module_key$module_id))
if (length(unknown_modules) > 0L) {
  stop(
    "The occurrence input contains modules outside the corrected nine: ",
    paste(unknown_modules, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 3. Construct complete gene-by-module plotting grids
# -----------------------------------------------------------------------------

make_matrix_data <- function(selected_genes) {
  gene_order <- selected_genes |>
    arrange(desc(occurrences), gene_symbol) |>
    mutate(gene_label = paste0(gene_symbol, "  (", occurrences, "/9)"))

  plot_grid <- crossing(
    gene_symbol = gene_order$gene_symbol,
    module_id = as.character(module_key$module_id)
  ) |>
    left_join(
      occurrence |>
        distinct(gene_symbol, module_id) |>
        mutate(present = TRUE),
      by = c("gene_symbol", "module_id")
    ) |>
    mutate(present = replace_na(present, FALSE)) |>
    left_join(
      module_key |>
        mutate(module_id = as.character(module_id)),
      by = "module_id"
    ) |>
    left_join(
      gene_order |>
        select(
          gene_symbol, gene_label, occurrences, copper,
          broad_iron, core_iron
        ),
      by = "gene_symbol"
    ) |>
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_order$gene_label)),
      module_label = factor(
        as.character(module_label),
        levels = levels(module_key$module_label)
      ),
      disease_group = factor(
        disease_group,
        levels = levels(module_key$disease_group)
      )
    )

  count_check <- plot_grid |>
    filter(present) |>
    count(gene_symbol, name = "observed_occurrences") |>
    right_join(
      gene_order |>
        select(gene_symbol, expected_occurrences = occurrences),
      by = "gene_symbol"
    )

  if (any(count_check$observed_occurrences != count_check$expected_occurrences)) {
    print(filter(count_check, observed_occurrences != expected_occurrences))
    stop("At least one plotted gene has an incorrect module-occurrence count.")
  }

  plot_grid
}

# -----------------------------------------------------------------------------
# 4. Shared plot design
# -----------------------------------------------------------------------------

matrix_theme <- function(base_size = 12, y_size = 10.5, x_size = 10.5) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(
        size = base_size,
        colour = "#222222",
        margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        size = base_size,
        colour = "#222222",
        margin = margin(r = 12)
      ),
      axis.text.x = element_text(
        size = x_size,
        angle = 45,
        hjust = 1,
        vjust = 1,
        colour = "#222222",
        lineheight = 0.9
      ),
      axis.text.y = element_text(
        size = y_size,
        face = "italic",
        colour = "#222222",
        margin = margin(r = 6)
      ),
      strip.text = element_text(
        size = 12.5,
        face = "bold",
        colour = "#1F1F1F",
        margin = margin(6, 4, 6, 4)
      ),
      strip.background = element_blank(),
      panel.spacing.x = grid::unit(0.35, "cm"),
      legend.position = "top",
      legend.justification = "right",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10.5),
      plot.margin = margin(12, 18, 16, 12)
    )
}

base_matrix_layers <- function(plot_data, tile_linewidth = 0.35) {
  list(
    geom_tile(
      fill = "#F6F7F8",
      colour = "white",
      linewidth = tile_linewidth,
      width = 0.96,
      height = 0.92
    ),
    facet_grid(
      . ~ disease_group,
      scales = "free_x",
      space = "free_x"
    ),
    scale_fill_manual(values = module_colours, guide = "none"),
    coord_cartesian(clip = "off")
  )
}

make_copper_plot <- function(selected_genes) {
  plot_data <- make_matrix_data(selected_genes)

  ggplot(plot_data, aes(module_label, gene_label)) +
    base_matrix_layers(plot_data) +
    geom_point(
      data = filter(plot_data, present),
      aes(fill = module_label),
      shape = 21,
      size = 4.0,
      stroke = 0.65,
      colour = "#3F3F3F"
    ) +
    labs(
      x = "Retained metal-enriched co-methylation module",
      y = "Recurrent copper-associated gene (occurrence count across nine modules)"
    ) +
    matrix_theme(base_size = 12, y_size = 11, x_size = 10.5)
}

make_iron_plot <- function(selected_genes, compact = FALSE) {
  plot_data <- make_matrix_data(selected_genes) |>
    mutate(
      iron_class = if_else(
        core_iron,
        "Core iron (also broad iron)",
        "Broad iron only"
      ),
      iron_class = factor(
        iron_class,
        levels = c("Broad iron only", "Core iron (also broad iron)")
      )
    )

  ggplot(plot_data, aes(module_label, gene_label)) +
    base_matrix_layers(
      plot_data,
      tile_linewidth = if (compact) 0.22 else 0.35
    ) +
    geom_point(
      data = filter(plot_data, present),
      aes(fill = module_label, shape = iron_class),
      size = if (compact) 2.8 else 3.8,
      stroke = 0.6,
      colour = "#3F3F3F"
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "Broad iron only" = 22,
        "Core iron (also broad iron)" = 21
      ),
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          fill = "#666666",
          colour = "#3F3F3F",
          size = 4
        )
      )
    ) +
    labs(
      x = "Retained metal-enriched co-methylation module",
      y = "Recurrent broad-iron candidate gene (occurrence count across nine modules)"
    ) +
    matrix_theme(
      base_size = if (compact) 11 else 12,
      y_size = if (compact) 7.7 else 10.2,
      x_size = if (compact) 9.5 else 10.5
    )
}

# -----------------------------------------------------------------------------
# 5. Select corrected recurrent sets and generate plots
# -----------------------------------------------------------------------------

copper_genes <- recurrence |>
  filter(copper) |>
  arrange(desc(occurrences), gene_symbol)

broad_iron_genes <- recurrence |>
  filter(broad_iron) |>
  arrange(desc(occurrences), gene_symbol)

if (nrow(copper_genes) != 13L) {
  stop("Figure 7 requires exactly 13 recurrent copper-associated genes.")
}
if (nrow(broad_iron_genes) != 77L) {
  stop("The complete broad-iron matrix requires exactly 77 recurrent genes.")
}

# Ties are resolved alphabetically after ranking by occurrence count.
prioritised_broad_iron_genes <- broad_iron_genes |>
  slice_head(n = 25L)

figure_7_data <- make_matrix_data(copper_genes)
figure_8_data <- make_matrix_data(prioritised_broad_iron_genes)
supplementary_s1_data <- make_matrix_data(broad_iron_genes)

figure_7 <- make_copper_plot(copper_genes)
figure_8 <- make_iron_plot(prioritised_broad_iron_genes, compact = FALSE)
supplementary_figure_s1 <- make_iron_plot(broad_iron_genes, compact = TRUE)

# -----------------------------------------------------------------------------
# 6. Save high-resolution figures and their source data
# -----------------------------------------------------------------------------

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white"
  )
}

save_plot_set(
  figure_7,
  "Figure_7_recurrent_copper_occurrence_matrix",
  width = 12.5,
  height = 8.2
)
save_plot_set(
  figure_8,
  "Figure_8_prioritised_broad_iron_occurrence_matrix",
  width = 13.2,
  height = 12.2
)
save_plot_set(
  supplementary_figure_s1,
  "Supplementary_Figure_S1_complete_broad_iron_occurrence_matrix",
  width = 13.2,
  height = 25.5
)

write_plot_source <- function(plot_data, path, include_iron_class = FALSE) {
  source_table <- plot_data |>
    filter(present) |>
    transmute(
      gene_symbol,
      occurrences,
      recurrence_label = paste0(occurrences, "/9"),
      module_id,
      module_order,
      module_label = as.character(module_label),
      disease_group = as.character(disease_group),
      copper,
      broad_iron,
      core_iron
    )

  if (include_iron_class) {
    source_table <- source_table |>
      mutate(
        iron_class = if_else(
          core_iron,
          "Core iron (also broad iron)",
          "Broad iron only"
        )
      )
  }

  write_csv(source_table, path)
}

write_plot_source(
  figure_7_data,
  file.path(output_dirs$source_data, "Figure_7_source_data.csv")
)
write_plot_source(
  figure_8_data,
  file.path(output_dirs$source_data, "Figure_8_source_data.csv"),
  include_iron_class = TRUE
)
write_plot_source(
  supplementary_s1_data,
  file.path(output_dirs$source_data, "Supplementary_Figure_S1_source_data.csv"),
  include_iron_class = TRUE
)

write_csv(
  copper_genes |>
    select(
      gene_symbol, occurrences, recurrence_label,
      retained_modules, copper, broad_iron, core_iron
    ),
  file.path(output_dirs$tables, "Figure_7_all_13_recurrent_copper_genes.csv")
)
write_csv(
  prioritised_broad_iron_genes |>
    select(
      gene_symbol, occurrences, recurrence_label,
      retained_modules, broad_iron, core_iron
    ),
  file.path(output_dirs$tables, "Figure_8_prioritised_25_broad_iron_genes.csv")
)
write_csv(
  broad_iron_genes |>
    select(
      gene_symbol, occurrences, recurrence_label,
      retained_modules, broad_iron, core_iron
    ),
  file.path(output_dirs$tables, "Supplementary_Figure_S1_all_77_broad_iron_genes.csv")
)

message("R/08_plot_occurrence_matrices.R completed successfully.")
message("Figure 7: 13 recurrent copper-associated genes.")
message("Figure 8: 25 prioritised recurrent broad-iron genes.")
message("Supplementary Figure S1: all 77 recurrent broad-iron genes.")
