
# Candidate prioritisation using module membership and gene significance
#
# Reproduces:
#   - Figure 6 (four representative retained metal-enriched modules)
#   - Supplementary Figure S2 (the other five retained modules)
#   - Spearman correlation and candidate-ranking tables
#
# Run from the project root, preferably after R/02_curate_metal_gene_sets.R
# and R/04_metal_gene_set_enrichment.R.

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "patchwork", "purrr",
  "readr", "stringr", "tibble"
)

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
  library(ggrepel)
  library(patchwork)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

set.seed(2026)

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
  tables = file.path(project_root, "outputs", "tables"),
  source_data = file.path(project_root, "outputs", "figure_source_data"),
  derived_data = file.path(project_root, "outputs", "derived_data")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# -----------------------------------------------------------------------------
# 2. Corrected nine-module specification
# -----------------------------------------------------------------------------

module_specification <- tribble(
  ~module_order, ~module_id,              ~object_name,       ~module,       ~mm_column,       ~gs_column,        ~metal_set,              ~plot_colour, ~main_panel, ~panel,
  1L,            "AD_ERC_greenyellow",   "geneInfo_ERC",     "greenyellow", "MM.greenyellow", "GS.sample_group", "Copper",                "#78A800",    TRUE,        "A",
  2L,            "AD_ERC_darkred",       "geneInfo_ERC",     "darkred",     "MM.darkred",     "GS.sample_group", "Broad iron",            "#8B1A1A",    FALSE,       "A",
  3L,            "AD_HIPPO_yellow",      "geneInfo_HIPPO",   "yellow",      "MM.yellow",      "GS.sample_group", "Core iron",             "#D8A900",    FALSE,       "B",
  4L,            "FTLD1_red",            "geneInfo_FTLD1",   "red",         "MM.red",         "GS.FTLD",         "Copper",                "#D73027",    TRUE,        "B",
  5L,            "FTLD2_darkmagenta",    "geneInfo_FTLD2",   "darkmagenta", "MMdarkmagenta",  "GS.sample_group", "Broad and core iron",   "#9E006E",    TRUE,        "C",
  6L,            "FTLD3_black",          "geneInfo_FTLD3",   "black",       "MM.black",       "GS.sample_group", "Broad and core iron",   "#2B2B2B",    FALSE,       "C",
  7L,            "FTLD3_turquoise",      "geneInfo_FTLD3",   "turquoise",   "MM.turquoise",   "GS.sample_group", "Broad and core iron",   "#109C95",    TRUE,        "D",
  8L,            "MSA_violet",           "MSAgeneInfo0K",    "violet",      "MM.violet",      "GS.MSA",          "Core iron",             "#7655A5",    FALSE,       "D",
  9L,            "PSP_pink",             "PSPgeneInfo0K",    "pink",        "MM.pink",        "GS.PSP",          "Core iron",             "#D65A91",    FALSE,       "E"
)

write_csv(
  module_specification,
  file.path(output_dirs$tables, "corrected_nine_module_MM_GS_specification.csv")
)

# -----------------------------------------------------------------------------
# 3. Load curated metal-gene sets produced by Script 02
# -----------------------------------------------------------------------------

standardise_symbols <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    str_to_upper() |>
    (\(z) z[!is.na(z) & z != ""])() |>
    unique() |>
    sort()
}

gene_set_rds_candidates <- c(
  file.path(project_root, "outputs", "derived_data", "curated_metal_gene_sets.rds"),
  file.path(project_root, "data", "derived", "curated_metal_gene_sets.rds"),
  file.path(project_root, "gene_lists", "curated_metal_gene_sets.rds")
)
gene_set_rds <- gene_set_rds_candidates[file.exists(gene_set_rds_candidates)][1]

if (is.na(gene_set_rds)) {
  stop(
    "The curated gene-set object was not found. Run ",
    "R/02_curate_metal_gene_sets.R first. The expected RDS object is a named ",
    "list containing copper, broad_iron and core_iron vectors."
  )
}

metal_gene_sets <- readRDS(gene_set_rds)
required_sets <- c("copper", "broad_iron", "core_iron")
if (!all(required_sets %in% names(metal_gene_sets))) {
  stop("The curated gene-set RDS must contain: ", paste(required_sets, collapse = ", "))
}

metal_gene_sets <- lapply(metal_gene_sets[required_sets], standardise_symbols)

expected_counts <- c(copper = 47L, broad_iron = 199L, core_iron = 136L)
observed_counts <- vapply(metal_gene_sets, length, integer(1))
if (!identical(unname(observed_counts), unname(expected_counts))) {
  stop(
    "Curated gene-set totals differ from the established analysis. Observed: ",
    paste(names(observed_counts), observed_counts, sep = " = ", collapse = "; ")
  )
}
if (!all(metal_gene_sets$core_iron %in% metal_gene_sets$broad_iron)) {
  stop("The core-iron set must be a subset of the broad-iron set.")
}

# -----------------------------------------------------------------------------
# 4. Load the published network gene-information objects
# -----------------------------------------------------------------------------

# A neutral combined RDS is preferred for a public repository. Its elements must
# be named as in `object_name` above. For compatibility with the completed local
# project, the script can also locate RData files by their contained object names.
combined_network_rds <- file.path(
  project_root, "data", "published_results", "network_gene_information.rds"
)
public_metrics_file <- file.path(
  project_root, "data", "derived", "nine_module_MM_GS_annotated_records.csv"
)

if (!file.exists(public_metrics_file) && file.exists(combined_network_rds)) {
  network_objects <- readRDS(combined_network_rds)
} else if (!file.exists(public_metrics_file)) {
  search_directories <- c(
    project_root,
    file.path(project_root, "data_processed"),
    file.path(project_root, "data", "published_results")
  )
  search_directories <- search_directories[dir.exists(search_directories)]
  rdata_files <- unique(unlist(lapply(
    search_directories,
    list.files,
    pattern = "[.](RData|rda)$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )))

  required_objects <- unique(module_specification$object_name)
  network_objects <- list()

  for (rdata_file in rdata_files) {
    temporary_environment <- new.env(parent = emptyenv())
    loaded_names <- tryCatch(
      load(rdata_file, envir = temporary_environment),
      error = function(e) character()
    )
    matches <- intersect(required_objects, loaded_names)
    for (object_name in matches) {
      if (is.null(network_objects[[object_name]])) {
        network_objects[[object_name]] <- get(object_name, envir = temporary_environment)
      }
    }
  }
}

if (!file.exists(public_metrics_file)) {
  missing_objects <- setdiff(unique(module_specification$object_name), names(network_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Missing published network gene-information objects: ",
      paste(missing_objects, collapse = ", "), ". Place the source RData files in ",
      "the project root/data_processed directory, create ",
      "data/published_results/network_gene_information.rds, or provide the ",
      "documented derived input data/derived/nine_module_MM_GS_annotated_records.csv."
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Extract absolute MM and GS values and annotate metal-set membership
# -----------------------------------------------------------------------------

candidate_rule <- function(gene_symbol, metal_set) {
  if (metal_set == "Copper") {
    gene_symbol %in% metal_gene_sets$copper
  } else if (metal_set == "Broad iron") {
    gene_symbol %in% metal_gene_sets$broad_iron
  } else if (metal_set == "Core iron") {
    gene_symbol %in% metal_gene_sets$core_iron
  } else {
    gene_symbol %in% metal_gene_sets$broad_iron
  }
}

extract_module_metrics <- function(specification_row) {
  specification_row <- as.list(specification_row)
  network_data <- network_objects[[specification_row$object_name]]

  required_columns <- c(
    "geneSymbol", "moduleColor", specification_row$mm_column,
    specification_row$gs_column
  )
  absent_columns <- setdiff(required_columns, colnames(network_data))
  if (length(absent_columns) > 0L) {
    stop(
      specification_row$module_id, " is missing columns: ",
      paste(absent_columns, collapse = ", ")
    )
  }

  extracted <- network_data |>
    transmute(
      gene_symbol = standardise_symbols(geneSymbol)[match(
        str_to_upper(str_trim(as.character(geneSymbol))),
        standardise_symbols(geneSymbol)
      )],
      module_colour = str_to_lower(str_trim(as.character(moduleColor))),
      module_membership = abs(suppressWarnings(as.numeric(.data[[specification_row$mm_column]]))),
      gene_significance = abs(suppressWarnings(as.numeric(.data[[specification_row$gs_column]])))
    ) |>
    filter(
      module_colour == specification_row$module,
      !is.na(gene_symbol), gene_symbol != "",
      is.finite(module_membership), is.finite(gene_significance)
    ) |>
    mutate(
      module_id = specification_row$module_id,
      module_order = specification_row$module_order,
      metal_enrichment = specification_row$metal_set,
      plot_colour = specification_row$plot_colour,
      main_panel = specification_row$main_panel,
      panel = specification_row$panel,
      copper_gene = gene_symbol %in% metal_gene_sets$copper,
      broad_iron_gene = gene_symbol %in% metal_gene_sets$broad_iron,
      core_iron_gene = gene_symbol %in% metal_gene_sets$core_iron,
      highlighted_candidate = candidate_rule(gene_symbol, specification_row$metal_set),
      priority_score = module_membership * gene_significance,
      record_id = row_number(),
      .before = 1
    )

  if (nrow(extracted) == 0L) {
    stop("No valid MM–GS records were found for ", specification_row$module_id, ".")
  }
  extracted
}

if (file.exists(public_metrics_file)) {
  message("Using the documented derived MM–GS input: ", public_metrics_file)
  module_metrics <- read_csv(public_metrics_file, show_col_types = FALSE)
  required_metric_columns <- c(
    "record_id", "module_id", "module_order", "gene_symbol",
    "module_membership", "gene_significance", "metal_enrichment",
    "plot_colour", "main_panel", "panel", "copper_gene",
    "broad_iron_gene", "core_iron_gene", "highlighted_candidate",
    "priority_score"
  )
  missing_metric_columns <- setdiff(required_metric_columns, names(module_metrics))
  if (length(missing_metric_columns) > 0L) {
    stop(
      "The derived MM–GS input is missing columns: ",
      paste(missing_metric_columns, collapse = ", ")
    )
  }
  if (!setequal(unique(module_metrics$module_id), module_specification$module_id)) {
    stop("The derived MM–GS input does not contain exactly the corrected nine modules.")
  }
} else {
  module_metrics <- map_dfr(
    seq_len(nrow(module_specification)),
    ~ extract_module_metrics(module_specification[.x, , drop = FALSE])
  )
}

# Preserve all annotated CpG–gene records for the scatterplots and correlations.
# A gene can therefore appear more than once when multiple CpGs map to it.
write_csv(
  module_metrics,
  file.path(output_dirs$derived_data, "nine_module_MM_GS_annotated_records.csv")
)

# -----------------------------------------------------------------------------
# 6. Spearman correlations and candidate ranking
# -----------------------------------------------------------------------------

format_p_value <- function(p_value) {
  if (is.na(p_value)) return(NA_character_)
  if (p_value < 2.2e-16) return("< 2.2 x 10^-16")
  format.pval(p_value, digits = 2, eps = 2.2e-16)
}

correlation_results <- module_metrics |>
  group_by(module_order, module_id, metal_enrichment) |>
  group_modify(~ {
    test <- suppressWarnings(cor.test(
      .x$module_membership,
      .x$gene_significance,
      method = "spearman",
      exact = FALSE
    ))
    tibble(
      annotated_records = nrow(.x),
      unique_genes = n_distinct(.x$gene_symbol),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      p_value_display = format_p_value(test$p.value)
    )
  }) |>
  ungroup() |>
  arrange(module_order)

write_csv(
  correlation_results,
  file.path(output_dirs$tables, "MM_GS_spearman_correlations_nine_modules.csv")
)

# Candidate ranking uses one representative record per gene: the record with the
# highest |MM| x |GS| score. This prevents genes with several annotated CpGs from
# receiving several ranks while retaining their strongest module-level signal.
candidate_rankings <- module_metrics |>
  filter(highlighted_candidate) |>
  group_by(module_order, module_id, metal_enrichment, gene_symbol) |>
  slice_max(priority_score, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  group_by(module_order, module_id) |>
  arrange(desc(priority_score), desc(module_membership), desc(gene_significance), gene_symbol, .by_group = TRUE) |>
  mutate(candidate_rank = row_number()) |>
  ungroup() |>
  select(
    module_order, module_id, metal_enrichment, candidate_rank, gene_symbol,
    module_membership, gene_significance, priority_score,
    copper_gene, broad_iron_gene, core_iron_gene
  )

write_csv(
  candidate_rankings,
  file.path(output_dirs$tables, "metal_candidate_MM_GS_rankings.csv")
)

# Confirm that the four principal panels reproduce the submitted correlations.
expected_main_rho <- c(
  AD_ERC_greenyellow = 0.66,
  FTLD1_red = 0.49,
  FTLD2_darkmagenta = 0.45,
  FTLD3_turquoise = 0.47
)
observed_main_rho <- correlation_results |>
  filter(module_id %in% names(expected_main_rho)) |>
  select(module_id, spearman_rho) |>
  tibble::deframe()

if (!all(round(observed_main_rho[names(expected_main_rho)], 2) == expected_main_rho)) {
  stop("The four principal-panel Spearman coefficients do not match the established results.")
}

# -----------------------------------------------------------------------------
# 7. Figure construction
# -----------------------------------------------------------------------------

main_label_genes <- list(
  AD_ERC_greenyellow = c("DAXX", "TP53", "PARK7", "COMMD1", "COX19", "S100A13"),
  FTLD1_red = c("DAXX", "PARK7", "TP53", "COMMD1", "COX19", "S100A13", "MUC2"),
  FTLD2_darkmagenta = c("BMP6", "SLC39A1", "AOX1", "LMTK2", "BCL2", "FTH1", "NUBP1"),
  FTLD3_turquoise = c("SLC39A8", "BCL2", "NFU1", "NDUFS4", "FTH1", "ACO1", "PARK7", "SOD1")
)

module_title <- function(module_id) {
  str_replace_all(module_id, "_", " ")
}

make_mm_gs_plot <- function(module_id, label_mode = c("specified", "top"), top_n = 12L) {
  label_mode <- match.arg(label_mode)
  plot_data <- module_metrics |>
    filter(.data$module_id == !!module_id)
  ranking_data <- candidate_rankings |>
    filter(.data$module_id == !!module_id)

  if (label_mode == "specified") {
    label_data <- ranking_data |>
      filter(gene_symbol %in% main_label_genes[[.env$module_id]])
  } else {
    label_data <- ranking_data |>
      slice_head(n = top_n)
  }

  rho_row <- correlation_results |>
    filter(.data$module_id == !!module_id)
  colour <- plot_data$plot_colour[[1]]

  ggplot(plot_data, aes(module_membership, gene_significance)) +
    geom_point(colour = "#A8C9EB", alpha = 0.58, size = 1.35) +
    geom_point(
      data = filter(plot_data, highlighted_candidate),
      shape = 8, colour = colour, size = 3.0, stroke = 0.75
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = gene_symbol),
      seed = 2026,
      size = 3.1,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = Inf,
      colour = "black"
    ) +
    annotate(
      "label",
      x = Inf, y = -Inf,
      hjust = 1.08, vjust = -0.45,
      label = sprintf("Spearman rho = %.2f\nP %s", rho_row$spearman_rho, rho_row$p_value_display),
      size = 3.2,
      linewidth = 0.2,
      fill = "white"
    ) +
    labs(
      title = paste0(module_title(module_id), " module"),
      x = "Absolute module membership",
      y = "Absolute gene significance"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 10.5),
      plot.margin = margin(8, 12, 8, 8)
    )
}

main_ids <- module_specification |>
  filter(main_panel) |>
  arrange(panel) |>
  pull(module_id)

figure_6 <- wrap_plots(
  map(main_ids, make_mm_gs_plot, label_mode = "specified"),
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A")

supplementary_ids <- module_specification |>
  filter(!main_panel) |>
  arrange(panel) |>
  pull(module_id)

supplementary_figure_s2 <- wrap_plots(
  map(supplementary_ids, make_mm_gs_plot, label_mode = "top", top_n = 12L),
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A")

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".png")),
    plot, width = width, height = height, units = "in", dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".pdf")),
    plot, width = width, height = height, units = "in", device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "white"
  )
  ggsave(
    file.path(output_dirs$figures, paste0(stem, ".tiff")),
    plot, width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
}

save_plot_set(figure_6, "Figure_6_MM_GS", width = 13, height = 10)
save_plot_set(
  supplementary_figure_s2,
  "Supplementary_Figure_S2_additional_MM_GS",
  width = 13,
  height = 15
)

write_csv(
  filter(module_metrics, main_panel),
  file.path(output_dirs$source_data, "Figure_6_MM_GS_source_data.csv")
)
write_csv(
  filter(module_metrics, !main_panel),
  file.path(output_dirs$source_data, "Supplementary_Figure_S2_MM_GS_source_data.csv")
)

message("R/06_candidate_prioritisation_MM_GS.R completed successfully.")
message("Figures: ", output_dirs$figures)
message("Candidate rankings and correlations: ", output_dirs$tables)
message("Figure source data: ", output_dirs$source_data)
