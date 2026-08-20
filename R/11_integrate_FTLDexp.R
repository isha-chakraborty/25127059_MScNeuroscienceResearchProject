
# Integration of corrected FTLD methylation-module candidates with FTLDexp
# differential-expression results
#
# The four retained FTLD modules are:
#   FTLD1 red, FTLD2 darkmagenta, FTLD3 black and FTLD3 turquoise.
#
# Every curated metal-associated gene assigned to each module is evaluated in
# FTLDexp. Copper, broad-iron and core-iron membership are recorded separately so
# genes contributing directly to the module's enrichment can be distinguished
# from other metal-associated genes in the same module.
#
# Generates Figure 12, Table 6 and Supplementary Table S7.

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "readr", "stringr", "tibble", "tidyr"
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
  library(patchwork)
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
  tables = file.path(project_root, "outputs", "tables"),
  supplementary = file.path(project_root, "outputs", "supplementary_tables"),
  derived_data = file.path(project_root, "outputs", "derived_data"),
  source_data = file.path(project_root, "outputs", "figure_source_data")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# -----------------------------------------------------------------------------
# 2. Four corrected FTLD modules
# -----------------------------------------------------------------------------

module_key <- tribble(
  ~module_order, ~module_id,            ~module_label,          ~enriched_metal_set,     ~module_colour,
  1L,            "FTLD1_red",          "FTLD1 red",            "Copper",                "#D73027",
  2L,            "FTLD2_darkmagenta",  "FTLD2 darkmagenta",    "Broad and core iron",   "#9E006E",
  3L,            "FTLD3_black",        "FTLD3 black",          "Broad and core iron",   "#2B2B2B",
  4L,            "FTLD3_turquoise",    "FTLD3 turquoise",      "Broad and core iron",   "#109C95"
)

candidate_occurrence_file <- file.path(
  project_root,
  "outputs", "derived_data", "nine_module_candidate_occurrence_long.csv"
)
if (!file.exists(candidate_occurrence_file)) {
  stop("Run R/07_recurrence_analysis_9_modules.R before FTLDexp integration.")
}

candidate_occurrence <- read_csv(
  candidate_occurrence_file,
  show_col_types = FALSE
)

required_candidate_columns <- c(
  "module_id", "gene_symbol", "copper", "broad_iron", "core_iron"
)
if (!all(required_candidate_columns %in% names(candidate_occurrence))) {
  stop("The candidate-occurrence file does not contain the expected Script 07 columns.")
}

# Candidate records are the union of curated copper and broad/core-iron genes in
# each module. Core iron is a subset of broad iron.
ftld_candidates <- candidate_occurrence |>
  filter(module_id %in% module_key$module_id) |>
  transmute(
    module_id,
    gene_symbol = str_to_upper(str_trim(as.character(gene_symbol))),
    copper = as.logical(copper),
    broad_iron = as.logical(broad_iron),
    core_iron = as.logical(core_iron)
  ) |>
  filter(copper | broad_iron | core_iron) |>
  distinct() |>
  left_join(module_key, by = "module_id") |>
  mutate(
    contributed_to_enrichment = case_when(
      module_id == "FTLD1_red" & copper ~ TRUE,
      module_id != "FTLD1_red" & broad_iron ~ TRUE,
      TRUE ~ FALSE
    ),
    candidate_basis = "Curated metal-associated gene assigned to the retained module"
  ) |>
  arrange(module_order, gene_symbol)

expected_candidate_occurrences <- c(
  FTLD1_red = 20L,
  FTLD2_darkmagenta = 21L,
  FTLD3_black = 72L,
  FTLD3_turquoise = 75L
)
observed_candidate_occurrences <- ftld_candidates |>
  count(module_id) |>
  tibble::deframe()

if (!identical(
  as.integer(observed_candidate_occurrences[names(expected_candidate_occurrences)]),
  as.integer(expected_candidate_occurrences)
)) {
  print(ftld_candidates |> count(module_id))
  stop("The four-module candidate counts differ from the corrected analysis.")
}
if (nrow(ftld_candidates) != 188L || n_distinct(ftld_candidates$gene_symbol) != 117L) {
  stop("Expected 188 candidate occurrences representing 117 unique genes.")
}

# -----------------------------------------------------------------------------
# 3. Load the published FTLDexp differential-expression results
# -----------------------------------------------------------------------------

# A neutral CSV is preferred for the public repository. Expected columns are
# gene_symbol, logFC, P.Value and adj.P.Val. The completed local analysis can also
# be reconstructed from the saved differential-expression and gene-mapping RData
# objects by inspecting their contents rather than relying on personal paths.
expression_csv_candidates <- c(
  file.path(
    project_root, "data", "published_results",
    "FTLDexp_differential_expression.csv"
  ),
  file.path(
    project_root, "data", "FTLDexp",
    "FTLDexp_differential_expression.csv"
  )
)
expression_csv <- expression_csv_candidates[file.exists(expression_csv_candidates)][1]

standardise_expression_table <- function(data) {
  names(data) <- str_trim(names(data))

  if (all(c(
    "gene_symbol", "log2_fold_change", "nominal_p", "adjusted_p"
  ) %in% names(data))) {
    return(
      data |>
        transmute(
          gene_symbol = str_to_upper(str_trim(as.character(gene_symbol))),
          log2_fold_change = suppressWarnings(as.numeric(log2_fold_change)),
          nominal_p = suppressWarnings(as.numeric(nominal_p)),
          adjusted_p = suppressWarnings(as.numeric(adjusted_p))
        ) |>
        filter(!is.na(gene_symbol), gene_symbol != "") |>
        arrange(adjusted_p, nominal_p) |>
        distinct(gene_symbol, .keep_all = TRUE)
    )
  }

  symbol_candidates <- c(
    "gene_symbol", "geneSymbol", "gene_symbols", "symbol", "SYMBOL"
  )
  symbol_column <- symbol_candidates[symbol_candidates %in% names(data)][1]
  if (is.na(symbol_column)) {
    stop("The expression table does not contain a gene-symbol column.")
  }

  required_statistics <- c("logFC", "P.Value", "adj.P.Val")
  if (!all(required_statistics %in% names(data))) {
    stop("The expression table lacks logFC, P.Value or adj.P.Val.")
  }

  data |>
    transmute(
      gene_symbol = str_to_upper(str_trim(as.character(.data[[symbol_column]]))),
      log2_fold_change = suppressWarnings(as.numeric(logFC)),
      nominal_p = suppressWarnings(as.numeric(P.Value)),
      adjusted_p = suppressWarnings(as.numeric(adj.P.Val))
    ) |>
    filter(!is.na(gene_symbol), gene_symbol != "") |>
    arrange(adjusted_p, nominal_p) |>
    distinct(gene_symbol, .keep_all = TRUE)
}

if (!is.na(expression_csv)) {
  expression_results <- read_csv(expression_csv, show_col_types = FALSE) |>
    standardise_expression_table()
} else {
  search_directories <- c(
    file.path(project_root, "FTLDexp"),
    file.path(project_root, "data", "FTLDexp"),
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

  differential_table <- NULL
  gene_mapping_candidates <- list()

  for (rdata_file in rdata_files) {
    temporary_environment <- new.env(parent = emptyenv())
    loaded_names <- tryCatch(
      load(rdata_file, envir = temporary_environment),
      error = function(error) character()
    )

    for (object_name in loaded_names) {
      object <- get(object_name, envir = temporary_environment)
      if (!is.data.frame(object)) next

      if (all(c("logFC", "P.Value", "adj.P.Val") %in% names(object))) {
        differential_table <- object
      }
      if (all(c("gene", "gene_symbols") %in% names(object))) {
        gene_mapping_candidates[[paste(rdata_file, object_name, sep = "::")]] <- object |>
          transmute(
            gene_id = as.character(gene),
            gene_symbol = as.character(gene_symbols)
          ) |>
          distinct()
      }
    }
  }

  if (is.null(differential_table) || length(gene_mapping_candidates) == 0L) {
    stop(
      "FTLDexp differential-expression or gene-mapping objects were not found. ",
      "Provide data/published_results/FTLDexp_differential_expression.csv."
    )
  }

  if (!"gene" %in% names(differential_table)) {
    differential_table$gene <- rownames(differential_table)
  }

  # Select the mapping with the greatest overlap with the differential-expression
  # table. This identifies the disease-matched FTLDexp mapping when other saved
  # expression objects are present in the same directory.
  differential_gene_ids <- unique(as.character(differential_table$gene))
  mapping_overlap <- vapply(
    gene_mapping_candidates,
    function(mapping) sum(differential_gene_ids %in% mapping$gene_id),
    integer(1)
  )
  gene_mapping <- gene_mapping_candidates[[which.max(mapping_overlap)]]

  expression_results <- differential_table |>
    transmute(
      gene_id = as.character(gene),
      logFC,
      P.Value,
      adj.P.Val
    ) |>
    left_join(gene_mapping, by = "gene_id") |>
    standardise_expression_table()
}

write_csv(
  expression_results,
  file.path(output_dirs$derived_data, "FTLDexp_differential_expression_clean.csv")
)

# -----------------------------------------------------------------------------
# 4. Join module-derived candidates to FTLDexp by gene symbol
# -----------------------------------------------------------------------------

integration <- ftld_candidates |>
  left_join(expression_results, by = "gene_symbol") |>
  mutate(
    matched_in_FTLDexp = !is.na(adjusted_p),
    significant = matched_in_FTLDexp & adjusted_p < 0.05,
    direction = case_when(
      significant & log2_fold_change > 0 ~ "Upregulated",
      significant & log2_fold_change < 0 ~ "Downregulated",
      significant ~ "No directional change",
      TRUE ~ NA_character_
    ),
    minus_log10_adjusted_p = if_else(
      matched_in_FTLDexp,
      -log10(adjusted_p),
      NA_real_
    )
  ) |>
  arrange(module_order, adjusted_p, gene_symbol)

significant_occurrences <- integration |>
  filter(significant)

expected_significant_by_module <- c(
  FTLD1_red = 3L,
  FTLD2_darkmagenta = 4L,
  FTLD3_black = 18L,
  FTLD3_turquoise = 15L
)
observed_significant_by_module <- significant_occurrences |>
  count(module_id) |>
  tibble::deframe()

if (!identical(
  as.integer(observed_significant_by_module[names(expected_significant_by_module)]),
  as.integer(expected_significant_by_module)
)) {
  print(significant_occurrences |> count(module_id))
  stop("Module-level significant FTLDexp counts do not match the established results.")
}
if (nrow(significant_occurrences) != 40L ||
    n_distinct(significant_occurrences$gene_symbol) != 25L) {
  stop("Expected 40 significant module–gene occurrences representing 25 unique genes.")
}

# The FTLD1 red module contained seven copper genes responsible for its module
# enrichment. Six mapped to FTLDexp, none was significantly differentially
# expressed, and MUC2 was unmatched.
ftld1_copper_audit <- integration |>
  filter(module_id == "FTLD1_red", copper)
if (nrow(ftld1_copper_audit) != 7L ||
    sum(ftld1_copper_audit$matched_in_FTLDexp) != 6L ||
    sum(ftld1_copper_audit$significant) != 0L ||
    !identical(
      sort(ftld1_copper_audit$gene_symbol[!ftld1_copper_audit$matched_in_FTLDexp]),
      "MUC2"
    )) {
  stop("The FTLD1 red copper-contributor audit does not match the final analysis.")
}

# -----------------------------------------------------------------------------
# 5. Table 6 and Supplementary Table S7
# -----------------------------------------------------------------------------

collapse_genes <- function(
    gene_symbol, direction, log2_fold_change, target_direction) {
  keep <- direction == target_direction & !is.na(direction)
  selected <- gene_symbol[keep]
  selected_fold_change <- log2_fold_change[keep]
  if (length(selected) == 0L) return("—")
  selected <- selected[order(
    selected_fold_change,
    decreasing = identical(target_direction, "Upregulated")
  )]
  paste(selected, collapse = ", ")
}

table_6 <- integration |>
  group_by(
    module_order, module_id, module_label,
    enriched_metal_set
  ) |>
  summarise(
    candidate_occurrences_evaluated = n(),
    matched_in_FTLDexp = sum(matched_in_FTLDexp),
    significant_genes_n = sum(significant),
    upregulated_n = sum(direction == "Upregulated", na.rm = TRUE),
    downregulated_n = sum(direction == "Downregulated", na.rm = TRUE),
    upregulated_genes = collapse_genes(
      gene_symbol[significant], direction[significant],
      log2_fold_change[significant], "Upregulated"
    ),
    downregulated_genes = collapse_genes(
      gene_symbol[significant], direction[significant],
      log2_fold_change[significant], "Downregulated"
    ),
    .groups = "drop"
  ) |>
  arrange(module_order) |>
  select(-module_order, -module_id)

write_csv(table_6, file.path(output_dirs$tables, "Table_6_FTLDexp_summary.csv"))

# S7 retains every evaluated module–gene occurrence, including unmatched and
# non-significant candidates, so the matching and filtering steps are auditable.
supplementary_table_s7 <- integration |>
  transmute(
    Module = module_label,
    `Gene symbol` = gene_symbol,
    `Copper-set member` = copper,
    `Broad-iron-set member` = broad_iron,
    `Core-iron-set member` = core_iron,
    `Contributed to module enrichment` = contributed_to_enrichment,
    `Matched in FTLDexp` = matched_in_FTLDexp,
    `Log2 fold change` = log2_fold_change,
    `Nominal P` = nominal_p,
    `FDR-adjusted P` = adjusted_p,
    `Significantly differentially expressed` = significant,
    Direction = direction
  )

write_csv(
  supplementary_table_s7,
  file.path(
    output_dirs$supplementary,
    "Supplementary_Table_S7_complete_FTLDexp_candidate_integration.csv"
  )
)
# The thesis display contained the 40 significant module–gene occurrences. Keep
# that submitted subset alongside the complete 188-row audit table so the two
# scopes are explicit rather than conflated.
write_csv(
  supplementary_table_s7 |>
    filter(`Significantly differentially expressed`) |>
    arrange(Module, `FDR-adjusted P`, `Gene symbol`),
  file.path(
    output_dirs$supplementary,
    "Supplementary_Table_S7_submitted_40_significant_occurrences.csv"
  )
)
write_csv(
  significant_occurrences,
  file.path(
    output_dirs$derived_data,
    "FTLDexp_40_significant_module_gene_occurrences.csv"
  )
)
write_csv(
  significant_occurrences |>
    arrange(adjusted_p, gene_symbol) |>
    distinct(gene_symbol, .keep_all = TRUE) |>
    select(
      gene_symbol, log2_fold_change, nominal_p,
      adjusted_p, direction
    ),
  file.path(
    output_dirs$derived_data,
    "FTLDexp_25_unique_significant_genes.csv"
  )
)

# Optional editable Word table.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(table_6) |>
    flextable::theme_booktabs() |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::autofit()

  document <- officer::read_docx() |>
    officer::body_add_par(
      "Table 6. Differential expression of metal-associated candidates from four retained FTLD co-methylation modules.",
      style = "Table Caption"
    ) |>
    flextable::body_add_flextable(ft)

  print(document, target = file.path(output_dirs$tables, "Table_6_FTLDexp_summary.docx"))
}

# -----------------------------------------------------------------------------
# 6. Figure 12: four-module differential-expression lollipop plot
# -----------------------------------------------------------------------------

lighten_colour <- function(colour, amount = 0.48) {
  rgb_values <- col2rgb(colour) / 255
  light_values <- rgb_values + (1 - rgb_values) * amount
  rgb(light_values[1, ], light_values[2, ], light_values[3, ])
}

module_colours <- setNames(module_key$module_colour, module_key$module_id)
direction_colours <- unlist(lapply(names(module_colours), function(module_id) {
  c(
    setNames(module_colours[[module_id]], paste0(module_id, "__Upregulated")),
    setNames(
      lighten_colour(module_colours[[module_id]]),
      paste0(module_id, "__Downregulated")
    )
  )
}))

figure_data <- significant_occurrences |>
  mutate(
    colour_key = paste(module_id, direction, sep = "__"),
    gene_panel_id = paste(module_id, gene_symbol, sep = "___")
  ) |>
  group_by(module_order) |>
  arrange(log2_fold_change, .by_group = TRUE) |>
  ungroup()

gene_levels <- figure_data |>
  arrange(module_order, log2_fold_change, gene_symbol) |>
  pull(gene_panel_id)
figure_data$gene_panel_id <- factor(figure_data$gene_panel_id, levels = gene_levels)

make_module_panel <- function(module_id, show_x_title = FALSE) {
  panel_data <- figure_data |>
    filter(.data$module_id == !!module_id)
  module_information <- module_key |>
    filter(.data$module_id == !!module_id)

  ggplot(panel_data, aes(log2_fold_change, gene_panel_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#8A8A8A", linewidth = 0.5) +
    geom_segment(
      aes(
        x = 0, xend = log2_fold_change,
        y = gene_panel_id, yend = gene_panel_id,
        colour = colour_key
      ),
      linewidth = 0.9
    ) +
    geom_point(
      aes(size = minus_log10_adjusted_p, colour = colour_key),
      alpha = 0.98
    ) +
    scale_colour_manual(values = direction_colours, guide = "none") +
    scale_size_continuous(
      name = expression(-log[10]("FDR-adjusted P")),
      range = c(2.5, 6.5)
    ) +
    scale_y_discrete(
      labels = function(labels) str_replace(labels, "^.*___", "")
    ) +
    scale_x_continuous(
      limits = c(-2.1, 2.1),
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    labs(
      title = module_information$module_label,
      x = if (show_x_title) expression("FTLDexp RNA-seq " * log[2] * " fold change") else NULL,
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        colour = module_information$module_colour,
        size = 12
      ),
      axis.text.y = element_text(face = "italic", size = 8.5),
      axis.title.x = element_text(size = 10.5),
      panel.border = element_rect(colour = "#C6CFD8", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 8)
    )
}

figure_12 <- wrap_plots(
  make_module_panel("FTLD1_red", show_x_title = FALSE),
  make_module_panel("FTLD2_darkmagenta", show_x_title = FALSE),
  make_module_panel("FTLD3_black", show_x_title = TRUE),
  make_module_panel("FTLD3_turquoise", show_x_title = TRUE),
  ncol = 2,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dirs$figures, "Figure_12_FTLDexp_four_modules.png"),
  figure_12, width = 12, height = 12, dpi = 600, bg = "white"
)
ggsave(
  file.path(output_dirs$figures, "Figure_12_FTLDexp_four_modules.tiff"),
  figure_12, width = 12, height = 12, dpi = 600,
  compression = "lzw", bg = "white"
)
ggsave(
  file.path(output_dirs$figures, "Figure_12_FTLDexp_four_modules.pdf"),
  figure_12, width = 12, height = 12,
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

write_csv(
  figure_data |>
    transmute(
      module = module_label,
      gene_symbol,
      copper,
      broad_iron,
      core_iron,
      contributed_to_enrichment,
      log2_fold_change,
      nominal_p,
      adjusted_p,
      direction,
      minus_log10_adjusted_p
    ),
  file.path(output_dirs$source_data, "Figure_12_source_data.csv")
)

message("R/11_integrate_FTLDexp.R completed successfully.")
message("Evaluated 188 module–gene occurrences representing 117 unique candidates.")
message("Verified 40 significant occurrences representing 25 unique genes.")
message("Figure 12, Table 6 and Supplementary Table S7 were generated.")
