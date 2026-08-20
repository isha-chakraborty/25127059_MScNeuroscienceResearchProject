# ==============================================================================
# 05_integrate_published_EWCE.R
# Integration of published expression-weighted cell-type enrichment results
#
# Purpose:
#   Add cellular context to the nine metal-enriched co-methylation modules by
#   integrating the EWCE findings reported in the original methylation studies.
#
# Important:
#   This script does not rerun EWCE. It integrates published source-study
#   results while retaining their original reference data and significance
#   criteria.
#
# Published EWCE framework:
#   - Zeisel et al. (2015) mouse nervous-system single-cell reference
#   - Mouse-to-human homologue mapping
#   - Source-study empirical P-values and multiple-testing corrections
#
# Output:
#   Table 4 and a complete nine-module EWCE integration table
# ==============================================================================

# Load shared setup -------------------------------------------------------------

if (!exists("paths", inherits = TRUE)) {
  source(file.path("R", "00_setup.R"))
}

check_packages(
  c(
    "dplyr",
    "readr",
    "stringr",
    "tibble"
  ),
  analysis_name = "published EWCE integration"
)

# Load final metal-enrichment results -------------------------------------------

metal_enrichment_file <- file.path(
  paths$results,
  "metal_enriched_modules_and_gene_sets.rds"
)

check_input_files(
  metal_enrichment_file,
  description = "metal-enrichment result file"
)

metal_enrichment_results <- readRDS(
  metal_enrichment_file
)

required_enrichment_columns <- c(
  "module_id",
  "network_id",
  "module",
  "gene_set",
  "disease_group",
  "disease_or_comparison",
  "dataset_or_network",
  "source_study"
)

missing_enrichment_columns <- setdiff(
  required_enrichment_columns,
  names(metal_enrichment_results)
)

if (length(missing_enrichment_columns) > 0L) {
  stop(
    paste0(
      "The metal-enrichment results are missing: ",
      paste(
        missing_enrichment_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# Collapse multiple metal results to one row per module -------------------------

describe_metal_enrichment <- function(gene_sets) {
  gene_sets <- unique(
    as.character(gene_sets)
  )

  if (all(
    c(
      "Broad iron",
      "Core iron"
    ) %in% gene_sets
  )) {
    return("Broad and core iron")
  }

  if ("Copper" %in% gene_sets) {
    return("Copper")
  }

  if ("Core iron" %in% gene_sets) {
    return("Core iron")
  }

  if ("Broad iron" %in% gene_sets) {
    return("Broad iron")
  }

  stop(
    "An unexpected metal-enrichment classification was detected.",
    call. = FALSE
  )
}

module_display_name <- function(
    network_id,
    module) {

  dplyr::case_when(
    network_id == "AD_ERC" ~
      paste("AD ERC", module),

    network_id == "AD_HIPPO" ~
      paste("AD HIPPO", module),

    network_id == "AD_DLPFC" ~
      paste("AD DLPFC", module),

    network_id %in% c(
      "FTLD1",
      "FTLD2",
      "FTLD3",
      "MSA",
      "PD",
      "PSP"
    ) ~
      paste(network_id, module),

    TRUE ~
      paste(network_id, module)
  )
}

metal_enriched_modules <- metal_enrichment_results |>
  dplyr::group_by(
    module_id,
    network_id,
    module,
    disease_group,
    disease_or_comparison,
    dataset_or_network,
    source_study
  ) |>
  dplyr::summarise(
    Metal_enrichment =
      describe_metal_enrichment(gene_set),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Module = module_display_name(
      network_id,
      module
    )
  )

if (nrow(metal_enriched_modules) != 9L) {
  stop(
    paste0(
      "Expected nine metal-enriched modules, but found ",
      nrow(metal_enriched_modules),
      "."
    ),
    call. = FALSE
  )
}

# Curated summary of published EWCE results -------------------------------------

# These entries reproduce the significant cell-type findings reported in the
# source methylation studies for the nine modules retained in the final thesis.
#
# A missing cell-type entry means that the source study did not report a
# statistically significant EWCE result for that module under its corresponding
# multiple-testing procedure.

published_ewce_reference <- tibble::tribble(
  ~module_id,
  ~published_EWCE_significant,
  ~significant_EWCE_cell_types,
  ~EWCE_source_study,
  ~source_significance_criterion,

  "AD_ERC_greenyellow",
  FALSE,
  NA_character_,
  "Fodder et al. (2026)",
  "Significance as reported in the source study",

  "AD_ERC_darkred",
  FALSE,
  NA_character_,
  "Fodder et al. (2026)",
  "Significance as reported in the source study",

  "AD_HIPPO_yellow",
  FALSE,
  NA_character_,
  "Fodder et al. (2026)",
  "Significance as reported in the source study",

  "FTLD1_red",
  TRUE,
  "CA1 pyramidal neurons",
  "Fodder et al. (2023)",
  "Bonferroni-adjusted P < 0.05",

  "FTLD2_darkmagenta",
  FALSE,
  NA_character_,
  "Fodder et al. (2023)",
  "Bonferroni-adjusted P < 0.05",

  "FTLD3_black",
  TRUE,
  "Somatosensory and CA1 pyramidal neurons",
  "Fodder et al. (2023)",
  "Bonferroni-adjusted P < 0.05",

  "FTLD3_turquoise",
  TRUE,
  paste(
    "Somatosensory and CA1 pyramidal neurons;",
    "oligodendrocytes; interneurons"
  ),
  "Fodder et al. (2023)",
  "Bonferroni-adjusted P < 0.05",

  "MSA_violet",
  FALSE,
  NA_character_,
  "Murthy et al. (2024)",
  "Bonferroni-adjusted P < 0.05",

  "PSP_pink",
  FALSE,
  NA_character_,
  "Murthy et al. (2024)",
  "Bonferroni-adjusted P < 0.05"
) |>
  dplyr::mutate(
    EWCE_reference_dataset =
      paste(
        "Zeisel et al. (2015) mouse nervous-system",
        "single-cell reference with human-homologue mapping"
      ),
    integration_type =
      "Published EWCE result"
  )

# Validate module correspondence ------------------------------------------------

expected_module_ids <- sort(
  published_ewce_reference$module_id
)

observed_module_ids <- sort(
  metal_enriched_modules$module_id
)

if (!identical(
  expected_module_ids,
  observed_module_ids
)) {
  missing_ewce_records <- setdiff(
    observed_module_ids,
    expected_module_ids
  )

  unexpected_ewce_records <- setdiff(
    expected_module_ids,
    observed_module_ids
  )

  stop(
    paste0(
      "The EWCE reference table does not match the nine ",
      "metal-enriched modules.\n",
      "Missing EWCE records: ",
      paste(missing_ewce_records, collapse = ", "),
      "\nUnexpected EWCE records: ",
      paste(unexpected_ewce_records, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Integrate EWCE and metal-enrichment results -----------------------------------

ewce_integration <- metal_enriched_modules |>
  dplyr::left_join(
    published_ewce_reference,
    by = "module_id"
  ) |>
  dplyr::mutate(
    Published_EWCE_result = dplyr::if_else(
      published_EWCE_significant,
      significant_EWCE_cell_types,
      "No significant enrichment reported"
    )
  )

if (anyNA(
  ewce_integration$published_EWCE_significant
)) {
  stop(
    "At least one metal-enriched module lacks an EWCE classification.",
    call. = FALSE
  )
}

if (sum(
  ewce_integration$published_EWCE_significant
) != 3L) {
  stop(
    "Expected significant published EWCE results for three modules.",
    call. = FALSE
  )
}

expected_significant_modules <- c(
  "FTLD1_red",
  "FTLD3_black",
  "FTLD3_turquoise"
)

observed_significant_modules <-
  ewce_integration |>
  dplyr::filter(
    published_EWCE_significant
  ) |>
  dplyr::pull(module_id)

if (!setequal(
  expected_significant_modules,
  observed_significant_modules
)) {
  stop(
    "The significant EWCE modules differ from the submitted thesis.",
    call. = FALSE
  )
}

# Set final module order --------------------------------------------------------

final_module_order <- c(
  "AD_ERC_greenyellow",
  "AD_ERC_darkred",
  "AD_HIPPO_yellow",
  "FTLD1_red",
  "FTLD2_darkmagenta",
  "FTLD3_black",
  "FTLD3_turquoise",
  "MSA_violet",
  "PSP_pink"
)

ewce_integration <- ewce_integration |>
  dplyr::mutate(
    module_order = match(
      module_id,
      final_module_order
    )
  ) |>
  dplyr::arrange(module_order)

# Prepare final Table 4 ---------------------------------------------------------

table_4_module_order <- c(
  "FTLD1_red",
  "FTLD3_black",
  "FTLD3_turquoise"
)

table_4 <- ewce_integration |>
  dplyr::filter(
    published_EWCE_significant
  ) |>
  dplyr::mutate(
    table_order = match(
      module_id,
      table_4_module_order
    )
  ) |>
  dplyr::arrange(table_order) |>
  dplyr::transmute(
    Module,
    `Metal enrichment` =
      Metal_enrichment,
    `Significant published EWCE cell type(s)` =
      significant_EWCE_cell_types
  )

if (nrow(table_4) != 3L) {
  stop(
    "Table 4 should contain exactly three modules.",
    call. = FALSE
  )
}

# Save source-study EWCE manifest -----------------------------------------------

readr::write_csv(
  published_ewce_reference,
  file.path(
    paths$published_results,
    "published_EWCE_summary.csv"
  ),
  na = ""
)

# Save complete integration results ---------------------------------------------

write_output_csv(
  ewce_integration |>
    dplyr::select(
      module_id,
      Module,
      disease_group,
      dataset_or_network,
      Metal_enrichment,
      published_EWCE_significant,
      significant_EWCE_cell_types,
      Published_EWCE_result,
      EWCE_source_study,
      source_significance_criterion,
      EWCE_reference_dataset,
      integration_type
    ),
  filename =
    "published_EWCE_integration_all_9_modules.csv"
)

# Save Table 4 ------------------------------------------------------------------

write_output_csv(
  table_4,
  filename =
    "Table_4_published_EWCE_results.csv",
  output_directory = paths$tables
)

saveRDS(
  ewce_integration,
  file.path(
    paths$results,
    "published_EWCE_integration.rds"
  )
)

# Optional Word version of Table 4 ----------------------------------------------

if (
  requireNamespace(
    "flextable",
    quietly = TRUE
  ) &&
  requireNamespace(
    "officer",
    quietly = TRUE
  )
) {
  table_4_flextable <-
    flextable::flextable(table_4) |>
    flextable::theme_booktabs() |>
    flextable::bg(
      part = "header",
      bg = "#374151"
    ) |>
    flextable::color(
      part = "header",
      color = "white"
    ) |>
    flextable::bold(
      part = "header"
    ) |>
    flextable::font(
      fontname = "Times New Roman",
      part = "all"
    ) |>
    flextable::fontsize(
      size = 9.5,
      part = "body"
    ) |>
    flextable::fontsize(
      size = 10,
      part = "header"
    ) |>
    flextable::align(
      align = "center",
      part = "header"
    ) |>
    flextable::align(
      j = c(
        "Module",
        "Metal enrichment"
      ),
      align = "center",
      part = "body"
    ) |>
    flextable::align(
      j =
        "Significant published EWCE cell type(s)",
      align = "left",
      part = "body"
    ) |>
    flextable::valign(
      valign = "center",
      part = "all"
    ) |>
    flextable::padding(
      padding.top = 6,
      padding.bottom = 6,
      padding.left = 6,
      padding.right = 6,
      part = "all"
    ) |>
    flextable::width(
      j = "Module",
      width = 1.35
    ) |>
    flextable::width(
      j = "Metal enrichment",
      width = 1.50
    ) |>
    flextable::width(
      j =
        "Significant published EWCE cell type(s)",
      width = 3.80
    ) |>
    flextable::set_table_properties(
      layout = "fixed",
      width = 1
    )

  # Copper-associated row
  table_4_flextable <-
    flextable::bg(
      table_4_flextable,
      i = 1,
      part = "body",
      bg = "#FFF4E6"
    ) |>
    flextable::color(
      i = 1,
      j = "Metal enrichment",
      part = "body",
      color = "#B45309"
    )

  # Iron-associated rows
  table_4_flextable <-
    flextable::bg(
      table_4_flextable,
      i = 2:3,
      part = "body",
      bg = "#EAF2F8"
    ) |>
    flextable::color(
      i = 2:3,
      j = "Metal enrichment",
      part = "body",
      color = "#1F4E78"
    ) |>
    flextable::bold(
      j = c(
        "Module",
        "Metal enrichment"
      ),
      part = "body"
    )

  table_4_title <- paste0(
    "Table 4. Neuronal and glial expression signatures ",
    "identified in three metal-enriched FTLD ",
    "co-methylation modules."
  )

  table_4_legend <- paste0(
    "Significant expression-weighted cell-type enrichment results ",
    "reported by Fodder et al. (2023) for three of the nine ",
    "metal-enriched disease-associated modules identified in the ",
    "present study. EWCE assesses whether genes assigned to a module ",
    "are expressed more highly in a reference cell type than expected ",
    "from matched random gene sets. Metal enrichment refers to the ",
    "copper, broad-iron or core-iron gene set overrepresented within ",
    "each module. Only modules with significant published EWCE ",
    "results are shown."
  )

  table_4_note <- paste0(
    "Abbreviations: EWCE, expression-weighted cell-type enrichment; ",
    "FTLD, frontotemporal lobar degeneration."
  )

  table_4_document <- officer::read_docx()

  table_4_document <- officer::body_add_par(
    table_4_document,
    value = table_4_title,
    style = "heading 2"
  )

  table_4_document <- officer::body_add_par(
    table_4_document,
    value = table_4_legend,
    style = "Normal"
  )

  table_4_document <-
    flextable::body_add_flextable(
      table_4_document,
      value = table_4_flextable
    )

  table_4_document <- officer::body_add_par(
    table_4_document,
    value = table_4_note,
    style = "Normal"
  )

  table_4_word_file <- file.path(
    paths$tables,
    "Table_4_published_EWCE_results.docx"
  )

  print(
    table_4_document,
    target = table_4_word_file
  )

  message(
    "Saved Word table: ",
    table_4_word_file
  )
} else {
  message(
    paste0(
      "Table 4 CSV was saved. Install the optional officer and ",
      "flextable packages to generate the Word version."
    )
  )
}

# Integration summary -----------------------------------------------------------

ewce_summary <- tibble::tibble(
  measure = c(
    "Metal-enriched modules evaluated",
    "Modules with significant published EWCE",
    "Copper-enriched modules with significant EWCE",
    "Iron-enriched modules with significant EWCE"
  ),
  value = c(
    nrow(ewce_integration),
    sum(
      ewce_integration$
        published_EWCE_significant
    ),
    ewce_integration |>
      dplyr::filter(
        published_EWCE_significant,
        Metal_enrichment == "Copper"
      ) |>
      nrow(),
    ewce_integration |>
      dplyr::filter(
        published_EWCE_significant,
        stringr::str_detect(
          Metal_enrichment,
          "iron"
        )
      ) |>
      nrow()
  )
)

write_output_csv(
  ewce_summary,
  filename =
    "published_EWCE_integration_summary.csv"
)

capture.output(
  sessionInfo(),
  file = file.path(
    paths$logs,
    "05_published_EWCE_integration_session_info.txt"
  )
)

message("Published EWCE integration completed successfully.")
message("Metal-enriched modules evaluated: 9")
message("Modules with significant published EWCE: 3")
message(
  "Significant modules: FTLD1 red, FTLD3 black and FTLD3 turquoise"
)
