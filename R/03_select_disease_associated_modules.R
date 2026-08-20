# ==============================================================================
# 03_select_disease_associated_modules.R
# Selection of disease-associated co-methylation modules
#
# Purpose:
#   Record and validate the disease-associated WGCNA modules retained from
#   published AD, FTLD and movement-disorder methylation networks.
#
# Selection criterion:
#   Within each independently constructed network, a module was retained when
#   its disease-status association met:
#
#       module-trait P <= 0.05 / number of modules tested
#
# Important distinction:
#   The modules selected here are the pre-enrichment disease-associated modules.
#   Script 04 subsequently tests these modules for overrepresentation of the
#   copper, broad-iron and core-iron gene sets.
#
# Outputs:
#   - Network-specific Bonferroni thresholds
#   - Complete long-format disease-associated module list
#   - Table 2 source data
# ==============================================================================

# Load shared setup -------------------------------------------------------------

if (!exists("paths", inherits = TRUE)) {
  source(file.path("R", "00_setup.R"))
}

check_packages(
  c(
    "dplyr",
    "purrr",
    "readr",
    "stringr",
    "tibble"
  ),
  analysis_name = "disease-associated module selection"
)

# Analysis settings -------------------------------------------------------------

familywise_alpha <- 0.05

# Network manifest --------------------------------------------------------------

# Each WGCNA network was constructed independently. Therefore, module colours
# identify modules only within their source network.

network_manifest <- tibble::tribble(
  ~network_order,
  ~network_id,
  ~disease_group,
  ~disease_or_comparison,
  ~dataset_or_network,
  ~disease_trait,
  ~modules_tested,
  ~source_study,

  1L,
  "AD_DLPFC",
  "Alzheimer's disease",
  "AD",
  "DLPFC",
  "AD",
  16L,
  "Fodder et al. (2026)",

  2L,
  "AD_ERC",
  "Alzheimer's disease",
  "AD",
  "ERC",
  "AD",
  25L,
  "Fodder et al. (2026)",

  3L,
  "AD_HIPPO",
  "Alzheimer's disease",
  "AD",
  "HIPPO",
  "AD",
  18L,
  "Fodder et al. (2026)",

  4L,
  "FTLD1",
  "Frontotemporal lobar degeneration",
  "FTLD",
  "FTLD1",
  "FTLD",
  33L,
  "Fodder et al. (2023)",

  5L,
  "FTLD2",
  "Frontotemporal lobar degeneration",
  "FTLD",
  "FTLD2",
  "FTLD",
  49L,
  "Fodder et al. (2023)",

  6L,
  "FTLD3",
  "Frontotemporal lobar degeneration",
  "FTLD",
  "FTLD3",
  "FTLD",
  14L,
  "Fodder et al. (2023)",

  7L,
  "MSA",
  "Movement disorders",
  "MSA",
  "Movement-disorder network",
  "MSA",
  32L,
  "Murthy et al. (2024)",

  8L,
  "PD",
  "Movement disorders",
  "PD",
  "Movement-disorder network",
  "PD",
  32L,
  "Murthy et al. (2024)",

  9L,
  "PSP",
  "Movement disorders",
  "PSP",
  "Movement-disorder network",
  "PSP",
  32L,
  "Murthy et al. (2024)"
) |>
  dplyr::mutate(
    bonferroni_threshold =
      familywise_alpha / modules_tested
  )

# Retained disease-associated modules ------------------------------------------

# These module lists reproduce final Table 2. They were obtained from the
# corresponding published module-trait results using the network-specific
# Bonferroni thresholds calculated above.

retained_modules_by_network <- list(
  AD_DLPFC = c(
    "purple",
    "yellow",
    "midnightblue",
    "green",
    "greenyellow",
    "salmon"
  ),

  AD_ERC = c(
    "purple",
    "lightgreen",
    "pink",
    "black",
    "lightyellow",
    "turquoise",
    "midnightblue",
    "tan",
    "greenyellow",
    "lightcyan",
    "darkred",
    "grey60",
    "orange",
    "blue",
    "cyan",
    "darkgrey",
    "royalblue"
  ),

  AD_HIPPO = c(
    "midnightblue",
    "grey60",
    "salmon",
    "magenta",
    "lightgreen",
    "yellow",
    "green",
    "brown"
  ),

  FTLD1 = c(
    "darkturquoise",
    "skyblue",
    "white",
    "brown",
    "violet",
    "grey60",
    "red",
    "blue",
    "darkgreen"
  ),

  FTLD2 = c(
    "darkorange2",
    "plum2",
    "lightcyan1",
    "mediumpurple3",
    "sienna3",
    "darkmagenta",
    "red",
    "ivory",
    "darkgrey",
    "darkorange",
    "tan",
    "midnightblue",
    "orangered4",
    "steelblue",
    "blue",
    "paleturquoise"
  ),

  FTLD3 = c(
    "purple",
    "cyan",
    "salmon",
    "black",
    "turquoise",
    "blue",
    "magenta",
    "pink",
    "yellow",
    "brown"
  ),

  MSA = c(
    "darkturquoise",
    "violet",
    "darkgreen",
    "saddlebrown",
    "lightcyan",
    "white",
    "midnightblue"
  ),

  PD = c(
    "darkturquoise",
    "darkgreen",
    "lightcyan",
    "darkgrey",
    "white",
    "darkred"
  ),

  PSP = c(
    "darkorange",
    "lightgreen",
    "black",
    "skyblue",
    "lightcyan",
    "pink",
    "darkgrey",
    "darkred",
    "steelblue"
  )
)

# Confirm that every network has a retained-module list -------------------------

missing_network_lists <- setdiff(
  network_manifest$network_id,
  names(retained_modules_by_network)
)

unexpected_network_lists <- setdiff(
  names(retained_modules_by_network),
  network_manifest$network_id
)

if (length(missing_network_lists) > 0L) {
  stop(
    paste0(
      "Retained-module lists are missing for: ",
      paste(missing_network_lists, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (length(unexpected_network_lists) > 0L) {
  stop(
    paste0(
      "Unexpected retained-module lists were supplied for: ",
      paste(unexpected_network_lists, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Convert retained modules to long format ---------------------------------------

capitalise_module_label <- function(module) {
  paste0(
    stringr::str_to_upper(
      stringr::str_sub(module, 1, 1)
    ),
    stringr::str_sub(module, 2)
  )
}

disease_associated_modules <- purrr::imap_dfr(
  retained_modules_by_network,
  function(modules, network_name) {
    tibble::tibble(
      network_id = network_name,
      module_order = seq_along(modules),
      module = stringr::str_to_lower(
        stringr::str_trim(modules)
      )
    )
  }
) |>
  dplyr::left_join(
    network_manifest,
    by = "network_id"
  ) |>
  dplyr::mutate(
    module_label = capitalise_module_label(module),
    module_id = paste(
      network_id,
      module,
      sep = "_"
    ),
    selection_method =
      "Network-specific Bonferroni correction",
    selection_threshold =
      bonferroni_threshold
  ) |>
  dplyr::arrange(
    network_order,
    module_order
  ) |>
  dplyr::select(
    network_order,
    module_order,
    module_id,
    disease_group,
    disease_or_comparison,
    dataset_or_network,
    network_id,
    disease_trait,
    module,
    module_label,
    modules_tested,
    bonferroni_threshold,
    selection_method,
    selection_threshold,
    source_study
  )

# Validate counts reported in final Table 2 -------------------------------------

expected_retained_counts <- c(
  AD_DLPFC = 6L,
  AD_ERC = 17L,
  AD_HIPPO = 8L,
  FTLD1 = 9L,
  FTLD2 = 16L,
  FTLD3 = 10L,
  MSA = 7L,
  PD = 6L,
  PSP = 9L
)

observed_retained_counts <- disease_associated_modules |>
  dplyr::count(
    network_id,
    name = "retained_modules"
  ) |>
  tibble::deframe()

incorrect_counts <- names(expected_retained_counts)[
  observed_retained_counts[names(expected_retained_counts)] !=
    expected_retained_counts
]

if (length(incorrect_counts) > 0L) {
  stop(
    paste0(
      "Unexpected retained-module counts for: ",
      paste(incorrect_counts, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (anyDuplicated(disease_associated_modules$module_id) > 0L) {
  stop(
    "Duplicated network-specific module identifiers were detected.",
    call. = FALSE
  )
}

if (any(
  disease_associated_modules$modules_tested <
  expected_retained_counts[
    disease_associated_modules$network_id
  ]
)) {
  stop(
    "A retained-module count exceeds the number of modules tested.",
    call. = FALSE
  )
}

# Optional validation using a complete machine-readable P-value table -----------

# If a complete table of published module-trait P-values is available, place it
# at:
#
# data/published_results/published_module_trait_results.csv
#
# Required columns:
#   network_id
#   module
#   module_trait_p
#
# The file should contain every tested module, not only the retained modules.

published_p_value_file <- file.path(
  paths$published_results,
  "published_module_trait_results.csv"
)

if (file.exists(published_p_value_file)) {
  published_module_traits <- readr::read_csv(
    published_p_value_file,
    show_col_types = FALSE
  )

  required_p_value_columns <- c(
    "network_id",
    "module",
    "module_trait_p"
  )

  missing_p_value_columns <- setdiff(
    required_p_value_columns,
    names(published_module_traits)
  )

  if (length(missing_p_value_columns) > 0L) {
    stop(
      paste0(
        "The published module-trait table is missing: ",
        paste(missing_p_value_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  published_module_traits <- published_module_traits |>
    dplyr::mutate(
      network_id = as.character(network_id),
      module = stringr::str_to_lower(
        stringr::str_trim(module)
      ),
      module_trait_p = as.numeric(module_trait_p)
    ) |>
    dplyr::left_join(
      network_manifest |>
        dplyr::select(
          network_id,
          modules_tested,
          bonferroni_threshold
        ),
      by = "network_id"
    )

  if (anyNA(published_module_traits$module_trait_p)) {
    stop(
      "Missing or non-numeric module-trait P-values were detected.",
      call. = FALSE
    )
  }

  if (anyNA(published_module_traits$bonferroni_threshold)) {
    stop(
      "Unknown network identifiers were detected in the P-value table.",
      call. = FALSE
    )
  }

  tested_counts <- published_module_traits |>
    dplyr::count(
      network_id,
      name = "observed_modules_tested"
    ) |>
    dplyr::left_join(
      network_manifest |>
        dplyr::select(
          network_id,
          expected_modules_tested = modules_tested
        ),
      by = "network_id"
    )

  incorrect_tested_counts <- tested_counts |>
    dplyr::filter(
      observed_modules_tested !=
        expected_modules_tested
    )

  if (nrow(incorrect_tested_counts) > 0L) {
    stop(
      paste0(
        "The complete P-value table does not contain the expected number ",
        "of tested modules for every network."
      ),
      call. = FALSE
    )
  }

  modules_selected_from_p_values <- published_module_traits |>
    dplyr::filter(
      module_trait_p <= bonferroni_threshold
    ) |>
    dplyr::transmute(
      module_id = paste(
        network_id,
        module,
        sep = "_"
      )
    ) |>
    dplyr::pull(module_id)

  recorded_module_ids <-
    disease_associated_modules$module_id

  missing_from_p_value_selection <- setdiff(
    recorded_module_ids,
    modules_selected_from_p_values
  )

  additional_passing_modules <- setdiff(
    modules_selected_from_p_values,
    recorded_module_ids
  )

  if (length(missing_from_p_value_selection) > 0L ||
      length(additional_passing_modules) > 0L) {
    stop(
      paste0(
        "Modules selected from the complete P-value table do not match ",
        "the audited Table 2 module lists."
      ),
      call. = FALSE
    )
  }

  disease_associated_modules <-
    disease_associated_modules |>
    dplyr::left_join(
      published_module_traits |>
        dplyr::select(
          network_id,
          module,
          module_trait_p
        ),
      by = c(
        "network_id",
        "module"
      )
    )

  message(
    "Complete published module-trait P-values were validated successfully."
  )
} else {
  message(
    paste0(
      "Using the audited retained-module lists recorded in final Table 2. ",
      "A complete machine-readable source P-value table was not required ",
      "for downstream analysis."
    )
  )
}

# Produce the Table 2 summary ---------------------------------------------------

table_2_source_data <- disease_associated_modules |>
  dplyr::arrange(
    network_order,
    module_order
  ) |>
  dplyr::group_by(
    network_order,
    disease_or_comparison,
    dataset_or_network,
    modules_tested,
    bonferroni_threshold,
    source_study
  ) |>
  dplyr::summarise(
    retained_modules_n = dplyr::n(),
    Retained_disease_associated_modules = paste(
      module_label,
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(network_order) |>
  dplyr::transmute(
    `Disease / comparison` =
      disease_or_comparison,
    `Dataset or network` =
      dataset_or_network,
    `Modules tested, n` =
      modules_tested,
    `Bonferroni threshold` =
      bonferroni_threshold,
    `Retained modules, n` =
      retained_modules_n,
    `Retained disease-associated modules` =
      Retained_disease_associated_modules,
    Source = source_study
  )

# Save analysis-ready outputs ---------------------------------------------------

write_output_csv(
  network_manifest |>
    dplyr::select(
      network_id,
      disease_group,
      disease_or_comparison,
      dataset_or_network,
      disease_trait,
      modules_tested,
      bonferroni_threshold,
      source_study
    ),
  filename =
    "network_specific_bonferroni_thresholds.csv"
)

write_output_csv(
  disease_associated_modules |>
    dplyr::select(-network_order),
  filename =
    "disease_associated_modules_long.csv"
)

write_output_csv(
  table_2_source_data,
  filename =
    "Table_2_disease_associated_modules.csv",
  output_directory = paths$tables
)

saveRDS(
  disease_associated_modules,
  file.path(
    paths$results,
    "disease_associated_modules.rds"
  )
)

# Selection summary -------------------------------------------------------------

selection_summary <- tibble::tibble(
  measure = c(
    "Source networks",
    "Total modules tested across source networks",
    "Disease-associated modules retained",
    "Family-wise alpha"
  ),
  value = c(
    nrow(network_manifest),
    sum(network_manifest$modules_tested),
    nrow(disease_associated_modules),
    familywise_alpha
  )
)

write_output_csv(
  selection_summary,
  filename =
    "disease_associated_module_selection_summary.csv"
)

capture.output(
  sessionInfo(),
  file = file.path(
    paths$logs,
    "03_module_selection_session_info.txt"
  )
)

message("Disease-associated module selection completed successfully.")
message("Networks represented: ", nrow(network_manifest))
message(
  "Disease-associated modules retained for enrichment testing: ",
  nrow(disease_associated_modules)
)
message(
  "These modules will be tested for metal-gene enrichment in Script 04."
)
