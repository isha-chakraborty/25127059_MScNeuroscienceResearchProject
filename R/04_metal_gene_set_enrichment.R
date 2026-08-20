# ==============================================================================
# 04_metal_gene_set_enrichment.R
# Copper- and iron-gene-set enrichment in co-methylation modules
#
# Purpose:
#   Test whether copper, broad-iron or core-iron genes are overrepresented in
#   each WGCNA module relative to the corresponding network gene background.
#
# Statistical analysis:
#   - Unique annotated gene symbols form each network background
#   - Unique genes are used within each module
#   - Curated genes absent from the network are excluded
#   - Two-sided Fisher's exact test
#   - Benjamini-Hochberg FDR correction across all modules within each
#     network and metal-gene-set analysis
#   - Positive enrichment: odds ratio > 1 and nominal P < 0.05
#
# Final outputs:
#   - Nine metal-enriched disease-associated modules
#   - Twelve module-gene-set enrichment results
#   - Tables 3A and 3B
#   - Figure 5
#   - Supplementary Table S4
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
    "tibble",
    "ggplot2",
    "patchwork",
    "scales"
  ),
  analysis_name = "metal gene-set enrichment"
)

# Required inputs ---------------------------------------------------------------

gene_set_file <- file.path(
  paths$gene_sets,
  "metal_gene_sets.rds"
)

disease_module_file <- file.path(
  paths$results,
  "disease_associated_modules_long.csv"
)

check_input_files(
  c(
    gene_set_file,
    disease_module_file
  ),
  description = "gene-set enrichment inputs"
)

metal_gene_sets <- readRDS(gene_set_file)

required_gene_sets <- c(
  "copper",
  "broad_iron",
  "core_iron"
)

if (!all(required_gene_sets %in% names(metal_gene_sets))) {
  stop(
    paste0(
      "The metal-gene-set object must contain: ",
      paste(required_gene_sets, collapse = ", ")
    ),
    call. = FALSE
  )
}

metal_gene_sets <- lapply(
  metal_gene_sets[required_gene_sets],
  function(genes) {
    genes |>
      standardise_gene_symbols() |>
      unique() |>
      sort()
  }
)

disease_associated_modules <- readr::read_csv(
  disease_module_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(
    network_id = as.character(network_id),
    module = stringr::str_to_lower(
      stringr::str_trim(module)
    )
  )

required_module_columns <- c(
  "network_id",
  "module",
  "module_id",
  "disease_group",
  "disease_or_comparison",
  "dataset_or_network"
)

missing_module_columns <- setdiff(
  required_module_columns,
  names(disease_associated_modules)
)

if (length(missing_module_columns) > 0L) {
  stop(
    paste0(
      "The disease-associated module file is missing: ",
      paste(missing_module_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Network annotation inputs -----------------------------------------------------

# Preferred public input files contain only:
#   geneSymbol
#   moduleColor
#
# This avoids distributing large working objects containing unrelated columns.

network_annotation_directory <- file.path(
  paths$published_results,
  "network_gene_annotations"
)

dir.create(
  network_annotation_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

network_input_files <- c(
  AD_DLPFC = file.path(
    network_annotation_directory,
    "AD_DLPFC_gene_modules.csv"
  ),
  AD_ERC = file.path(
    network_annotation_directory,
    "AD_ERC_gene_modules.csv"
  ),
  AD_HIPPO = file.path(
    network_annotation_directory,
    "AD_HIPPO_gene_modules.csv"
  ),
  FTLD1 = file.path(
    network_annotation_directory,
    "FTLD1_gene_modules.csv"
  ),
  FTLD2 = file.path(
    network_annotation_directory,
    "FTLD2_gene_modules.csv"
  ),
  FTLD3 = file.path(
    network_annotation_directory,
    "FTLD3_gene_modules.csv"
  ),
  MSA = file.path(
    network_annotation_directory,
    "MSA_gene_modules.csv"
  ),
  PD = file.path(
    network_annotation_directory,
    "PD_gene_modules.csv"
  ),
  PSP = file.path(
    network_annotation_directory,
    "PSP_gene_modules.csv"
  )
)

# Helper for cleaning network annotations --------------------------------------

clean_network_annotation <- function(gene_information) {
  required_columns <- c(
    "geneSymbol",
    "moduleColor"
  )

  missing_columns <- setdiff(
    required_columns,
    names(gene_information)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "A network annotation object is missing: ",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gene_information |>
    dplyr::transmute(
      geneSymbol = standardise_gene_symbols(geneSymbol),
      module = stringr::str_to_lower(
        stringr::str_trim(
          as.character(moduleColor)
        )
      )
    ) |>
    dplyr::filter(
      !is.na(geneSymbol),
      geneSymbol != "",
      !is.na(module),
      module != ""
    ) |>
    dplyr::distinct(
      geneSymbol,
      module
    ) |>
    dplyr::arrange(
      module,
      geneSymbol
    )
}

# Convert existing working objects into neutral public inputs ------------------

# This block supports the current local project without including personal
# filenames or computer paths in the public script. If the neutral CSV files
# already exist, this conversion is skipped.

if (!all(file.exists(network_input_files))) {
  message(
    "Neutral network annotation files are absent; checking existing working objects."
  )

  load_to_environment <- function(filename) {
    object_environment <- new.env(
      parent = emptyenv()
    )

    load(
      filename,
      envir = object_environment
    )

    object_environment
  }

  # Locate an AD object file by its contents rather than a personal filename.
  ad_object_candidates <- list.files(
    project_root,
    pattern = "^GeneInfos_.*\\.RData$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  ad_object_candidates <- ad_object_candidates[
    !stringr::str_detect(
      basename(ad_object_candidates),
      stringr::regex(
        "copy",
        ignore_case = TRUE
      )
    )
  ]

  ad_environment <- NULL

  for (candidate_file in ad_object_candidates) {
    candidate_environment <- load_to_environment(
      candidate_file
    )

    required_ad_objects <- c(
      "geneInfo_DLPFC",
      "geneInfo_ERC",
      "geneInfo_HIPPO"
    )

    if (all(
      required_ad_objects %in%
      ls(candidate_environment)
    )) {
      ad_environment <- candidate_environment
      break
    }
  }

  ftld_object_file <- file.path(
    project_root,
    "FTLD_geneInfoObjects.RData"
  )

  movement_object_file <- file.path(
    project_root,
    "MsaPspPdGeneInfos.RData"
  )

  if (is.null(ad_environment) ||
      !file.exists(ftld_object_file) ||
      !file.exists(movement_object_file)) {
    stop(
      paste0(
        "Network gene annotations could not be located.\n",
        "Provide the nine neutral gene-module CSV files in:\n",
        network_annotation_directory
      ),
      call. = FALSE
    )
  }

  ftld_environment <- load_to_environment(
    ftld_object_file
  )

  movement_environment <- load_to_environment(
    movement_object_file
  )

  required_ftld_objects <- c(
    "geneInfo_FTLD1",
    "geneInfo_FTLD2",
    "geneInfo_FTLD3"
  )

  required_movement_objects <- c(
    "MSAgeneInfo0K",
    "PDgeneInfo0K",
    "PSPgeneInfo0K"
  )

  if (!all(
    required_ftld_objects %in%
    ls(ftld_environment)
  )) {
    stop(
      "The FTLD working object lacks one or more required networks.",
      call. = FALSE
    )
  }

  if (!all(
    required_movement_objects %in%
    ls(movement_environment)
  )) {
    stop(
      "The movement-disorder working object lacks one or more required networks.",
      call. = FALSE
    )
  }

  working_networks <- list(
    AD_DLPFC =
      ad_environment[["geneInfo_DLPFC"]],
    AD_ERC =
      ad_environment[["geneInfo_ERC"]],
    AD_HIPPO =
      ad_environment[["geneInfo_HIPPO"]],
    FTLD1 =
      ftld_environment[["geneInfo_FTLD1"]],
    FTLD2 =
      ftld_environment[["geneInfo_FTLD2"]],
    FTLD3 =
      ftld_environment[["geneInfo_FTLD3"]],
    MSA =
      movement_environment[["MSAgeneInfo0K"]],
    PD =
      movement_environment[["PDgeneInfo0K"]],
    PSP =
      movement_environment[["PSPgeneInfo0K"]]
  )

  for (network_name in names(working_networks)) {
    neutral_annotation <- clean_network_annotation(
      working_networks[[network_name]]
    ) |>
      dplyr::rename(
        moduleColor = module
      )

    readr::write_csv(
      neutral_annotation,
      network_input_files[[network_name]]
    )

    message(
      "Created neutral network input: ",
      network_input_files[[network_name]]
    )
  }

  rm(
    ad_environment,
    ftld_environment,
    movement_environment,
    working_networks
  )
}

check_input_files(
  unname(network_input_files),
  description = "network gene-module annotation files"
)

# Read and validate the nine networks -------------------------------------------

network_annotations <- lapply(
  network_input_files,
  function(filename) {
    readr::read_csv(
      filename,
      show_col_types = FALSE
    ) |>
      clean_network_annotation()
  }
)

expected_module_counts <- c(
  AD_DLPFC = 16L,
  AD_ERC = 25L,
  AD_HIPPO = 18L,
  FTLD1 = 33L,
  FTLD2 = 49L,
  FTLD3 = 14L,
  MSA = 32L,
  PD = 32L,
  PSP = 32L
)

observed_module_counts <- vapply(
  network_annotations,
  function(annotation) {
    dplyr::n_distinct(annotation$module)
  },
  FUN.VALUE = integer(1)
)

incorrect_module_counts <- names(
  expected_module_counts
)[
  observed_module_counts[
    names(expected_module_counts)
  ] != expected_module_counts
]

if (length(incorrect_module_counts) > 0L) {
  stop(
    paste0(
      "Unexpected module counts were detected for: ",
      paste(incorrect_module_counts, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Fisher enrichment function ---------------------------------------------------

run_fisher_enrichment <- function(
    network_name,
    network_annotation,
    gene_set_name,
    curated_genes) {

  background_genes <- network_annotation |>
    dplyr::distinct(geneSymbol) |>
    dplyr::pull(geneSymbol)

  available_gene_set <- intersect(
    curated_genes,
    background_genes
  )

  network_modules <- sort(
    unique(network_annotation$module)
  )

  module_results <- purrr::map_dfr(
    network_modules,
    function(module_name) {
      module_genes <- network_annotation |>
        dplyr::filter(
          module == module_name
        ) |>
        dplyr::distinct(geneSymbol) |>
        dplyr::pull(geneSymbol)

      genes_in_set_and_module <- length(
        intersect(
          module_genes,
          available_gene_set
        )
      )

      genes_not_in_set_in_module <- length(
        setdiff(
          module_genes,
          available_gene_set
        )
      )

      genes_in_set_outside_module <- length(
        setdiff(
          available_gene_set,
          module_genes
        )
      )

      genes_not_in_set_outside_module <- length(
        setdiff(
          background_genes,
          union(
            module_genes,
            available_gene_set
          )
        )
      )

      contingency_table <- matrix(
        c(
          genes_in_set_and_module,
          genes_not_in_set_in_module,
          genes_in_set_outside_module,
          genes_not_in_set_outside_module
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          Module_location = c(
            "Inside_module",
            "Outside_module"
          ),
          Metal_set_membership = c(
            "Metal_gene",
            "Other_gene"
          )
        )
      )

      fisher_result <- stats::fisher.test(
        contingency_table,
        alternative = "two.sided"
      )

      tibble::tibble(
        network_id = network_name,
        module = module_name,
        gene_set = gene_set_name,
        overlap = genes_in_set_and_module,
        module_size = length(module_genes),
        network_background_size =
          length(background_genes),
        available_gene_set_size =
          length(available_gene_set),
        odds_ratio =
          unname(fisher_result$estimate),
        p_value =
          fisher_result$p.value
      )
    }
  )

  module_results |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        p_value,
        method = "BH"
      )
    )
}

# Run every network-gene-set combination ----------------------------------------

gene_set_labels <- c(
  copper = "Copper",
  broad_iron = "Broad iron",
  core_iron = "Core iron"
)

message("Running Fisher enrichment tests.")

all_module_enrichment_results <- purrr::imap_dfr(
  network_annotations,
  function(network_annotation, network_name) {
    purrr::imap_dfr(
      metal_gene_sets,
      function(curated_genes, internal_gene_set_name) {
        run_fisher_enrichment(
          network_name = network_name,
          network_annotation =
            network_annotation,
          gene_set_name =
            gene_set_labels[[
              internal_gene_set_name
            ]],
          curated_genes = curated_genes
        )
      }
    )
  }
) |>
  dplyr::arrange(
    network_id,
    gene_set,
    p_value
  )

# Restrict interpretation to disease-associated modules ------------------------

disease_module_lookup <- disease_associated_modules |>
  dplyr::select(
    network_id,
    module,
    module_id,
    disease_group,
    disease_or_comparison,
    dataset_or_network,
    source_study
  ) |>
  dplyr::distinct()

disease_module_enrichment_results <-
  all_module_enrichment_results |>
  dplyr::inner_join(
    disease_module_lookup,
    by = c(
      "network_id",
      "module"
    )
  ) |>
  dplyr::mutate(
    positive_enrichment =
      p_value < analysis_constants$nominal_p_threshold &
      odds_ratio > 1,
    FDR_significant =
      FDR < analysis_constants$fdr_threshold
  )

positive_enrichment_results <-
  disease_module_enrichment_results |>
  dplyr::filter(
    positive_enrichment
  )

# Validate final reported result identities -------------------------------------

expected_result_keys <- c(
  "AD_ERC|greenyellow|Copper",
  "FTLD1|red|Copper",
  "AD_ERC|darkred|Broad iron",
  "AD_HIPPO|yellow|Core iron",
  "FTLD2|darkmagenta|Core iron",
  "FTLD2|darkmagenta|Broad iron",
  "FTLD3|black|Core iron",
  "FTLD3|turquoise|Core iron",
  "FTLD3|black|Broad iron",
  "FTLD3|turquoise|Broad iron",
  "MSA|violet|Core iron",
  "PSP|pink|Core iron"
)

positive_enrichment_results <-
  positive_enrichment_results |>
  dplyr::mutate(
    result_key = paste(
      network_id,
      module,
      gene_set,
      sep = "|"
    ),
    result_order = match(
      result_key,
      expected_result_keys
    )
  )

observed_result_keys <-
  positive_enrichment_results$result_key

if (!setequal(
  observed_result_keys,
  expected_result_keys
)) {
  missing_results <- setdiff(
    expected_result_keys,
    observed_result_keys
  )

  unexpected_results <- setdiff(
    observed_result_keys,
    expected_result_keys
  )

  stop(
    paste0(
      "Enrichment results do not match the submitted analysis.\n",
      "Missing: ",
      paste(missing_results, collapse = ", "),
      "\nUnexpected: ",
      paste(unexpected_results, collapse = ", ")
    ),
    call. = FALSE
  )
}

positive_enrichment_results <-
  positive_enrichment_results |>
  dplyr::arrange(result_order)

if (nrow(positive_enrichment_results) != 12L) {
  stop(
    "Expected 12 positive module-gene-set enrichment results.",
    call. = FALSE
  )
}

metal_enriched_modules <- positive_enrichment_results |>
  dplyr::distinct(
    module_id,
    .keep_all = TRUE
  )

if (nrow(metal_enriched_modules) != 9L) {
  stop(
    "Expected nine metal-enriched disease-associated modules.",
    call. = FALSE
  )
}

# Validate key numerical results ------------------------------------------------

expected_statistics <- tibble::tribble(
  ~result_key, ~overlap, ~module_size, ~odds_ratio, ~p_value, ~FDR,

  "AD_ERC|greenyellow|Copper",
  6L, 863L, 3.65, 0.0107, 0.267,

  "FTLD1|red|Copper",
  7L, 1598L, 2.72, 0.0272, 0.448,

  "AD_ERC|darkred|Broad iron",
  24L, 1823L, 1.64, 0.0418, 0.467,

  "AD_HIPPO|yellow|Core iron",
  22L, 2433L, 1.66, 0.0445, 0.400,

  "FTLD2|darkmagenta|Core iron",
  17L, 1063L, 2.89, 0.000327, 0.0160,

  "FTLD2|darkmagenta|Broad iron",
  21L, 1063L, 2.35, 0.00120, 0.0587,

  "FTLD3|black|Core iron",
  47L, 5371L, 1.66, 0.0146, 0.106,

  "FTLD3|turquoise|Core iron",
  48L, 5506L, 1.67, 0.0151, 0.106,

  "FTLD3|black|Broad iron",
  63L, 5371L, 1.46, 0.0274, 0.243,

  "FTLD3|turquoise|Broad iron",
  64L, 5506L, 1.44, 0.0347, 0.243,

  "MSA|violet|Core iron",
  16L, 1461L, 1.85, 0.0362, 0.582,

  "PSP|pink|Core iron",
  18L, 1729L, 1.77, 0.0364, 0.582
)

numerical_validation <- positive_enrichment_results |>
  dplyr::select(
    result_key,
    observed_overlap = overlap,
    observed_module_size = module_size,
    observed_odds_ratio = odds_ratio,
    observed_p_value = p_value,
    observed_FDR = FDR
  ) |>
  dplyr::left_join(
    expected_statistics,
    by = "result_key"
  ) |>
  dplyr::mutate(
    matches_expected =
      observed_overlap == overlap &
      observed_module_size == module_size &
      abs(observed_odds_ratio - odds_ratio) < 0.01 &
      abs(observed_p_value - p_value) < 0.0001 &
      abs(observed_FDR - FDR) < 0.001
  )

if (!all(numerical_validation$matches_expected)) {
  stop(
    "One or more enrichment statistics differ from the submitted results.",
    call. = FALSE
  )
}

# Prepare Table 3A: copper ------------------------------------------------------

capitalise_module_label <- function(x) {
  paste0(toupper(substr(x, 1L, 1L)), substring(x, 2L))
}

table_3a <- positive_enrichment_results |>
  dplyr::filter(
    gene_set == "Copper"
  ) |>
  dplyr::mutate(
    No. = dplyr::row_number()
  ) |>
  dplyr::transmute(
    `No.` = No.,
    `Disease group` = disease_group,
    `Dataset/cohort` = dataset_or_network,
    Module = capitalise_module_label(module),
    `Gene set` = gene_set,
    Overlap = overlap,
    `Module size` = module_size,
    `Odds ratio` = odds_ratio,
    `P-value` = p_value,
    FDR = FDR
  )

# Prepare Table 3B: iron --------------------------------------------------------

table_3b <- positive_enrichment_results |>
  dplyr::filter(
    gene_set %in% c(
      "Broad iron",
      "Core iron"
    )
  ) |>
  dplyr::mutate(
    No. = dplyr::row_number()
  ) |>
  dplyr::transmute(
    `No.` = No.,
    `Disease group` = disease_group,
    `Dataset/cohort` = dataset_or_network,
    Module = capitalise_module_label(module),
    `Gene set` = gene_set,
    Overlap = overlap,
    `Module size` = module_size,
    `Odds ratio` = odds_ratio,
    `P-value` = p_value,
    FDR = FDR
  )

# Prepare Supplementary Table S4 -----------------------------------------------

supplementary_table_s4 <-
  positive_enrichment_results |>
  dplyr::mutate(
    No. = dplyr::row_number()
  ) |>
  dplyr::transmute(
    `No.` = No.,
    `Disease group` = disease_group,
    `Dataset/cohort` = dataset_or_network,
    Module = capitalise_module_label(module),
    `Gene set` = gene_set,
    Overlap = overlap,
    `Module size` = module_size,
    `Network background size` =
      network_background_size,
    `Available gene-set size` =
      available_gene_set_size,
    `Odds ratio` = odds_ratio,
    `P-value` = p_value,
    `FDR-adjusted P-value` = FDR
  )

# Figure 5 source data ----------------------------------------------------------

module_colour_values <- c(
  "AD ERC greenyellow" = "#66A61E",
  "AD ERC darkred" = "#9E2A2B",
  "AD HIPPO yellow" = "#D99A00",
  "FTLD1 red" = "#E53935",
  "FTLD2 darkmagenta" = "#9C0075",
  "FTLD3 black" = "#222222",
  "FTLD3 turquoise" = "#119C95",
  "MSA violet" = "#7655A5",
  "PSP pink" = "#D95F9F"
)

figure_5_data <- positive_enrichment_results |>
  dplyr::mutate(
    panel = dplyr::case_when(
      disease_or_comparison == "AD" ~
        "Alzheimer's disease",
      disease_or_comparison == "FTLD" ~
        "Frontotemporal lobar degeneration",
      TRUE ~ "Movement disorders"
    ),
    module_display = dplyr::case_when(
      network_id == "AD_ERC" ~
        paste("AD ERC", module),
      network_id == "AD_HIPPO" ~
        paste("AD HIPPO", module),
      network_id %in% c(
        "FTLD1",
        "FTLD2",
        "FTLD3",
        "MSA",
        "PSP"
      ) ~
        paste(network_id, module),
      TRUE ~
        paste(network_id, module)
    ),
    gene_set = factor(
      gene_set,
      levels = c(
        "Copper",
        "Broad iron",
        "Core iron"
      )
    ),
    minus_log10_p = -log10(p_value),
    FDR_significant =
      FDR < analysis_constants$fdr_threshold
  )

# Figure construction -----------------------------------------------------------

make_axis_labels <- function(module_names) {
  colours <- module_colour_values[module_names]

  labels <- paste0(
    "<span style='color:",
    colours,
    "'><b>",
    module_names,
    "</b></span>"
  )

  stats::setNames(
    labels,
    module_names
  )
}

make_enrichment_panel <- function(
    panel_data,
    panel_letter,
    panel_title,
    module_order_top_to_bottom) {

  panel_data <- panel_data |>
    dplyr::mutate(
      module_display = factor(
        module_display,
        levels = rev(
          module_order_top_to_bottom
        )
      )
    )

  axis_labels <- make_axis_labels(
    module_order_top_to_bottom
  )

  plot_object <- ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = gene_set,
      y = module_display
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = odds_ratio,
        fill = minus_log10_p
      ),
      shape = 21,
      colour = "#555555",
      stroke = 0.8
    ) +
    ggplot2::geom_point(
      data = panel_data |>
        dplyr::filter(FDR_significant),
      ggplot2::aes(
        size = odds_ratio,
        fill = minus_log10_p
      ),
      shape = 21,
      colour = "black",
      stroke = 2
    ) +
    ggplot2::geom_text(
      data = panel_data |>
        dplyr::filter(FDR_significant),
      ggplot2::aes(
        label = "*"
      ),
      nudge_x = 0.17,
      nudge_y = 0.08,
      size = 7,
      fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      labels = axis_labels,
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      name = paste0(
        "Odds ratio\n",
        "(enrichment strength)"
      ),
      range = c(5, 12),
      breaks = c(
        1.5,
        2.0,
        2.5,
        3.0,
        3.5
      ),
      limits = c(1.4, 3.7)
    ) +
    ggplot2::scale_fill_gradient(
      name = expression(
        -log[10](
          plain("nominal ") * italic(P)
        )
      ),
      low = "#FDE4CC",
      high = "#CF0000"
    ) +
    ggplot2::labs(
      title = paste(
        panel_letter,
        panel_title
      ),
      x = "Curated metal-associated gene set",
      y = NULL
    ) +
    ggplot2::theme_bw(
      base_size = 14
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 17,
        hjust = 0
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      axis.text.y = ggplot2::element_text(
        size = 13,
        hjust = 1,
        margin = ggplot2::margin(r = 10)
      ),
      axis.title.x = ggplot2::element_text(
        size = 13,
        margin = ggplot2::margin(t = 8)
      ),
      panel.grid.major.x =
        ggplot2::element_line(
          colour = "#E5E5E5"
        ),
      panel.grid.major.y =
        ggplot2::element_blank(),
      panel.grid.minor =
        ggplot2::element_blank(),
      panel.border =
        ggplot2::element_rect(
          colour = "#9E9E9E",
          linewidth = 0.8
        ),
      axis.line =
        ggplot2::element_line(
          colour = "black"
        ),
      legend.position = "bottom"
    )

  plot_object
}

figure_5_ad <- make_enrichment_panel(
  panel_data = figure_5_data |>
    dplyr::filter(
      panel == "Alzheimer's disease"
    ),
  panel_letter = "A",
  panel_title = "Alzheimer's disease",
  module_order_top_to_bottom = c(
    "AD ERC greenyellow",
    "AD ERC darkred",
    "AD HIPPO yellow"
  )
)

figure_5_ftld <- make_enrichment_panel(
  panel_data = figure_5_data |>
    dplyr::filter(
      panel ==
        "Frontotemporal lobar degeneration"
    ),
  panel_letter = "B",
  panel_title =
    "Frontotemporal lobar degeneration",
  module_order_top_to_bottom = c(
    "FTLD1 red",
    "FTLD2 darkmagenta",
    "FTLD3 black",
    "FTLD3 turquoise"
  )
)

figure_5_movement <- make_enrichment_panel(
  panel_data = figure_5_data |>
    dplyr::filter(
      panel == "Movement disorders"
    ),
  panel_letter = "C",
  panel_title = "Movement disorders",
  module_order_top_to_bottom = c(
    "MSA violet",
    "PSP pink"
  )
)

figure_5 <- (
  figure_5_ad /
    figure_5_ftld /
    figure_5_movement
) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(
      1.0,
      1.25,
      0.75
    )
  ) +
  patchwork::plot_annotation(
    title =
      "Copper- and iron-gene enrichment across diseases",
    caption =
      "* Black outline indicates FDR-adjusted P < 0.05.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 20,
        hjust = 0
      ),
      plot.caption = ggplot2::element_text(
        size = 10,
        hjust = 1
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

# Save tables and source data ---------------------------------------------------

write_output_csv(
  all_module_enrichment_results,
  filename =
    "all_network_module_metal_enrichment_results.csv"
)

write_output_csv(
  disease_module_enrichment_results,
  filename =
    "disease_associated_module_metal_enrichment_results.csv"
)

write_output_csv(
  positive_enrichment_results |>
    dplyr::select(-result_order),
  filename =
    "positive_metal_enrichment_results.csv"
)

write_output_csv(
  table_3a,
  filename =
    "Table_3A_copper_gene_set_enrichment.csv",
  output_directory = paths$tables
)

write_output_csv(
  table_3b,
  filename =
    "Table_3B_iron_gene_set_enrichment.csv",
  output_directory = paths$tables
)

write_output_csv(
  supplementary_table_s4,
  filename =
    "Supplementary_Table_S4_module_level_enrichment.csv",
  output_directory =
    paths$supplementary_tables
)

write_output_csv(
  figure_5_data |>
    dplyr::select(
      panel,
      module_display,
      gene_set,
      overlap,
      module_size,
      odds_ratio,
      p_value,
      FDR,
      minus_log10_p,
      FDR_significant
    ),
  filename =
    "Figure_5_source_data.csv",
  output_directory = paths$figure_data
)

saveRDS(
  positive_enrichment_results,
  file.path(
    paths$results,
    "metal_enriched_modules_and_gene_sets.rds"
  )
)

# Save Figure 5 -----------------------------------------------------------------

save_thesis_figure(
  plot = figure_5,
  filename =
    "Figure_5_metal_gene_set_enrichment",
  width = 10,
  height = 14,
  output_directory = paths$figures
)

ggplot2::ggsave(
  filename = file.path(
    paths$figures,
    "Figure_5_metal_gene_set_enrichment.tiff"
  ),
  plot = figure_5,
  width = 10,
  height = 14,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Analysis summary --------------------------------------------------------------

enrichment_summary <- tibble::tibble(
  measure = c(
    "Disease-associated modules tested",
    "Metal-enriched modules",
    "Positive module-gene-set results",
    "Copper-enriched modules",
    "Iron-enriched modules",
    "FDR-significant enrichment results"
  ),
  value = c(
    dplyr::n_distinct(
      disease_associated_modules$module_id
    ),
    dplyr::n_distinct(
      positive_enrichment_results$module_id
    ),
    nrow(positive_enrichment_results),
    positive_enrichment_results |>
      dplyr::filter(gene_set == "Copper") |>
      dplyr::distinct(module_id) |>
      nrow(),
    positive_enrichment_results |>
      dplyr::filter(
        gene_set %in% c(
          "Broad iron",
          "Core iron"
        )
      ) |>
      dplyr::distinct(module_id) |>
      nrow(),
    sum(
      positive_enrichment_results$FDR <
        analysis_constants$fdr_threshold
    )
  )
)

write_output_csv(
  enrichment_summary,
  filename =
    "metal_gene_set_enrichment_summary.csv"
)

capture.output(
  sessionInfo(),
  file = file.path(
    paths$logs,
    "04_metal_gene_set_enrichment_session_info.txt"
  )
)

message("Metal gene-set enrichment analysis completed successfully.")
message(
  "Metal-enriched disease-associated modules: ",
  dplyr::n_distinct(
    positive_enrichment_results$module_id
  )
)
message(
  "Positive module-gene-set results: ",
  nrow(positive_enrichment_results)
)
message(
  "FDR-significant results: ",
  sum(
    positive_enrichment_results$FDR <
      analysis_constants$fdr_threshold
  )
)
