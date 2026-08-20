
# Functional enrichment of corrected metal-enriched modules and recurrent genes
#
# Module analysis:
#   Each complete module gene list is tested against all unique annotated genes
#   represented in the corresponding published network.
#
# Recurrent-set analysis:
#   Recurrent copper, broad-iron and core-iron genes are tested separately against
#   their complete curated metal-gene sets (47, 199 and 136 genes, respectively).
#
# g:Profiler settings:
#   Homo sapiens; unranked queries; GO Biological Process, Reactome and
#   WikiPathways; custom annotated backgrounds; FDR-adjusted P < 0.05.
#
# Generates Figure 11, Table 5 and Supplementary Table S6.

required_packages <- c(
  "dplyr", "ggplot2", "gprofiler2", "patchwork", "purrr",
  "readr", "scales", "stringr", "tibble", "tidyr"
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
  library(gprofiler2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
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
  source_data = file.path(project_root, "outputs", "figure_source_data"),
  logs = file.path(project_root, "outputs", "logs")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Set REUSE_GPROFILER_RESULTS=1 to regenerate figures and tables from a completed
# saved run without submitting the nine module and three recurrent queries again.
reuse_saved_results <- identical(Sys.getenv("REUSE_GPROFILER_RESULTS"), "1")

# -----------------------------------------------------------------------------
# 2. Corrected nine-module specification
# -----------------------------------------------------------------------------

module_specification <- tribble(
  ~module_order, ~module_id,             ~module_label,           ~object_name,      ~annotation_file,              ~module_colour, ~metal_enrichment,       ~panel,    ~module_hex,
  1L,            "AD_ERC_greenyellow",  "AD ERC greenyellow",   "geneInfo_ERC",    "AD_ERC_gene_modules.csv",    "greenyellow", "Copper",                "Copper", "#7CAE00",
  2L,            "AD_ERC_darkred",      "AD ERC darkred",       "geneInfo_ERC",    "AD_ERC_gene_modules.csv",    "darkred",     "Broad iron",            "Iron",   "#8B1A1A",
  3L,            "AD_HIPPO_yellow",     "AD HIPPO yellow",      "geneInfo_HIPPO",  "AD_HIPPO_gene_modules.csv",  "yellow",      "Core iron",             "Iron",   "#E6B800",
  4L,            "FTLD1_red",           "FTLD1 red",            "geneInfo_FTLD1",  "FTLD1_gene_modules.csv",     "red",         "Copper",                "Copper", "#D73027",
  5L,            "FTLD2_darkmagenta",   "FTLD2 darkmagenta",    "geneInfo_FTLD2",  "FTLD2_gene_modules.csv",     "darkmagenta", "Broad and core iron", "Iron",   "#8B006B",
  6L,            "FTLD3_black",         "FTLD3 black",          "geneInfo_FTLD3",  "FTLD3_gene_modules.csv",     "black",       "Broad and core iron", "Iron",   "#222222",
  7L,            "FTLD3_turquoise",     "FTLD3 turquoise",      "geneInfo_FTLD3",  "FTLD3_gene_modules.csv",     "turquoise",   "Broad and core iron", "Iron",   "#109D92",
  8L,            "MSA_violet",          "MSA violet",           "MSAgeneInfo0K",   "MSA_gene_modules.csv",       "violet",      "Core iron",             "Iron",   "#7452A3",
  9L,            "PSP_pink",            "PSP pink",             "PSPgeneInfo0K",   "PSP_gene_modules.csv",       "pink",        "Core iron",             "Iron",   "#D95F9D"
)

clean_gene_symbols <- function(x) {
  x <- str_to_upper(str_trim(as.character(x)))
  sort(unique(x[!is.na(x) & x != "" & x != "NA"]))
}

# -----------------------------------------------------------------------------
# 3. Load published network gene-information objects
# -----------------------------------------------------------------------------

combined_network_rds <- file.path(
  project_root, "data", "published_results", "network_gene_information.rds"
)
annotation_directory <- file.path(
  project_root, "data", "published_results", "network_gene_annotations"
)
annotation_specification <- module_specification |>
  distinct(object_name, annotation_file)
annotation_paths <- file.path(
  annotation_directory, annotation_specification$annotation_file
)

if (all(file.exists(annotation_paths))) {
  message("Using neutral published-network annotation CSV files.")
  network_objects <- setNames(
    map(annotation_paths, ~ read_csv(.x, show_col_types = FALSE)),
    annotation_specification$object_name
  )
} else if (file.exists(combined_network_rds)) {
  network_objects <- readRDS(combined_network_rds)
} else {
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

missing_objects <- setdiff(
  unique(module_specification$object_name),
  names(network_objects)
)
if (length(missing_objects) > 0L) {
  stop(
    "Missing network gene-information objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

extract_module_input <- function(specification_row) {
  specification_row <- as.list(specification_row)
  network_data <- network_objects[[specification_row$object_name]]

  if (!all(c("geneSymbol", "moduleColor") %in% names(network_data))) {
    stop(specification_row$object_name, " lacks geneSymbol or moduleColor.")
  }

  background_genes <- clean_gene_symbols(network_data$geneSymbol)
  module_genes <- network_data |>
    filter(
      str_to_lower(str_trim(as.character(moduleColor))) ==
        specification_row$module_colour
    ) |>
    pull(geneSymbol) |>
    clean_gene_symbols()

  if (!all(module_genes %in% background_genes)) {
    stop(specification_row$module_id, " is not a subset of its network background.")
  }

  list(
    summary = tibble(
      module_order = specification_row$module_order,
      module_id = specification_row$module_id,
      module_label = specification_row$module_label,
      metal_enrichment = specification_row$metal_enrichment,
      module_genes = length(module_genes),
      background_genes = length(background_genes)
    ),
    query = module_genes,
    background = background_genes
  )
}

module_inputs <- map(
  seq_len(nrow(module_specification)),
  ~ extract_module_input(module_specification[.x, , drop = FALSE])
)
names(module_inputs) <- module_specification$module_id

module_input_summary <- map_dfr(module_inputs, "summary") |>
  arrange(module_order)

expected_module_sizes <- c(863L, 1823L, 2433L, 1598L, 1063L, 5371L, 5506L, 1461L, 1729L)
expected_background_sizes <- c(14429L, 14429L, 13606L, 15242L, 15193L, 15037L, 15037L, 15579L, 15579L)

if (!identical(module_input_summary$module_genes, expected_module_sizes) ||
    !identical(module_input_summary$background_genes, expected_background_sizes)) {
  print(module_input_summary)
  stop(
    "Module or background sizes differ from the submitted analysis. ",
    "Confirm that the final published network objects are being used."
  )
}

write_csv(
  module_input_summary,
  file.path(output_dirs$derived_data, "nine_module_gprofiler_input_summary.csv")
)

# -----------------------------------------------------------------------------
# 4. Load recurrent sets and their complete curated backgrounds
# -----------------------------------------------------------------------------

gene_set_rds_candidates <- c(
  file.path(project_root, "outputs", "derived_data", "curated_metal_gene_sets.rds"),
  file.path(project_root, "data", "derived", "curated_metal_gene_sets.rds"),
  file.path(project_root, "gene_lists", "curated_metal_gene_sets.rds")
)
gene_set_rds <- gene_set_rds_candidates[file.exists(gene_set_rds_candidates)][1]
if (is.na(gene_set_rds)) {
  stop("Run R/02_curate_metal_gene_sets.R before functional enrichment.")
}

curated_sets <- readRDS(gene_set_rds)
required_sets <- c("copper", "broad_iron", "core_iron")
if (!all(required_sets %in% names(curated_sets))) {
  stop("The curated gene-set RDS must contain copper, broad_iron and core_iron.")
}
curated_sets <- lapply(curated_sets[required_sets], clean_gene_symbols)

recurrent_query_file <- file.path(
  project_root, "outputs", "derived_data", "recurrent_gene_queries.csv"
)
if (!file.exists(recurrent_query_file)) {
  stop("Run R/07_recurrence_analysis_9_modules.R before functional enrichment.")
}

recurrent_query_table <- read_csv(recurrent_query_file, show_col_types = FALSE)
if (!all(c("query", "gene_symbol") %in% names(recurrent_query_table))) {
  stop("The recurrent query file does not contain query and gene_symbol columns.")
}

recurrent_queries <- split(
  recurrent_query_table$gene_symbol,
  recurrent_query_table$query
) |>
  lapply(clean_gene_symbols)

recurrent_query_names <- c(
  "Recurrent_copper", "Recurrent_broad_iron", "Recurrent_core_iron"
)
if (!all(recurrent_query_names %in% names(recurrent_queries))) {
  stop("One or more corrected recurrent query sets are missing.")
}
recurrent_queries <- recurrent_queries[recurrent_query_names]

recurrent_backgrounds <- list(
  Recurrent_copper = curated_sets$copper,
  Recurrent_broad_iron = curated_sets$broad_iron,
  Recurrent_core_iron = curated_sets$core_iron
)

expected_query_counts <- c(13L, 77L, 57L)
expected_background_counts <- c(47L, 199L, 136L)
observed_query_counts <- map_int(recurrent_queries, length)
observed_background_counts <- map_int(recurrent_backgrounds, length)

if (!identical(unname(observed_query_counts), expected_query_counts) ||
    !identical(unname(observed_background_counts), expected_background_counts)) {
  stop("Recurrent query or curated-background counts do not match the final analysis.")
}
if (!all(map2_lgl(
  recurrent_queries,
  recurrent_backgrounds,
  ~ all(.x %in% .y)
))) {
  stop("Every recurrent query must be a subset of its curated metal background.")
}

# -----------------------------------------------------------------------------
# 5. Run g:Profiler
# -----------------------------------------------------------------------------

empty_gost_result <- function() {
  tibble(
    source = character(),
    term_id = character(),
    term_name = character(),
    p_value = double(),
    term_size = integer(),
    query_size = integer(),
    intersection_size = integer(),
    effective_domain_size = integer(),
    precision = double(),
    recall = double(),
    intersection = character()
  )
}

flatten_list_columns <- function(data) {
  list_columns <- names(data)[vapply(data, is.list, logical(1))]
  for (column in list_columns) {
    data[[column]] <- vapply(
      data[[column]],
      function(value) paste(as.character(value), collapse = "; "),
      character(1)
    )
  }
  data
}

run_gost_query <- function(query, background, query_name, analysis_type) {
  query <- intersect(clean_gene_symbols(query), clean_gene_symbols(background))
  background <- clean_gene_symbols(background)

  message(
    "Running ", query_name, ": ", length(query),
    " query genes against ", length(background), " background genes"
  )

  result <- tryCatch(
    gost(
      query = query,
      organism = "hsapiens",
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      user_threshold = 0.05,
      correction_method = "fdr",
      sources = c("GO:BP", "REAC", "WP"),
      custom_bg = background,
      domain_scope = "custom_annotated",
      evcodes = TRUE,
      highlight = TRUE
    ),
    error = function(error) {
      stop("g:Profiler failed for ", query_name, ": ", conditionMessage(error))
    }
  )

  if (is.null(result) || is.null(result$result) || nrow(result$result) == 0L) {
    return(empty_gost_result() |>
      mutate(
        analysis_type = analysis_type,
        query_name = query_name,
        query_genes_n = length(query),
        background_genes_n = length(background),
        .before = 1
      ))
  }

  as_tibble(result$result) |>
    flatten_list_columns() |>
    mutate(
      analysis_type = analysis_type,
      query_name = query_name,
      query_genes_n = length(query),
      background_genes_n = length(background),
      .before = 1
    )
}

module_results_rds <- file.path(
  output_dirs$derived_data, "nine_module_gprofiler_significant_results.rds"
)
recurrent_results_rds <- file.path(
  output_dirs$derived_data, "recurrent_sets_gprofiler_significant_results.rds"
)

module_result_candidates <- c(
  module_results_rds,
  file.path(
    project_root, "data", "derived",
    "submitted_nine_module_gprofiler_significant_results.rds"
  )
)
recurrent_result_candidates <- c(
  recurrent_results_rds,
  file.path(
    project_root, "data", "derived",
    "submitted_recurrent_sets_gprofiler_significant_results.rds"
  )
)
saved_module_results <- module_result_candidates[file.exists(module_result_candidates)][1]
saved_recurrent_results <- recurrent_result_candidates[file.exists(recurrent_result_candidates)][1]

# Earlier completed runs used descriptive module labels rather than the stable
# module identifiers used by this repository. Standardise either schema before
# rebuilding the final figures and tables.
standardise_saved_module_results <- function(data) {
  data <- as_tibble(data)
  if (!"query_name" %in% names(data)) {
    if (!"module" %in% names(data)) {
      stop("Saved module enrichment results lack query_name or module.")
    }
    data <- data |>
      left_join(
        module_specification |> select(module_id, module_label),
        by = c("module" = "module_label")
      ) |>
      mutate(
        analysis_type = "Metal-enriched module",
        query_name = module_id,
        query_genes_n = input_genes,
        background_genes_n = custom_background_genes,
        .before = 1
      ) |>
      select(-module_id)
  }
  data
}

standardise_saved_recurrent_results <- function(data) {
  data <- as_tibble(data)
  rows <- nrow(data)
  if (!"query_name" %in% names(data)) {
    query_name <- if ("analysis" %in% names(data)) {
      recode(
        as.character(data$analysis),
        "Copper" = "recurrent_copper",
        "Broad iron" = "recurrent_broad_iron",
        "Core iron" = "recurrent_core_iron",
        .default = as.character(data$analysis)
      )
    } else {
      rep(NA_character_, rows)
    }
    data$query_name <- query_name
  }
  if (!"analysis_type" %in% names(data)) {
    data$analysis_type <- rep("Recurrent metal-gene set", rows)
  }
  if (!"query_genes_n" %in% names(data)) {
    data$query_genes_n <- if ("query_genes_n" %in% names(data)) {
      data$query_genes_n
    } else rep(NA_integer_, rows)
  }
  if (!"background_genes_n" %in% names(data)) {
    data$background_genes_n <- if ("background_genes_n" %in% names(data)) {
      data$background_genes_n
    } else rep(NA_integer_, rows)
  }
  data
}

if (reuse_saved_results &&
    !is.na(saved_module_results) &&
    !is.na(saved_recurrent_results)) {
  message("Reusing the completed saved g:Profiler results.")
  module_results <- readRDS(saved_module_results) |>
    standardise_saved_module_results()
  recurrent_results <- readRDS(saved_recurrent_results) |>
    standardise_saved_recurrent_results()
} else {
  module_results <- imap_dfr(
    module_inputs,
    ~ run_gost_query(
      query = .x$query,
      background = .x$background,
      query_name = .y,
      analysis_type = "Metal-enriched module"
    )
  )

  recurrent_results <- imap_dfr(
    recurrent_queries,
    ~ run_gost_query(
      query = .x,
      background = recurrent_backgrounds[[.y]],
      query_name = .y,
      analysis_type = "Recurrent metal-gene set"
    )
  )

  saveRDS(module_results, module_results_rds)
  saveRDS(recurrent_results, recurrent_results_rds)
}

# Save canonical repository-schema copies even when a submitted legacy object
# was used as the input.
saveRDS(module_results, module_results_rds)
saveRDS(recurrent_results, recurrent_results_rds)

write_csv(
  module_results,
  file.path(output_dirs$derived_data, "nine_module_gprofiler_significant_results.csv")
)
write_csv(
  recurrent_results,
  file.path(output_dirs$derived_data, "recurrent_sets_gprofiler_significant_results.csv")
)

# The submitted run yielded 771 significant module–term results from eight of
# nine modules. g:Profiler is a live resource, so later database releases may
# alter exact counts; differences are reported rather than silently overwritten.
module_result_count <- nrow(module_results)
modules_with_results <- n_distinct(module_results$query_name)
if (module_result_count != 771L || modules_with_results != 8L) {
  warning(
    "The current g:Profiler release returned ", module_result_count,
    " significant module–term results across ", modules_with_results,
    " modules; the submitted analysis returned 771 results across eight modules."
  )
}

# -----------------------------------------------------------------------------
# 6. Query summaries and complete Supplementary Table S6
# -----------------------------------------------------------------------------

module_summary <- module_input_summary |>
  left_join(
    module_results |>
      count(query_name, name = "significant_terms"),
    by = c("module_id" = "query_name")
  ) |>
  mutate(
    significant_terms = replace_na(significant_terms, 0L),
    outcome = if_else(
      significant_terms == 0L,
      "No term reached FDR-adjusted P < 0.05",
      paste0(significant_terms, " significant term(s)")
    )
  )

recurrent_summary <- tibble(
  query_name = names(recurrent_queries),
  recurrent_set = c("Copper", "Broad iron", "Core iron"),
  query_genes = unname(observed_query_counts),
  background_genes = unname(observed_background_counts)
) |>
  left_join(
    recurrent_results |>
      count(query_name, name = "significant_terms"),
    by = "query_name"
  ) |>
  mutate(
    significant_terms = replace_na(significant_terms, 0L),
    outcome = if_else(
      significant_terms == 0L,
      "No GO Biological Process, Reactome or WikiPathways term reached FDR-adjusted P < 0.05",
      paste0(significant_terms, " significant term(s)")
    )
  )

write_csv(
  module_summary,
  file.path(output_dirs$tables, "functional_enrichment_module_query_summary.csv")
)
write_csv(
  recurrent_summary,
  file.path(output_dirs$tables, "functional_enrichment_recurrent_query_summary.csv")
)

supplementary_table_s6 <- bind_rows(module_results, recurrent_results) |>
  select(
    analysis_type, query_name, query_genes_n, background_genes_n,
    source, term_id, term_name, p_value, term_size, query_size,
    intersection_size, effective_domain_size, precision, recall,
    intersection, everything()
  )

write_csv(
  supplementary_table_s6,
  file.path(
    output_dirs$supplementary,
    "Supplementary_Table_S6_complete_functional_enrichment_results.csv"
  )
)

# -----------------------------------------------------------------------------
# 7. Representative terms for Figure 11
# -----------------------------------------------------------------------------

# These terms were selected from the complete significant output to show distinct
# biological themes while avoiding extensive display of closely related parent
# and child terms. Selection changes only the visual summary; Supplementary Table
# S6 retains every significant result.
representative_term_key <- tribble(
  ~panel,    ~term_order, ~term_name,                                               ~source,
  "Copper", 1L,          "cellular component organization",                       "GO:BP",
  "Iron",   1L,          "nervous system development",                            "GO:BP",
  "Iron",   2L,          "generation of neurons",                                 "GO:BP",
  "Iron",   3L,          "cell morphogenesis involved in neuron differentiation", "GO:BP",
  "Iron",   4L,          "nucleic acid metabolic process",                        "GO:BP",
  "Iron",   5L,          "signal transduction in response to DNA damage",         "GO:BP",
  "Iron",   6L,          "Metabolism of RNA",                                     "REAC",
  "Iron",   7L,          "Cell Cycle, Mitotic",                                   "REAC",
  "Iron",   8L,          "Nucleotide Excision Repair",                            "REAC",
  "Iron",   9L,          "Cellular responses to stress",                          "REAC",
  "Iron",   10L,         "regulation of TORC1 signaling",                         "GO:BP",
  "Iron",   11L,         "cell death",                                             "GO:BP",
  "Iron",   12L,         "signal transduction by p53 class mediator",              "GO:BP",
  "Iron",   13L,         "ATP biosynthetic process",                               "GO:BP",
  "Iron",   14L,         "EGF EGFR signaling",                                     "WP",
  "Iron",   15L,         "Microtubule cytoskeleton regulation",                    "WP",
  "Iron",   16L,         "regulation of programmed cell death",                    "GO:BP"
)

representative_results <- module_results |>
  left_join(
    module_specification |>
      select(module_id, panel, module_order, module_label, module_hex),
    by = c("query_name" = "module_id")
  ) |>
  inner_join(
    representative_term_key,
    by = c("panel", "term_name", "source")
  ) |>
  group_by(
    panel, term_order, term_name, source,
    query_name, module_order, module_label, module_hex
  ) |>
  summarise(
    adjusted_p = min(p_value),
    minus_log10_fdr = max(-log10(p_value)),
    intersection_size = max(intersection_size),
    .groups = "drop"
  )

missing_representative_terms <- representative_term_key |>
  anti_join(
    representative_results |>
      distinct(panel, term_name, source),
    by = c("panel", "term_name", "source")
  )
if (nrow(missing_representative_terms) > 0L) {
  warning(
    "Some prespecified representative terms were absent from the current ",
    "g:Profiler release; Figure 11 will display the available terms."
  )
}

source_display <- c("GO:BP" = "GO BP", "REAC" = "Reactome", "WP" = "WikiPathways")
representative_term_key <- representative_term_key |>
  mutate(
    term_label = str_wrap(
      paste0(term_name, " [", source_display[source], "]"),
      width = 47
    )
  )

prepare_figure_panel <- function(panel_name) {
  panel_modules <- module_specification |>
    filter(panel == panel_name) |>
    arrange(module_order)
  panel_terms <- representative_term_key |>
    filter(panel == panel_name) |>
    arrange(term_order)

  result_data <- representative_results |>
    filter(panel == panel_name) |>
    left_join(
      panel_terms |>
        select(term_name, source, term_label),
      by = c("term_name", "source")
    ) |>
    mutate(
      module_label = factor(module_label, levels = panel_modules$module_label),
      term_label = factor(term_label, levels = rev(panel_terms$term_label))
    )

  background_grid <- expand_grid(
    module_label = factor(
      panel_modules$module_label,
      levels = panel_modules$module_label
    ),
    term_label = factor(
      panel_terms$term_label,
      levels = rev(panel_terms$term_label)
    )
  )

  module_strip <- panel_modules |>
    mutate(
      module_label = factor(module_label, levels = panel_modules$module_label),
      y = 1
    )

  list(results = result_data, background = background_grid, strip = module_strip)
}

make_dotplot <- function(prepared, title, subtitle, low_colour, high_colour, max_size) {
  ggplot() +
    geom_tile(
      data = prepared$background,
      aes(module_label, term_label),
      fill = "#F5F5F5",
      colour = "white",
      linewidth = 0.6
    ) +
    geom_point(
      data = prepared$results,
      aes(
        module_label, term_label,
        size = intersection_size,
        fill = minus_log10_fdr
      ),
      shape = 21,
      colour = "#333333",
      stroke = 0.45,
      alpha = 0.95
    ) +
    scale_fill_gradient(
      low = low_colour,
      high = high_colour,
      name = expression(-log[10]("FDR-adjusted P")),
      labels = label_number(accuracy = 0.1)
    ) +
    scale_size_area(max_size = max_size, name = "Intersecting\ngenes") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 17, face = "bold", colour = "#202020"),
      plot.subtitle = element_text(
        size = 11.5, colour = "#4A4A4A", margin = margin(b = 8)
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10.5, colour = "#252525", lineheight = 0.95),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9.5),
      plot.margin = margin(8, 18, 8, 8)
    )
}

make_module_strip <- function(prepared) {
  ggplot(prepared$strip, aes(module_label, y)) +
    geom_point(
      aes(fill = module_hex),
      shape = 22,
      size = 7,
      colour = "#333333",
      stroke = 0.7
    ) +
    scale_fill_identity() +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(0.75, 1.25), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      axis.text.x = element_text(
        angle = 38, hjust = 1, vjust = 1,
        size = 10, face = "bold", colour = "#252525"
      ),
      plot.margin = margin(0, 117, 36, 8)
    )
}

copper_panel <- prepare_figure_panel("Copper")
iron_panel <- prepare_figure_panel("Iron")

figure_11 <- (
  make_dotplot(
    copper_panel,
    "A  Copper-enriched co-methylation modules",
    "Selected non-redundant significant term; blank cells indicate no significant result",
    "#FEE8C8", "#B33B00", 12
  ) /
    make_module_strip(copper_panel) /
    make_dotplot(
      iron_panel,
      "B  Iron-enriched co-methylation modules",
      "Selected non-redundant GO Biological Process, Reactome and WikiPathways terms",
      "#DEEBF7", "#084594", 10
    ) /
    make_module_strip(iron_panel)
) + plot_layout(heights = c(1.1, 0.22, 3.8, 0.22))

ggsave(
  file.path(output_dirs$figures, "Figure_11_functional_enrichment.png"),
  figure_11, width = 16, height = 17, dpi = 600, bg = "white"
)
ggsave(
  file.path(output_dirs$figures, "Figure_11_functional_enrichment.tiff"),
  figure_11, width = 16, height = 17, dpi = 600,
  compression = "lzw", bg = "white"
)
ggsave(
  file.path(output_dirs$figures, "Figure_11_functional_enrichment.pdf"),
  figure_11, width = 16, height = 17, device = grDevices::pdf,
  useDingbats = FALSE, bg = "white"
)

write_csv(
  representative_results,
  file.path(output_dirs$source_data, "Figure_11_source_data.csv")
)

# -----------------------------------------------------------------------------
# 8. Table 5 summaries
# -----------------------------------------------------------------------------

source_counts <- module_results |>
  count(query_name, source, name = "significant_terms") |>
  mutate(
    source = recode(
      source,
      "GO:BP" = "GO_BP",
      "REAC" = "Reactome",
      "WP" = "WikiPathways"
    )
  ) |>
  pivot_wider(
    names_from = source,
    values_from = significant_terms,
    values_fill = 0
  )

for (column in c("GO_BP", "Reactome", "WikiPathways")) {
  if (!column %in% names(source_counts)) source_counts[[column]] <- 0L
}

strongest_results <- module_results |>
  group_by(query_name) |>
  summarise(strongest_fdr = min(p_value), .groups = "drop")

representative_term_summary <- representative_results |>
  mutate(
    source_label = recode(
      source,
      "GO:BP" = "GO BP",
      "REAC" = "Reactome",
      "WP" = "WikiPathways"
    ),
    term_with_source = paste0(term_name, " [", source_label, "]")
  ) |>
  arrange(module_order, adjusted_p, desc(intersection_size)) |>
  group_by(query_name) |>
  summarise(
    representative_terms = paste(term_with_source, collapse = "; "),
    .groups = "drop"
  )

table_5a <- module_input_summary |>
  left_join(source_counts, by = c("module_id" = "query_name")) |>
  left_join(strongest_results, by = c("module_id" = "query_name")) |>
  left_join(representative_term_summary, by = c("module_id" = "query_name")) |>
  mutate(
    across(c(GO_BP, Reactome, WikiPathways), ~ replace_na(.x, 0L)),
    total_significant_terms = GO_BP + Reactome + WikiPathways,
    representative_terms = case_when(
      total_significant_terms == 0L ~ "No term reached FDR-adjusted P < 0.05",
      is.na(representative_terms) ~ "See Supplementary Table S6",
      TRUE ~ representative_terms
    )
  ) |>
  arrange(module_order) |>
  transmute(
    Module = module_label,
    `Metal enrichment` = metal_enrichment,
    `Module genes` = module_genes,
    `GO BP terms` = GO_BP,
    `Reactome pathways` = Reactome,
    `WikiPathways terms` = WikiPathways,
    `Total significant` = total_significant_terms,
    `Strongest FDR-adjusted P` = strongest_fdr,
    `Representative terms` = representative_terms
  )

table_5b <- recurrent_summary |>
  transmute(
    `Recurrent set` = recurrent_set,
    `Query genes` = query_genes,
    `Curated metal-set background` = background_genes,
    `Significant terms` = significant_terms,
    Outcome = outcome
  )

write_csv(table_5a, file.path(output_dirs$tables, "Table_5A_module_functional_enrichment.csv"))
write_csv(table_5b, file.path(output_dirs$tables, "Table_5B_recurrent_set_functional_enrichment.csv"))

# Optional editable Word version when the document packages are available.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  table_5a_word <- table_5a |>
    mutate(
      `Strongest FDR-adjusted P` = if_else(
        is.na(`Strongest FDR-adjusted P`),
        "—",
        formatC(`Strongest FDR-adjusted P`, format = "g", digits = 3)
      )
    )

  ft_a <- flextable::flextable(table_5a_word) |>
    flextable::theme_booktabs() |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::autofit()
  ft_b <- flextable::flextable(table_5b) |>
    flextable::theme_booktabs() |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::autofit()

  document <- officer::read_docx() |>
    officer::body_add_par(
      "Table 5. Functional enrichment of corrected metal-enriched modules and recurrent gene sets.",
      style = "Table Caption"
    ) |>
    officer::body_add_par("A. Module-level functional enrichment") |>
    flextable::body_add_flextable(ft_a) |>
    officer::body_add_par("B. Recurrent gene-set functional enrichment") |>
    flextable::body_add_flextable(ft_b)

  print(document, target = file.path(output_dirs$tables, "Table_5_functional_enrichment.docx"))
}

# -----------------------------------------------------------------------------
# 9. Record settings for reproducibility
# -----------------------------------------------------------------------------

settings <- c(
  "Functional-enrichment analysis settings",
  paste("Run date:", Sys.Date()),
  "Organism: Homo sapiens (hsapiens)",
  "Query type: unranked",
  "Sources: GO Biological Process (GO:BP), Reactome (REAC), WikiPathways (WP)",
  "Correction: false discovery rate (FDR)",
  "Threshold: FDR-adjusted P < 0.05",
  "Domain scope: custom annotated background",
  "Module background: all unique annotated genes in the corresponding network",
  "Recurrent background: complete corresponding curated metal-gene set",
  "Submitted module result: 771 significant module-term results across eight of nine modules",
  "Submitted recurrent result: zero significant terms for all three recurrent sets"
)
writeLines(
  settings,
  file.path(output_dirs$logs, "functional_enrichment_settings.txt")
)

message("R/10_functional_enrichment_gprofiler.R completed successfully.")
message(
  "Current module result: ", nrow(module_results),
  " significant terms across ", n_distinct(module_results$query_name), " modules."
)
message(
  "Current recurrent-set result: ", nrow(recurrent_results),
  " significant terms. Zero-result queries remain explicitly reported in Table 5B."
)
