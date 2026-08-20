
# Recurrence analysis across the corrected nine metal-enriched modules
#
# Recurrence is defined as representation of the same gene in at least two
# distinct retained module–disease comparisons. It is different from gene-set
# enrichment, which tests whether an entire curated metal set is overrepresented
# within one module.
#
# Reproduces the candidate totals and Supplementary Table S5 and prepares the
# source tables used for Figures 7–10 and Supplementary Figure S1.

required_packages <- c("dplyr", "readr", "stringr", "tibble", "tidyr")
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
  tables = file.path(project_root, "outputs", "tables"),
  supplementary = file.path(project_root, "outputs", "supplementary_tables"),
  derived_data = file.path(project_root, "outputs", "derived_data"),
  figure_source_data = file.path(project_root, "outputs", "figure_source_data")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# -----------------------------------------------------------------------------
# 2. Corrected module identities and order
# -----------------------------------------------------------------------------

module_key <- tribble(
  ~module_order, ~module_id,             ~disease_group,                        ~cohort, ~module_colour,
  1L,            "AD_ERC_greenyellow",  "Alzheimer's disease",                 "ERC",   "greenyellow",
  2L,            "AD_ERC_darkred",      "Alzheimer's disease",                 "ERC",   "darkred",
  3L,            "AD_HIPPO_yellow",     "Alzheimer's disease",                 "HIPPO", "yellow",
  4L,            "FTLD1_red",           "Frontotemporal lobar degeneration",   "FTLD1", "red",
  5L,            "FTLD2_darkmagenta",   "Frontotemporal lobar degeneration",   "FTLD2", "darkmagenta",
  6L,            "FTLD3_black",         "Frontotemporal lobar degeneration",   "FTLD3", "black",
  7L,            "FTLD3_turquoise",     "Frontotemporal lobar degeneration",   "FTLD3", "turquoise",
  8L,            "MSA_violet",          "Multiple system atrophy",             "MSA",   "violet",
  9L,            "PSP_pink",            "Progressive supranuclear palsy",       "PSP",   "pink"
)

# -----------------------------------------------------------------------------
# 3. Import candidate membership records from the corrected nine modules
# -----------------------------------------------------------------------------

# Script 06 creates this neutral, project-relative file from the published
# network objects and curated metal lists. Each row is an annotated CpG–gene
# record; this script collapses it to one gene–module occurrence before counting.
candidate_input_candidates <- c(
  file.path(
    project_root, "outputs", "derived_data",
    "nine_module_MM_GS_annotated_records.csv"
  ),
  file.path(
    project_root, "outputs", "derived_data",
    "nine_module_metal_candidate_records.csv"
  )
)
candidate_input <- candidate_input_candidates[file.exists(candidate_input_candidates)][1]

if (is.na(candidate_input)) {
  stop(
    "Corrected nine-module candidate records were not found. Run ",
    "R/06_candidate_prioritisation_MM_GS.R first."
  )
}

candidate_records <- read_csv(candidate_input, show_col_types = FALSE)

required_columns <- c(
  "module_id", "gene_symbol", "copper_gene", "broad_iron_gene", "core_iron_gene"
)
missing_columns <- setdiff(required_columns, names(candidate_records))
if (length(missing_columns) > 0L) {
  stop(
    "The candidate input is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  str_to_lower(str_trim(as.character(x))) %in% c("true", "t", "1", "yes")
}

candidate_records <- candidate_records |>
  transmute(
    module_id = as.character(module_id),
    gene_symbol = str_to_upper(str_trim(as.character(gene_symbol))),
    copper = as_logical_flag(copper_gene),
    broad_iron = as_logical_flag(broad_iron_gene),
    core_iron = as_logical_flag(core_iron_gene)
  ) |>
  filter(
    module_id %in% module_key$module_id,
    !is.na(gene_symbol), gene_symbol != "",
    copper | broad_iron | core_iron
  )

if (!all(candidate_records$core_iron <= candidate_records$broad_iron)) {
  stop("Every core-iron candidate must also be present in the broad-iron set.")
}

missing_modules <- setdiff(module_key$module_id, unique(candidate_records$module_id))
if (length(missing_modules) > 0L) {
  stop(
    "No metal-associated candidates were found for: ",
    paste(missing_modules, collapse = ", ")
  )
}

# One occurrence means that a gene is present in one retained module–disease
# comparison. Multiple CpGs mapping to the same gene within that module count once.
candidate_occurrence_long <- candidate_records |>
  group_by(module_id, gene_symbol) |>
  summarise(
    copper = any(copper),
    broad_iron = any(broad_iron),
    core_iron = any(core_iron),
    .groups = "drop"
  ) |>
  left_join(module_key, by = "module_id") |>
  arrange(module_order, gene_symbol)

# -----------------------------------------------------------------------------
# 4. Count recurrence and retain metal-set membership
# -----------------------------------------------------------------------------

candidate_recurrence <- candidate_occurrence_long |>
  group_by(gene_symbol) |>
  arrange(module_order, .by_group = TRUE) |>
  summarise(
    occurrences = n_distinct(module_id),
    retained_modules = paste(unique(module_id), collapse = "; "),
    copper = any(copper),
    broad_iron = any(broad_iron),
    core_iron = any(core_iron),
    .groups = "drop"
  ) |>
  mutate(
    status = if_else(occurrences >= 2L, "Recurrent", "Single module"),
    recurrence_label = paste0(occurrences, "/9")
  ) |>
  arrange(desc(occurrences), gene_symbol)

recurrent_candidates <- candidate_recurrence |>
  filter(occurrences >= 2L)

recurrent_copper <- recurrent_candidates |>
  filter(copper) |>
  arrange(desc(occurrences), gene_symbol)

recurrent_broad_iron <- recurrent_candidates |>
  filter(broad_iron) |>
  arrange(desc(occurrences), gene_symbol)

recurrent_core_iron <- recurrent_candidates |>
  filter(core_iron) |>
  arrange(desc(occurrences), gene_symbol)

single_module_candidates <- candidate_recurrence |>
  filter(occurrences == 1L) |>
  arrange(gene_symbol)

# -----------------------------------------------------------------------------
# 5. Validate the established corrected results
# -----------------------------------------------------------------------------

observed_totals <- c(
  retained_modules = n_distinct(candidate_occurrence_long$module_id),
  gene_module_occurrences = nrow(candidate_occurrence_long),
  unique_candidates = nrow(candidate_recurrence),
  recurrent_candidates = nrow(recurrent_candidates),
  single_module_candidates = nrow(single_module_candidates),
  recurrent_copper = nrow(recurrent_copper),
  recurrent_broad_iron = nrow(recurrent_broad_iron),
  recurrent_core_iron = nrow(recurrent_core_iron)
)

expected_totals <- c(
  retained_modules = 9L,
  gene_module_occurrences = 310L,
  unique_candidates = 140L,
  recurrent_candidates = 89L,
  single_module_candidates = 51L,
  recurrent_copper = 13L,
  recurrent_broad_iron = 77L,
  recurrent_core_iron = 57L
)

if (!identical(unname(observed_totals), unname(expected_totals))) {
  comparison <- tibble(
    result = names(expected_totals),
    expected = as.integer(expected_totals),
    observed = as.integer(observed_totals)
  )
  print(comparison)
  stop(
    "The recurrence counts do not match the established corrected results. ",
    "Check that Scripts 02, 04 and 06 used the final nine modules and gene sets."
  )
}

recurrence_summary <- tibble(
  result = c(
    "Retained metal-enriched modules",
    "Gene-module occurrences",
    "Unique metal-associated candidates",
    "Recurrent candidates (at least two modules)",
    "Single-module candidates",
    "Recurrent copper candidates",
    "Recurrent broad-iron candidates",
    "Recurrent core-iron candidates"
  ),
  n = as.integer(observed_totals)
)

# -----------------------------------------------------------------------------
# 6. Publication-ready Supplementary Table S5
# -----------------------------------------------------------------------------

describe_membership <- function(copper, broad_iron, core_iron) {
  case_when(
    copper & core_iron ~ "Copper; core iron (also broad iron)",
    copper & broad_iron ~ "Copper; broad iron",
    copper ~ "Copper",
    core_iron ~ "Core iron (also broad iron)",
    broad_iron ~ "Broad iron",
    TRUE ~ NA_character_
  )
}

supplementary_table_s5 <- candidate_recurrence |>
  mutate(
    `Metal-set membership` = describe_membership(copper, broad_iron, core_iron),
    No. = row_number()
  ) |>
  transmute(
    No.,
    `Gene symbol` = gene_symbol,
    Occurrences = occurrences,
    Status = status,
    `Metal-set membership`,
    `Retained module(s)` = retained_modules
  )

# -----------------------------------------------------------------------------
# 7. Occurrence matrices and downstream figure inputs
# -----------------------------------------------------------------------------

make_occurrence_long <- function(candidate_table) {
  candidate_occurrence_long |>
    semi_join(candidate_table, by = "gene_symbol") |>
    left_join(
      candidate_table |>
        select(gene_symbol, occurrences, recurrence_label),
      by = "gene_symbol"
    ) |>
    arrange(desc(occurrences), gene_symbol, module_order)
}

make_occurrence_matrix <- function(candidate_table) {
  gene_order <- candidate_table |>
    arrange(desc(occurrences), gene_symbol) |>
    pull(gene_symbol)

  candidate_occurrence_long |>
    semi_join(candidate_table, by = "gene_symbol") |>
    transmute(gene_symbol, module_id, present = 1L) |>
    complete(
      gene_symbol = gene_order,
      module_id = module_key$module_id,
      fill = list(present = 0L)
    ) |>
    mutate(
      gene_symbol = factor(gene_symbol, levels = gene_order),
      module_id = factor(module_id, levels = module_key$module_id)
    ) |>
    arrange(gene_symbol, module_id) |>
    pivot_wider(names_from = module_id, values_from = present) |>
    mutate(gene_symbol = as.character(gene_symbol))
}

top_20_broad_iron <- recurrent_broad_iron |>
  slice_head(n = 20L)

top_25_broad_iron <- recurrent_broad_iron |>
  slice_head(n = 25L)

copper_occurrence_long <- make_occurrence_long(recurrent_copper)
broad_iron_occurrence_long <- make_occurrence_long(recurrent_broad_iron)
core_iron_occurrence_long <- make_occurrence_long(recurrent_core_iron)

copper_occurrence_matrix <- make_occurrence_matrix(recurrent_copper)
broad_iron_occurrence_matrix <- make_occurrence_matrix(recurrent_broad_iron)
core_iron_occurrence_matrix <- make_occurrence_matrix(recurrent_core_iron)

# Long query table used by the recurrent-set functional-enrichment script.
recurrent_gene_queries <- bind_rows(
  transmute(recurrent_copper, query = "Recurrent_copper", gene_symbol),
  transmute(recurrent_broad_iron, query = "Recurrent_broad_iron", gene_symbol),
  transmute(recurrent_core_iron, query = "Recurrent_core_iron", gene_symbol)
) |>
  arrange(query, gene_symbol)

# -----------------------------------------------------------------------------
# 8. Save all outputs
# -----------------------------------------------------------------------------

write_csv(
  candidate_occurrence_long,
  file.path(output_dirs$derived_data, "nine_module_candidate_occurrence_long.csv")
)
write_csv(
  candidate_recurrence,
  file.path(output_dirs$tables, "nine_module_all_candidate_recurrence.csv")
)
write_csv(
  recurrent_candidates,
  file.path(output_dirs$tables, "nine_module_all_recurrent_candidates.csv")
)
write_csv(
  recurrent_copper,
  file.path(output_dirs$tables, "nine_module_recurrent_copper_candidates.csv")
)
write_csv(
  recurrent_broad_iron,
  file.path(output_dirs$tables, "nine_module_recurrent_broad_iron_candidates.csv")
)
write_csv(
  recurrent_core_iron,
  file.path(output_dirs$tables, "nine_module_recurrent_core_iron_candidates.csv")
)
write_csv(
  single_module_candidates,
  file.path(output_dirs$tables, "nine_module_single_module_candidates.csv")
)
write_csv(
  recurrence_summary,
  file.path(output_dirs$tables, "nine_module_recurrence_summary.csv")
)
write_csv(
  supplementary_table_s5,
  file.path(
    output_dirs$supplementary,
    "Supplementary_Table_S5_candidate_recurrence_and_occurrence.csv"
  )
)

write_csv(
  copper_occurrence_long,
  file.path(output_dirs$figure_source_data, "Figures_7_9_recurrent_copper_source_data.csv")
)
write_csv(
  broad_iron_occurrence_long,
  file.path(output_dirs$figure_source_data, "Figures_8_10_recurrent_broad_iron_source_data.csv")
)
write_csv(
  core_iron_occurrence_long,
  file.path(output_dirs$figure_source_data, "recurrent_core_iron_source_data.csv")
)
write_csv(
  copper_occurrence_matrix,
  file.path(output_dirs$derived_data, "recurrent_copper_occurrence_matrix.csv")
)
write_csv(
  broad_iron_occurrence_matrix,
  file.path(output_dirs$derived_data, "recurrent_broad_iron_occurrence_matrix.csv")
)
write_csv(
  core_iron_occurrence_matrix,
  file.path(output_dirs$derived_data, "recurrent_core_iron_occurrence_matrix.csv")
)
write_csv(
  top_20_broad_iron,
  file.path(output_dirs$derived_data, "top20_recurrent_broad_iron_candidates.csv")
)
write_csv(
  top_25_broad_iron,
  file.path(output_dirs$derived_data, "top25_recurrent_broad_iron_candidates.csv")
)
write_csv(
  recurrent_gene_queries,
  file.path(output_dirs$derived_data, "recurrent_gene_queries.csv")
)

message("R/07_recurrence_analysis_9_modules.R completed successfully.")
message("Verified 140 unique and 89 recurrent metal-associated candidates.")
message("Verified 13 recurrent copper, 77 broad-iron and 57 core-iron genes.")
message("Supplementary Table S5: ", output_dirs$supplementary)
