
# Final automated validation of the corrected nine-module analysis
#
# Checks numerical results reported in the thesis, confirms that required figures
# and tables exist, scans final text-based outputs for obsolete 16-module or
# six-FTLD-module wording, and records package/session information.

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
# 1. Project-relative paths and audit collectors
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
outputs_root <- file.path(project_root, "outputs")
validation_dir <- file.path(outputs_root, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()

record_check <- function(
    category, check_name, expected, observed, passed,
    severity = "ERROR", note = "") {
  checks[[length(checks) + 1L]] <<- tibble(
    category = category,
    check = check_name,
    expected = as.character(expected),
    observed = as.character(observed),
    passed = isTRUE(passed),
    severity = severity,
    note = note
  )
}

read_csv_if_present <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}

count_if_present <- function(data) {
  if (is.null(data)) NA_integer_ else nrow(data)
}

# -----------------------------------------------------------------------------
# 2. Confirm all clean analysis scripts exist
# -----------------------------------------------------------------------------

required_scripts <- sprintf("R/%02d_%s", 0:12, c(
  "setup.R",
  "preprocess_AD_methylation.R",
  "curate_metal_gene_sets.R",
  "select_disease_associated_modules.R",
  "metal_gene_set_enrichment.R",
  "integrate_published_EWCE.R",
  "candidate_prioritisation_MM_GS.R",
  "recurrence_analysis_9_modules.R",
  "plot_occurrence_matrices.R",
  "plot_circos_networks.R",
  "functional_enrichment_gprofiler.R",
  "integrate_FTLDexp.R",
  "validate_final_outputs.R"
))

required_scripts <- c("R/utils.R", required_scripts)

for (script in required_scripts) {
  present <- file.exists(file.path(project_root, script))
  record_check(
    "Scripts", script, "Present", ifelse(present, "Present", "Missing"), present
  )
}

repository_files <- c(".gitignore", "README.md", "renv.lock", "run_all.R")
for (repository_file in repository_files) {
  present <- file.exists(file.path(project_root, repository_file))
  record_check(
    "Repository",
    repository_file,
    "Present",
    ifelse(present, "Present", "Missing"),
    present
  )
}

# -----------------------------------------------------------------------------
# 3. Validate curated gene sets
# -----------------------------------------------------------------------------

gene_set_candidates <- c(
  file.path(outputs_root, "derived_data", "curated_metal_gene_sets.rds"),
  file.path(project_root, "data", "derived", "curated_metal_gene_sets.rds"),
  file.path(project_root, "gene_lists", "curated_metal_gene_sets.rds")
)
gene_set_file <- gene_set_candidates[file.exists(gene_set_candidates)][1]

if (is.na(gene_set_file)) {
  record_check(
    "Gene sets", "Curated metal-gene-set object", "Present", "Missing", FALSE
  )
} else {
  curated_sets <- readRDS(gene_set_file)
  expected_gene_set_counts <- c(copper = 47L, broad_iron = 199L, core_iron = 136L)
  for (set_name in names(expected_gene_set_counts)) {
    observed <- if (set_name %in% names(curated_sets)) {
      length(unique(curated_sets[[set_name]]))
    } else {
      NA_integer_
    }
    record_check(
      "Gene sets",
      paste0(set_name, " genes"),
      expected_gene_set_counts[[set_name]],
      observed,
      identical(as.integer(observed), expected_gene_set_counts[[set_name]])
    )
  }

  subset_valid <- all(
    unique(curated_sets$core_iron) %in% unique(curated_sets$broad_iron)
  )
  record_check(
    "Gene sets",
    "Core iron is a subset of broad iron",
    TRUE,
    subset_valid,
    subset_valid
  )
}

# -----------------------------------------------------------------------------
# 4. Validate corrected enrichment and recurrence totals
# -----------------------------------------------------------------------------

recurrence_summary_file <- file.path(
  outputs_root, "tables", "nine_module_recurrence_summary.csv"
)
recurrence_summary <- read_csv_if_present(recurrence_summary_file)

expected_recurrence <- c(
  "Retained metal-enriched modules" = 9L,
  "Gene-module occurrences" = 310L,
  "Unique metal-associated candidates" = 140L,
  "Recurrent candidates (at least two modules)" = 89L,
  "Single-module candidates" = 51L,
  "Recurrent copper candidates" = 13L,
  "Recurrent broad-iron candidates" = 77L,
  "Recurrent core-iron candidates" = 57L
)

if (is.null(recurrence_summary) ||
    !all(c("result", "n") %in% names(recurrence_summary))) {
  record_check(
    "Recurrence", "Recurrence summary", "Present and readable", "Missing", FALSE
  )
} else {
  observed_recurrence <- setNames(recurrence_summary$n, recurrence_summary$result)
  for (result_name in names(expected_recurrence)) {
    observed <- unname(observed_recurrence[result_name])
    record_check(
      "Recurrence",
      result_name,
      expected_recurrence[[result_name]],
      observed,
      length(observed) == 1L && !is.na(observed) &&
        as.integer(observed) == expected_recurrence[[result_name]]
    )
  }
}

# Validate the 12 significant module–gene-set results when the Script 04 output
# is available. Several accepted filenames are supported to keep the validator
# compatible with a clean repository layout.
enrichment_result_candidates <- c(
  file.path(outputs_root, "derived_data", "significant_metal_gene_set_enrichment.csv"),
  file.path(outputs_root, "tables", "complete_significant_metal_enrichment.csv"),
  file.path(outputs_root, "supplementary_tables", "Supplementary_Table_S4_complete_module_enrichment_results.csv")
)
enrichment_file <- enrichment_result_candidates[file.exists(enrichment_result_candidates)][1]

if (is.na(enrichment_file) && dir.exists(outputs_root)) {
  s4_candidates <- list.files(
    outputs_root,
    pattern = "Supplementary[_ ]Table[_ ]S4.*[.]csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(s4_candidates) > 0L) enrichment_file <- s4_candidates[1]
}

if (is.na(enrichment_file)) {
  record_check(
    "Metal enrichment",
    "Complete significant module–gene-set result file",
    "Present",
    "Missing",
    FALSE,
    note = "Expected output from Script 04/Supplementary Table S4"
  )
} else {
  enrichment_results <- read_csv(enrichment_file, show_col_types = FALSE)
  module_column <- c("module_id", "module", "Module")[
    c("module_id", "module", "Module") %in% names(enrichment_results)
  ][1]
  observed_rows <- nrow(enrichment_results)
  observed_modules <- if (!is.na(module_column)) {
    n_distinct(enrichment_results[[module_column]])
  } else {
    NA_integer_
  }
  record_check(
    "Metal enrichment", "Significant module–gene-set results",
    12L, observed_rows, observed_rows == 12L
  )
  record_check(
    "Metal enrichment", "Metal-enriched modules",
    9L, observed_modules, observed_modules == 9L
  )
}

# -----------------------------------------------------------------------------
# 5. Validate MM–GS correlations
# -----------------------------------------------------------------------------

correlation_file <- file.path(
  outputs_root, "tables", "MM_GS_spearman_correlations_nine_modules.csv"
)
correlations <- read_csv_if_present(correlation_file)

expected_rho <- c(
  AD_ERC_greenyellow = 0.66,
  FTLD1_red = 0.49,
  FTLD2_darkmagenta = 0.45,
  FTLD3_turquoise = 0.47
)

if (is.null(correlations) ||
    !all(c("module_id", "spearman_rho") %in% names(correlations))) {
  record_check(
    "MM-GS", "Spearman correlation table", "Present and readable", "Missing", FALSE
  )
} else {
  observed_rho <- setNames(correlations$spearman_rho, correlations$module_id)
  for (module_id in names(expected_rho)) {
    observed <- round(unname(observed_rho[module_id]), 2)
    record_check(
      "MM-GS", paste0(module_id, " Spearman rho"),
      expected_rho[[module_id]], observed,
      length(observed) == 1L && !is.na(observed) && observed == expected_rho[[module_id]]
    )
  }
}

# -----------------------------------------------------------------------------
# 6. Validate functional-enrichment totals
# -----------------------------------------------------------------------------

module_function_file <- file.path(
  outputs_root, "derived_data", "nine_module_gprofiler_significant_results.csv"
)
recurrent_function_file <- file.path(
  outputs_root, "derived_data", "recurrent_sets_gprofiler_significant_results.csv"
)
module_function <- read_csv_if_present(module_function_file)
recurrent_function <- read_csv_if_present(recurrent_function_file)

if (is.null(module_function)) {
  record_check(
    "Functional enrichment", "Module-level significant terms",
    771L, "Missing", FALSE
  )
} else {
  record_check(
    "Functional enrichment", "Module-level significant terms",
    771L, nrow(module_function), nrow(module_function) == 771L,
    note = "Exact submitted result; a later live g:Profiler release may differ"
  )
  module_name_column <- c("query_name", "module", "module_id")[
    c("query_name", "module", "module_id") %in% names(module_function)
  ][1]
  modules_with_terms <- if (!is.na(module_name_column)) {
    n_distinct(module_function[[module_name_column]])
  } else {
    NA_integer_
  }
  record_check(
    "Functional enrichment", "Modules with significant terms",
    8L, modules_with_terms, modules_with_terms == 8L
  )
}

if (is.null(recurrent_function)) {
  record_check(
    "Functional enrichment", "Recurrent-set significant terms",
    0L, "Missing", FALSE
  )
} else {
  record_check(
    "Functional enrichment", "Recurrent-set significant terms",
    0L, nrow(recurrent_function), nrow(recurrent_function) == 0L
  )
}

# -----------------------------------------------------------------------------
# 7. Validate FTLDexp integration
# -----------------------------------------------------------------------------

ftld_integration_file <- file.path(
  outputs_root,
  "supplementary_tables",
  "Supplementary_Table_S7_complete_FTLDexp_candidate_integration.csv"
)
ftld_integration <- read_csv_if_present(ftld_integration_file)

if (is.null(ftld_integration)) {
  record_check(
    "FTLDexp", "Complete candidate integration", "Present", "Missing", FALSE
  )
} else {
  required_ftld_columns <- c(
    "Module", "Gene symbol", "Matched in FTLDexp",
    "Significantly differentially expressed"
  )
  if (!all(required_ftld_columns %in% names(ftld_integration))) {
    record_check(
      "FTLDexp", "Complete candidate integration columns",
      paste(required_ftld_columns, collapse = "; "),
      paste(names(ftld_integration), collapse = "; "),
      FALSE
    )
  } else {
    significant_flag <- as.logical(
      ftld_integration$`Significantly differentially expressed`
    )
    matched_flag <- as.logical(ftld_integration$`Matched in FTLDexp`)
    significant_table <- ftld_integration[significant_flag %in% TRUE, , drop = FALSE]

    ftld_checks <- list(
      "Candidate module-gene occurrences evaluated" = c(188L, nrow(ftld_integration)),
      "Unique candidates evaluated" = c(117L, n_distinct(ftld_integration$`Gene symbol`)),
      "Matched module-gene occurrences" = c(185L, sum(matched_flag, na.rm = TRUE)),
      "Matched unique candidates" = c(
        114L,
        n_distinct(ftld_integration$`Gene symbol`[matched_flag %in% TRUE])
      ),
      "Significant module-gene occurrences" = c(40L, nrow(significant_table)),
      "Unique differentially expressed genes" = c(
        25L,
        n_distinct(significant_table$`Gene symbol`)
      )
    )

    for (check_name in names(ftld_checks)) {
      values <- ftld_checks[[check_name]]
      record_check(
        "FTLDexp", check_name, values[1], values[2], values[1] == values[2]
      )
    }

    expected_ftld_module_counts <- c(
      "FTLD1 red" = 3L,
      "FTLD2 darkmagenta" = 4L,
      "FTLD3 black" = 18L,
      "FTLD3 turquoise" = 15L
    )
    observed_ftld_module_counts <- table(significant_table$Module)
    for (module in names(expected_ftld_module_counts)) {
      observed <- unname(observed_ftld_module_counts[module])
      if (length(observed) == 0L || is.na(observed)) observed <- 0L
      record_check(
        "FTLDexp", paste0(module, " significant occurrences"),
        expected_ftld_module_counts[[module]], observed,
        observed == expected_ftld_module_counts[[module]]
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 8. Confirm required figures and tables exist
# -----------------------------------------------------------------------------

all_output_files <- if (dir.exists(outputs_root)) {
  list.files(outputs_root, recursive = TRUE, full.names = TRUE)
} else {
  character()
}
relative_output_files <- str_remove(
  all_output_files,
  paste0("^", stringr::fixed(normalizePath(outputs_root, mustWork = FALSE)), "/?")
)

check_pattern_exists <- function(category, label, pattern) {
  matches <- relative_output_files[str_detect(relative_output_files, regex(pattern, ignore_case = TRUE))]
  record_check(
    category,
    label,
    "At least one matching file",
    ifelse(length(matches) > 0L, paste(matches, collapse = "; "), "Missing"),
    length(matches) > 0L
  )
}

for (figure_number in 5:12) {
  check_pattern_exists(
    "Deliverables",
    paste0("Figure ", figure_number),
    paste0("(^|/)Figure[_ ]?", figure_number, ".*[.](png|tiff|pdf)$")
  )
}
check_pattern_exists(
  "Deliverables",
  "Supplementary Figure S1",
  "(^|/)Supplementary[_ ]Figure[_ ]S1.*[.](png|tiff|pdf)$"
)
check_pattern_exists(
  "Deliverables",
  "Supplementary Figure S2",
  "(^|/)Supplementary[_ ]Figure[_ ]S2.*[.](png|tiff|pdf)$"
)

for (table_label in c("3A", "3B", "4", "5", "6")) {
  check_pattern_exists(
    "Deliverables",
    paste0("Table ", table_label),
    paste0("(^|/)Table[_ ]?", table_label, ".*[.](csv|docx)$")
  )
}

for (supplementary_number in 1:7) {
  check_pattern_exists(
    "Deliverables",
    paste0("Supplementary Table S", supplementary_number),
    paste0(
      "(^|/)Supplementary[_ ]Table[_ ]S", supplementary_number,
      ".*[.](csv|docx|xlsx)$"
    )
  )
}

# -----------------------------------------------------------------------------
# 9. Scan final outputs for obsolete wording or removed FTLD modules
# -----------------------------------------------------------------------------

obsolete_patterns <- c(
  "(?i)\\b16\\s+(metal[- ]enriched\\s+)?modules\\b",
  "(?i)\\b16[- ]module\\b",
  "(?i)\\bacross\\s+16\\s+modules\\b",
  "(?i)/16\\b",
  "(?i)\\bsix\\s+(metal[- ]enriched\\s+)?FTLD\\s+modules\\b",
  "(?i)\\ball\\s+six\\s+FTLD\\s+modules\\b",
  "(?i)FTLD1[ _-]+darkgrey",
  "(?i)FTLD2[ _-]+skyblue"
)

text_extensions <- "[.](csv|tsv|txt|md|json|ya?ml|html|xml)$"
text_files <- all_output_files[str_detect(all_output_files, regex(text_extensions))]
text_files <- text_files[!str_starts(text_files, validation_dir)]

scan_hits <- list()
for (file in text_files) {
  lines <- tryCatch(
    readLines(file, warn = FALSE, encoding = "UTF-8"),
    error = function(error) character()
  )
  if (length(lines) == 0L) next

  for (pattern in obsolete_patterns) {
    matching_lines <- which(str_detect(lines, regex(pattern)))
    if (length(matching_lines) > 0L) {
      scan_hits[[length(scan_hits) + 1L]] <- tibble(
        file = str_remove(file, paste0("^", fixed(project_root), "/?")),
        line = matching_lines,
        pattern = pattern,
        text = str_squish(lines[matching_lines])
      )
    }
  }
}

# Scan editable Word tables by reading their XML contents.
docx_files <- all_output_files[str_detect(all_output_files, regex("[.]docx$", ignore_case = TRUE))]
docx_files <- docx_files[!str_starts(docx_files, validation_dir)]

read_docx_text <- function(file) {
  connection <- unz(file, "word/document.xml", open = "r")
  on.exit(close(connection), add = TRUE)
  xml_lines <- readLines(connection, warn = FALSE, encoding = "UTF-8")
  str_replace_all(paste(xml_lines, collapse = " "), "<[^>]+>", " ") |>
    str_squish()
}

for (file in docx_files) {
  xml_text <- tryCatch(
    read_docx_text(file),
    error = function(error) character()
  )

  if (length(xml_text) == 0L) next

  for (pattern in obsolete_patterns) {
    if (str_detect(xml_text, regex(pattern))) {
      scan_hits[[length(scan_hits) + 1L]] <- tibble(
        file = str_remove(file, paste0("^", fixed(project_root), "/?")),
        line = NA_integer_,
        pattern = pattern,
        text = str_sub(xml_text, 1, 500)
      )
    }
  }
}

obsolete_wording <- if (length(scan_hits) == 0L) {
  tibble(file = character(), line = integer(), pattern = character(), text = character())
} else {
  bind_rows(scan_hits) |>
    distinct()
}

write_csv(
  obsolete_wording,
  file.path(validation_dir, "obsolete_wording_audit.csv")
)
record_check(
  "Terminology",
  "Obsolete 16-module/six-FTLD wording and removed modules",
  "Zero matches",
  nrow(obsolete_wording),
  nrow(obsolete_wording) == 0L,
  note = "Scanned final text outputs and editable Word-table XML"
)

# -----------------------------------------------------------------------------
# 10. Record session, package versions and output manifest
# -----------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(validation_dir, "sessionInfo.txt")
)

packages_to_record <- unique(c(
  "ChAMP", "minfi", "WGCNA", "dplyr", "tidyr", "readr", "purrr",
  "ggplot2", "ggrepel", "patchwork", "circlize", "RColorBrewer",
  "gprofiler2", "limma", "flextable", "officer"
))

package_versions <- tibble(
  package = packages_to_record,
  installed = vapply(
    packages_to_record,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
) |>
  mutate(
    version = vapply(
      package,
      function(package_name) {
        if (!requireNamespace(package_name, quietly = TRUE)) return(NA_character_)
        as.character(utils::packageVersion(package_name))
      },
      character(1)
    )
  )

write_csv(
  package_versions,
  file.path(validation_dir, "package_versions.csv")
)

manifest_files <- unique(c(
  list.files(file.path(project_root, "R"), full.names = TRUE, recursive = TRUE),
  file.path(project_root, repository_files),
  all_output_files
))
manifest_files <- manifest_files[file.exists(manifest_files)]
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]

output_manifest <- tibble(
  file = str_remove(
    manifest_files,
    paste0("^", fixed(project_root), "/?")
  ),
  size_bytes = file.info(manifest_files)$size,
  modified = as.character(file.info(manifest_files)$mtime),
  md5 = unname(tools::md5sum(manifest_files))
) |>
  arrange(file)

write_csv(
  output_manifest,
  file.path(validation_dir, "final_output_manifest.csv")
)

# -----------------------------------------------------------------------------
# 11. Write final report and stop on required failures
# -----------------------------------------------------------------------------

validation_report <- bind_rows(checks) |>
  arrange(desc(severity), category, check)

write_csv(
  validation_report,
  file.path(validation_dir, "final_validation_report.csv")
)

validation_summary <- validation_report |>
  count(severity, passed, name = "checks")
write_csv(
  validation_summary,
  file.path(validation_dir, "final_validation_summary.csv")
)

failed_required <- validation_report |>
  filter(!passed, severity == "ERROR")

message("Final validation report: ", validation_dir)
message(
  "Checks passed: ", sum(validation_report$passed),
  "/", nrow(validation_report)
)

if (nrow(failed_required) > 0L) {
  print(failed_required |>
    select(category, check, expected, observed, note), n = Inf)
  stop(
    nrow(failed_required),
    " required validation check(s) failed. Review final_validation_report.csv."
  )
}

message("All required final-output validation checks passed.")
