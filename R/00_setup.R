# ==============================================================================
# 00_setup.R
# MSc thesis: Exploring the contributions of brain DNA methylation changes to copper and iron metabolism dysregulation across neurodegenerative diseases

#
# Purpose:
#   Define reproducible project settings, check essential packages, establish
#   project-relative paths and create the output directories used throughout
#   the analysis.
# ==============================================================================

# Reproducible settings ---------------------------------------------------------

set.seed(12345)

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  dplyr.summarise.inform = FALSE,
  warn = 1
)

# Essential packages -----------------------------------------------------------

essential_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "ggplot2"
)

missing_packages <- essential_packages[
  !vapply(
    essential_packages,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "The following required packages are not installed: ",
      paste(missing_packages, collapse = ", "),
      ".\nInstall them before continuing, preferably using renv."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

# Project root -----------------------------------------------------------------

# Locate the repository root without a user-specific path or optional helper
# package. This works from the project root, R/, or through run_all.R.
locate_project_root <- function(start = getwd()) {
  candidate <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(candidate, "run_all.R")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop(
    "Could not locate the repository root containing run_all.R and R/.",
    call. = FALSE
  )
}

project_root <- locate_project_root()

if (!file.exists(file.path(
  project_root,
  "ChAMP_CuFe_Methylation_Project.Rproj"
))) {
  warning(
    paste0(
      "The expected RStudio project file was not found in:\n",
      project_root,
      "\nOpen the project from its root directory before running the pipeline."
    ),
    call. = FALSE
  )
}

# Project directories ----------------------------------------------------------

paths <- list(
  root = project_root,

  # Scripts
  scripts = file.path(project_root, "R"),

  # Input data
  data = file.path(project_root, "data"),
  data_raw = file.path(project_root, "data", "raw"),
  data_processed = file.path(project_root, "data", "processed"),
  metadata = file.path(project_root, "data", "metadata"),
  gene_sets = file.path(project_root, "data", "gene_sets"),
  published_results = file.path(project_root, "data", "published_results"),
  ftldexp = file.path(project_root, "data", "FTLDexp"),

  # Analysis outputs
  outputs = file.path(project_root, "outputs"),
  derived_data = file.path(project_root, "outputs", "derived_data"),
  results = file.path(project_root, "outputs", "results"),
  tables = file.path(project_root, "outputs", "tables"),
  supplementary_tables = file.path(
    project_root,
    "outputs",
    "supplementary_tables"
  ),
  figures = file.path(project_root, "outputs", "figures"),
  supplementary_figures = file.path(
    project_root,
    "outputs",
    "supplementary_figures"
  ),
  figure_data = file.path(project_root, "outputs", "figure_source_data"),
  qc = file.path(project_root, "outputs", "quality_control"),
  logs = file.path(project_root, "outputs", "logs")
)

directories_to_create <- unname(unlist(paths[c(
  "scripts",
  "data",
  "data_raw",
  "data_processed",
  "metadata",
  "gene_sets",
  "published_results",
  "ftldexp",
  "outputs",
  "derived_data",
  "results",
  "tables",
  "supplementary_tables",
  "figures",
  "supplementary_figures",
  "figure_data",
  "qc",
  "logs"
)]))

invisible(
  lapply(
    directories_to_create,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# Package-checking helper -------------------------------------------------------

check_packages <- function(packages, analysis_name = "this analysis") {
  unavailable <- packages[
    !vapply(
      packages,
      requireNamespace,
      FUN.VALUE = logical(1),
      quietly = TRUE
    )
  ]

  if (length(unavailable) > 0L) {
    stop(
      paste0(
        "Packages required for ",
        analysis_name,
        " are not installed: ",
        paste(unavailable, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# File-checking helper ----------------------------------------------------------

check_input_files <- function(files, description = "input files") {
  missing_files <- files[!file.exists(files)]

  if (length(missing_files) > 0L) {
    stop(
      paste0(
        "The following ",
        description,
        " could not be found:\n",
        paste(paste0(" - ", missing_files), collapse = "\n")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Gene-symbol standardisation helper -------------------------------------------

standardise_gene_symbols <- function(symbols) {
  symbols |>
    as.character() |>
    stringr::str_trim() |>
    stringr::str_to_upper() |>
    dplyr::na_if("")
}

# Safe table export -------------------------------------------------------------

write_output_csv <- function(data, filename, output_directory = paths$results) {
  if (!is.data.frame(data)) {
    stop("The object supplied to write_output_csv() is not a data frame.")
  }

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_file <- file.path(output_directory, filename)

  readr::write_csv(
    data,
    output_file,
    na = ""
  )

  message("Saved: ", output_file)

  invisible(output_file)
}

# Consistent figure export ------------------------------------------------------

save_thesis_figure <- function(
    plot,
    filename,
    width,
    height,
    units = "in",
    dpi = 600,
    output_directory = paths$figures,
    save_pdf = TRUE) {

  if (!inherits(plot, c("gg", "ggplot"))) {
    stop("The object supplied to save_thesis_figure() is not a ggplot.")
  }

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  png_file <- file.path(
    output_directory,
    paste0(filename, ".png")
  )

  ggplot2::ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    bg = "white"
  )

  message("Saved: ", png_file)

  if (isTRUE(save_pdf)) {
    pdf_file <- file.path(
      output_directory,
      paste0(filename, ".pdf")
    )

    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = units,
      device = grDevices::pdf,
      useDingbats = FALSE,
      bg = "white"
    )

    message("Saved: ", pdf_file)
  }

  invisible(png_file)
}

# Analysis constants ------------------------------------------------------------

analysis_constants <- list(
  random_seed = 12345L,
  nominal_p_threshold = 0.05,
  fdr_threshold = 0.05,
  recurrence_threshold = 2L,

  expected_gene_set_sizes = c(
    copper = 47L,
    broad_iron = 199L,
    core_iron = 136L
  ),

  expected_final_results = c(
    metal_enriched_modules = 9L,
    module_gene_set_results = 12L,
    unique_candidates = 140L,
    recurrent_candidates = 89L,
    recurrent_copper = 13L,
    recurrent_broad_iron = 77L,
    recurrent_core_iron = 57L,
    ftldexp_module_gene_occurrences = 40L,
    ftldexp_unique_genes = 25L
  )
)

# Record initial setup information ---------------------------------------------

setup_information <- tibble::tibble(
  item = c(
    "Analysis date",
    "R version",
    "Project root",
    "Random seed"
  ),
  value = c(
    as.character(Sys.Date()),
    R.version.string,
    project_root,
    as.character(analysis_constants$random_seed)
  )
)

readr::write_csv(
  setup_information,
  file.path(paths$logs, "analysis_setup.csv")
)

message("Project setup completed successfully.")
message("Project root: ", project_root)
