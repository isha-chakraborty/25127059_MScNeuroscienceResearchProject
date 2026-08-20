
# Run the complete corrected analysis in order.
#
# Default behaviour:
#   - AD preprocessing is skipped because it is resource-intensive and is not
#     required to reconstruct the final cross-disease analysis from the published
#     network outputs.
#   - g:Profiler results are reused when saved result objects exist; otherwise the
#     live queries are run.
#
# Optional command-line switches:
#   Rscript run_all.R --run-ad-preprocessing
#   Rscript run_all.R --skip-ad-preprocessing
#   Rscript run_all.R --run-gprofiler
#   Rscript run_all.R --reuse-gprofiler
#
# Equivalent environment variables:
#   RUN_AD_PREPROCESSING = "run" or "skip"
#   RUN_GPROFILER = "run", "reuse" or "auto"

options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------------------------------------------------------
# 1. Locate and enter the project root
# -----------------------------------------------------------------------------

find_project_root <- function() {
  command_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- command_arguments[grepl("^--file=", command_arguments)]

  if (length(file_argument) > 0L) {
    script_path <- sub("^--file=", "", file_argument[1])
    candidate <- normalizePath(dirname(script_path), mustWork = TRUE)
  } else {
    candidate <- normalizePath(getwd(), mustWork = TRUE)
  }

  repeat {
    if (file.exists(file.path(candidate, "run_all.R")) &&
        dir.exists(file.path(candidate, "R"))) {
      return(candidate)
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop("Could not locate the project root containing run_all.R and the R directory.")
}

project_root <- find_project_root()
original_working_directory <- getwd()
setwd(project_root)
on.exit(setwd(original_working_directory), add = TRUE)

log_directory <- file.path(project_root, "outputs", "logs")
dir.create(log_directory, recursive = TRUE, showWarnings = FALSE)

utils_file <- file.path(project_root, "R", "utils.R")
if (!file.exists(utils_file)) {
  stop("Shared utility file is missing: R/utils.R")
}
source(utils_file, local = .GlobalEnv)

# -----------------------------------------------------------------------------
# 2. Parse safe execution choices
# -----------------------------------------------------------------------------

normalise_choice <- function(value, allowed, default) {
  value <- tolower(trimws(value))
  if (value == "") value <- default
  if (!value %in% allowed) {
    stop(
      "Invalid choice '", value, "'. Allowed values: ",
      paste(allowed, collapse = ", ")
    )
  }
  value
}

ad_preprocessing_choice <- normalise_choice(
  Sys.getenv("RUN_AD_PREPROCESSING", unset = "skip"),
  allowed = c("run", "skip"),
  default = "skip"
)

gprofiler_choice <- normalise_choice(
  Sys.getenv("RUN_GPROFILER", unset = "auto"),
  allowed = c("run", "reuse", "auto"),
  default = "auto"
)

user_arguments <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% user_arguments) {
  cat(
    "Corrected Cu/Fe methylation analysis\n\n",
    "Usage:\n",
    "  Rscript run_all.R [options]\n\n",
    "Options:\n",
    "  --run-ad-preprocessing    Run resource-intensive AD preprocessing\n",
    "  --skip-ad-preprocessing   Use published/saved downstream inputs (default)\n",
    "  --run-gprofiler           Submit new live g:Profiler queries\n",
    "  --reuse-gprofiler         Use completed saved g:Profiler results\n",
    "  --help                    Show this message\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

known_arguments <- c(
  "--run-ad-preprocessing", "--skip-ad-preprocessing",
  "--run-gprofiler", "--reuse-gprofiler", "--help"
)
unknown_arguments <- setdiff(user_arguments, known_arguments)
if (length(unknown_arguments) > 0L) {
  stop("Unknown argument(s): ", paste(unknown_arguments, collapse = ", "))
}

if (all(c("--run-ad-preprocessing", "--skip-ad-preprocessing") %in% user_arguments)) {
  stop("Choose either --run-ad-preprocessing or --skip-ad-preprocessing, not both.")
}
if (all(c("--run-gprofiler", "--reuse-gprofiler") %in% user_arguments)) {
  stop("Choose either --run-gprofiler or --reuse-gprofiler, not both.")
}

if ("--run-ad-preprocessing" %in% user_arguments) {
  ad_preprocessing_choice <- "run"
}
if ("--skip-ad-preprocessing" %in% user_arguments) {
  ad_preprocessing_choice <- "skip"
}
if ("--run-gprofiler" %in% user_arguments) {
  gprofiler_choice <- "run"
}
if ("--reuse-gprofiler" %in% user_arguments) {
  gprofiler_choice <- "reuse"
}

saved_gprofiler_groups <- list(
  module = c(
    file.path(
      project_root, "outputs", "derived_data",
      "nine_module_gprofiler_significant_results.rds"
    ),
    file.path(
      project_root, "data", "derived",
      "submitted_nine_module_gprofiler_significant_results.rds"
    )
  ),
  recurrent = c(
    file.path(
      project_root, "outputs", "derived_data",
      "recurrent_sets_gprofiler_significant_results.rds"
    ),
    file.path(
      project_root, "data", "derived",
      "submitted_recurrent_sets_gprofiler_significant_results.rds"
    )
  )
)
saved_gprofiler_available <- all(vapply(
  saved_gprofiler_groups,
  function(candidates) any(file.exists(candidates)),
  logical(1)
))

if (gprofiler_choice == "auto") {
  gprofiler_choice <- if (saved_gprofiler_available) "reuse" else "run"
}
if (gprofiler_choice == "reuse" && !saved_gprofiler_available) {
  stop(
    "Saved g:Profiler result objects were requested but are missing. ",
    "Use --run-gprofiler once to create them."
  )
}

previous_reuse_setting <- Sys.getenv("REUSE_GPROFILER_RESULTS", unset = NA_character_)
if (gprofiler_choice == "reuse") {
  Sys.setenv(REUSE_GPROFILER_RESULTS = "1")
} else {
  Sys.setenv(REUSE_GPROFILER_RESULTS = "0")
}
on.exit({
  if (is.na(previous_reuse_setting)) {
    Sys.unsetenv("REUSE_GPROFILER_RESULTS")
  } else {
    Sys.setenv(REUSE_GPROFILER_RESULTS = previous_reuse_setting)
  }
}, add = TRUE)

# -----------------------------------------------------------------------------
# 3. Ordered workflow
# -----------------------------------------------------------------------------

workflow <- data.frame(
  step = 0:12,
  script = file.path("R", c(
    "00_setup.R",
    "01_preprocess_AD_methylation.R",
    "02_curate_metal_gene_sets.R",
    "03_select_disease_associated_modules.R",
    "04_metal_gene_set_enrichment.R",
    "05_integrate_published_EWCE.R",
    "06_candidate_prioritisation_MM_GS.R",
    "07_recurrence_analysis_9_modules.R",
    "08_plot_occurrence_matrices.R",
    "09_plot_circos_networks.R",
    "10_functional_enrichment_gprofiler.R",
    "11_integrate_FTLDexp.R",
    "12_validate_final_outputs.R"
  )),
  description = c(
    "Set up packages, paths and reproducibility settings",
    "Preprocess and normalise AD methylation data",
    "Curate copper, broad-iron and core-iron gene sets",
    "Select disease-associated modules",
    "Test metal gene-set enrichment",
    "Integrate published cell-type enrichment results",
    "Prioritise candidates using module membership and gene significance",
    "Calculate recurrence across nine modules",
    "Generate occurrence matrices",
    "Generate Circos plots",
    "Run or reuse functional-enrichment results",
    "Integrate FTLDexp differential-expression results",
    "Validate final outputs"
  ),
  run = TRUE,
  stringsAsFactors = FALSE
)

if (ad_preprocessing_choice == "skip") {
  workflow$run[workflow$step == 1L] <- FALSE
}

missing_scripts <- workflow$script[!file.exists(file.path(project_root, workflow$script))]
if (length(missing_scripts) > 0L) {
  stop(
    "The following required scripts are missing:\n  ",
    paste(missing_scripts, collapse = "\n  ")
  )
}

# -----------------------------------------------------------------------------
# 4. Execute each script in an isolated environment and record status
# -----------------------------------------------------------------------------

run_status <- data.frame(
  step = integer(),
  script = character(),
  description = character(),
  status = character(),
  started = character(),
  finished = character(),
  elapsed_seconds = double(),
  message = character(),
  stringsAsFactors = FALSE
)

status_file <- file.path(log_directory, "run_all_status.csv")

write_status <- function() {
  utils::write.csv(run_status, status_file, row.names = FALSE, na = "")
}

cat(
  "\nCorrected nine-module Cu/Fe methylation workflow\n",
  "Project root: ", project_root, "\n",
  "AD preprocessing: ", ad_preprocessing_choice, "\n",
  "g:Profiler mode: ", gprofiler_choice, "\n\n",
  sep = ""
)

for (row_number in seq_len(nrow(workflow))) {
  current <- workflow[row_number, ]
  start_time <- Sys.time()

  if (!current$run) {
    cat(
      sprintf(
        "[%02d/12] SKIPPED: %s\n           %s\n",
        current$step, current$script,
        "Resource-intensive AD preprocessing was not requested."
      )
    )
    run_status <- rbind(
      run_status,
      data.frame(
        step = current$step,
        script = current$script,
        description = current$description,
        status = "SKIPPED",
        started = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
        finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        elapsed_seconds = 0,
        message = "AD preprocessing not requested",
        stringsAsFactors = FALSE
      )
    )
    write_status()
    next
  }

  cat(
    sprintf(
      "[%02d/12] RUNNING: %s\n           %s\n",
      current$step, current$script, current$description
    )
  )

  result <- tryCatch({
    script_environment <- new.env(parent = globalenv())
    sys.source(
      file.path(project_root, current$script),
      envir = script_environment,
      chdir = FALSE
    )
    list(success = TRUE, message = "Completed")
  }, error = function(error) {
    list(success = FALSE, message = conditionMessage(error))
  })

  finish_time <- Sys.time()
  elapsed <- as.numeric(difftime(finish_time, start_time, units = "secs"))

  run_status <- rbind(
    run_status,
    data.frame(
      step = current$step,
      script = current$script,
      description = current$description,
      status = if (result$success) "COMPLETED" else "FAILED",
      started = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
      finished = format(finish_time, "%Y-%m-%d %H:%M:%S %Z"),
      elapsed_seconds = round(elapsed, 2),
      message = result$message,
      stringsAsFactors = FALSE
    )
  )
  write_status()

  if (!result$success) {
    cat(sprintf("           FAILED after %.1f seconds\n", elapsed))
    stop(
      "Workflow stopped at ", current$script, ": ", result$message,
      "\nReview ", status_file
    )
  }

  cat(sprintf("           COMPLETED in %.1f seconds\n", elapsed))
}

cat(
  "\nAll requested workflow steps completed successfully.\n",
  "Run report: ", status_file, "\n",
  "Final validation: ",
  file.path(project_root, "outputs", "validation", "final_validation_report.csv"),
  "\n",
  sep = ""
)
