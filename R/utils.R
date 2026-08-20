# Shared utilities for the corrected Cu/Fe methylation analysis
#
# This file contains small, analysis-neutral helpers used across scripts. It does
# not run an analysis or modify data when sourced.

# -----------------------------------------------------------------------------
# Project paths and input checks
# -----------------------------------------------------------------------------

find_project_root <- function(start = getwd()) {
  candidate <- normalizePath(start, mustWork = TRUE)

  repeat {
    has_r_directory <- dir.exists(file.path(candidate, "R"))
    has_project_marker <-
      file.exists(file.path(candidate, "run_all.R")) ||
      length(list.files(candidate, pattern = "[.]Rproj$")) > 0L

    if (has_r_directory && has_project_marker) return(candidate)

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop(
    "Could not locate the project root. Run the script from inside the project ",
    "containing the R directory and project file."
  )
}

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

require_file <- function(path, description = "Required input") {
  if (!file.exists(path)) {
    stop(description, " was not found: ", path)
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

read_required_csv <- function(path, required_columns = character(), ...) {
  require_file(path, "Required CSV file")
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required to read: ", path)
  }

  data <- readr::read_csv(path, show_col_types = FALSE, ...)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      basename(path), " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  data
}

# -----------------------------------------------------------------------------
# Gene-symbol and logical-value standardisation
# -----------------------------------------------------------------------------

clean_gene_symbols <- function(x, sort_result = TRUE) {
  symbols <- toupper(trimws(as.character(x)))
  symbols <- unique(symbols[
    !is.na(symbols) & symbols != "" & symbols != "NA"
  ])
  if (sort_result) symbols <- sort(symbols)
  symbols
}

as_logical_flag <- function(x, missing_value = FALSE) {
  if (is.logical(x)) {
    x[is.na(x)] <- missing_value
    return(x)
  }

  normalised <- tolower(trimws(as.character(x)))
  output <- normalised %in% c("true", "t", "1", "yes", "y")
  output[is.na(normalised) | normalised == ""] <- missing_value
  output
}

# -----------------------------------------------------------------------------
# Stable module identifiers and display labels
# -----------------------------------------------------------------------------

standardise_module_name <- function(x) {
  original <- trimws(as.character(x))
  normalised <- gsub("[-[:space:]]+", "_", original)
  normalised <- gsub("_+", "_", normalised)
  normalised <- gsub("^_|_$", "", normalised)

  pieces <- strsplit(normalised, "_", fixed = TRUE)
  vapply(pieces, function(parts) {
    if (length(parts) == 0L) return(NA_character_)

    first <- toupper(parts[1])

    if (first == "AD" && length(parts) >= 3L) {
      return(paste("AD", toupper(parts[2]), tolower(parts[3]), sep = "_"))
    }
    if (grepl("^FTLD[123]$", first) && length(parts) >= 2L) {
      return(paste(first, tolower(parts[2]), sep = "_"))
    }
    if (first %in% c("MSA", "PD", "PSP", "ND") && length(parts) >= 2L) {
      return(paste(first, tolower(parts[2]), sep = "_"))
    }

    normalised_value <- paste(parts, collapse = "_")
    warning("Unrecognised module-name format retained as: ", normalised_value)
    normalised_value
  }, character(1))
}

display_module_name <- function(x) {
  gsub("_", " ", standardise_module_name(x), fixed = TRUE)
}

corrected_nine_module_ids <- function() {
  c(
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
}

# -----------------------------------------------------------------------------
# Fisher overrepresentation test
# -----------------------------------------------------------------------------

fisher_overrepresentation <- function(
    module_genes,
    gene_set,
    background_genes,
    alternative = "two.sided") {

  background <- clean_gene_symbols(background_genes)
  module <- intersect(clean_gene_symbols(module_genes), background)
  tested_set <- intersect(clean_gene_symbols(gene_set), background)

  overlap <- intersect(module, tested_set)
  module_only <- setdiff(module, tested_set)
  set_only <- setdiff(tested_set, module)
  neither <- setdiff(background, union(module, tested_set))

  contingency_table <- matrix(
    c(
      length(overlap),
      length(module_only),
      length(set_only),
      length(neither)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      module = c("Inside", "Outside"),
      metal_set = c("Member", "Not_member")
    )
  )

  test <- stats::fisher.test(
    contingency_table,
    alternative = alternative
  )

  data.frame(
    overlap = length(overlap),
    module_size = length(module),
    gene_set_size_in_background = length(tested_set),
    background_size = length(background),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value,
    overlap_genes = paste(sort(overlap), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

adjust_fdr <- function(p_values) {
  stats::p.adjust(p_values, method = "BH")
}

# -----------------------------------------------------------------------------
# Reproducible checks
# -----------------------------------------------------------------------------

assert_expected_count <- function(observed, expected, label) {
  observed <- as.integer(observed)
  expected <- as.integer(expected)
  if (length(observed) != 1L || is.na(observed) || observed != expected) {
    stop(label, ": expected ", expected, ", observed ", observed, ".")
  }
  invisible(TRUE)
}

assert_set_equal <- function(observed, expected, label) {
  observed <- unique(as.character(observed))
  expected <- unique(as.character(expected))
  if (!setequal(observed, expected)) {
    missing_values <- setdiff(expected, observed)
    unexpected_values <- setdiff(observed, expected)
    stop(
      label, " does not match the expected set. Missing: ",
      paste(missing_values, collapse = ", "), "; unexpected: ",
      paste(unexpected_values, collapse = ", "), "."
    )
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Colours and figure export
# -----------------------------------------------------------------------------

lighten_colour <- function(colour, amount = 0.5) {
  if (!is.numeric(amount) || length(amount) != 1L ||
      is.na(amount) || amount < 0 || amount > 1) {
    stop("amount must be a single number between 0 and 1.")
  }

  rgb_values <- grDevices::col2rgb(colour) / 255
  light_values <- rgb_values + (1 - rgb_values) * amount
  grDevices::rgb(light_values[1, ], light_values[2, ], light_values[3, ])
}

save_figure_set <- function(
    plot,
    output_directory,
    filename_stem,
    width,
    height,
    dpi = 600,
    units = "in",
    save_png = TRUE,
    save_tiff = TRUE,
    save_pdf = TRUE,
    background = "white") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to export figures.")
  }

  ensure_directory(output_directory)

  output_files <- character()

  if (save_png) {
    path <- file.path(output_directory, paste0(filename_stem, ".png"))
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height,
      units = units, dpi = dpi, bg = background
    )
    output_files <- c(output_files, path)
  }

  if (save_tiff) {
    path <- file.path(output_directory, paste0(filename_stem, ".tiff"))
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height,
      units = units, dpi = dpi, compression = "lzw", bg = background
    )
    output_files <- c(output_files, path)
  }

  if (save_pdf) {
    path <- file.path(output_directory, paste0(filename_stem, ".pdf"))
    pdf_device <- grDevices::pdf
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height,
      units = units, device = pdf_device, useDingbats = FALSE, bg = background
    )
    output_files <- c(output_files, path)
  }

  invisible(normalizePath(output_files, mustWork = TRUE))
}

# -----------------------------------------------------------------------------
# Table-safe conversion of list columns
# -----------------------------------------------------------------------------

flatten_list_columns <- function(data, separator = "; ") {
  list_columns <- names(data)[vapply(data, is.list, logical(1))]
  for (column in list_columns) {
    data[[column]] <- vapply(
      data[[column]],
      function(value) paste(as.character(value), collapse = separator),
      character(1)
    )
  }
  data
}
