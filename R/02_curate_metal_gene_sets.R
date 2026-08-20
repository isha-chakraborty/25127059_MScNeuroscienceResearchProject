# ==============================================================================
# 02_curate_metal_gene_sets.R
# Curation of copper, broad-iron and core-iron gene sets
#
# Purpose:
#   Reconstruct the three prespecified metal-associated gene sets from the
#   AmiGO annotation exports used in the thesis.
#
# Curation criteria:
#   - Homo sapiens annotations only
#   - Experimental or high-throughput experimental evidence
#   - Standardised, non-missing gene symbols
#   - Duplicate symbols removed
#   - Non-gene complexes, assemblies and microRNA entries removed
#
# Expected outputs:
#   47 copper genes
#   199 broad-iron genes
#   136 core-iron genes
#   Core iron must be a subset of broad iron
#
# Supplementary outputs:
#   Supplementary Table S1: copper gene-set curation
#   Supplementary Table S2: broad-iron gene-set curation
#   Supplementary Table S3: core-iron gene-set curation
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
  analysis_name = "metal gene-set curation"
)

# Locate the AmiGO source files -------------------------------------------------

# Preferred public-repository organisation:
#
# data/
# └── gene_sets/
#     ├── copper/
#     │   └── amigo/
#     └── iron/
#         └── amigo/
#
# The second location supports the existing project organisation during
# repository preparation.

candidate_gene_set_roots <- c(
  paths$gene_sets,
  file.path(project_root, "gene_lists")
)

is_valid_gene_set_root <- function(directory) {
  dir.exists(file.path(directory, "copper", "amigo")) &&
    dir.exists(file.path(directory, "iron", "amigo"))
}

valid_gene_set_roots <- candidate_gene_set_roots[
  vapply(
    candidate_gene_set_roots,
    is_valid_gene_set_root,
    FUN.VALUE = logical(1)
  )
]

if (length(valid_gene_set_roots) == 0L) {
  stop(
    paste0(
      "The AmiGO source folders could not be found.\n",
      "Expected the following folders beneath ",
      paths$gene_sets,
      ":\n",
      " - copper/amigo\n",
      " - iron/amigo"
    ),
    call. = FALSE
  )
}

gene_set_source_root <- valid_gene_set_roots[[1]]

copper_amigo_directory <- file.path(
  gene_set_source_root,
  "copper",
  "amigo"
)

iron_amigo_directory <- file.path(
  gene_set_source_root,
  "iron",
  "amigo"
)

message("Reading AmiGO files from: ", gene_set_source_root)

# Identify the final AmiGO exports ----------------------------------------------

# The "_3.txt" suffix identifies the final AmiGO downloads used in the thesis.

copper_amigo_files <- list.files(
  copper_amigo_directory,
  pattern = "_3\\.txt$",
  full.names = TRUE
) |>
  sort()

broad_iron_amigo_files <- list.files(
  iron_amigo_directory,
  pattern = "_3\\.txt$",
  full.names = TRUE
) |>
  sort()

if (length(copper_amigo_files) == 0L) {
  stop(
    "No final copper AmiGO files ending in '_3.txt' were found.",
    call. = FALSE
  )
}

if (length(broad_iron_amigo_files) == 0L) {
  stop(
    "No final iron AmiGO files ending in '_3.txt' were found.",
    call. = FALSE
  )
}

# The core-iron list uses terms directly related to iron binding, transport,
# import, homeostasis, cellular responses, ferroxidase activity and
# iron-sulfur-cluster assembly. Heme-related terms are retained in the broad
# list but are not included in the core list.

core_iron_source_filenames <- c(
  "cellular_response_to_iron_ion_3.txt",
  "ferric_iron_binding_3.txt",
  "ferrous_iron_binding_3.txt",
  "ferroxidase_activity_3.txt",
  "intracellular_iron_ion_homeostasis_3.txt",
  "iron_homeostasis_3.txt",
  "iron_import_into_cell_3.txt",
  "iron_ion_binding_3.txt",
  "iron_ion_transmembrane_transport_3.txt",
  "iron_ion_transmembrane_transporter_activity_3.txt",
  "iron_ion_transport_3.txt",
  "iron-sulfur_cluster_assembly_3.txt",
  "multicellular_organismal-level_iron_ion_homeostasis_3.txt",
  "response_to_iron_ion_3.txt"
)

core_iron_amigo_files <- file.path(
  iron_amigo_directory,
  core_iron_source_filenames
)

check_input_files(
  core_iron_amigo_files,
  description = "core-iron AmiGO source files"
)

# Evidence-code definitions -----------------------------------------------------

experimental_evidence_codes <- c(
  "EXP",
  "IDA",
  "IMP",
  "IGI",
  "IEP",
  "IPI",
  "HDA"
)

evidence_code_information <- tibble::tribble(
  ~evidence_code, ~description,
  "EXP", "Inferred from experiment",
  "IDA", "Inferred from direct assay",
  "IMP", "Inferred from mutant phenotype",
  "IGI", "Inferred from genetic interaction",
  "IEP", "Inferred from expression pattern",
  "IPI", "Inferred from physical interaction",
  "HDA", "Inferred from high-throughput direct assay"
)

# AmiGO/Gene Ontology annotation reader ----------------------------------------

# The downloaded files follow Gene Ontology annotation-file column order.
# Columns required here are:
#   column 3:  gene symbol
#   column 5:  Gene Ontology identifier
#   column 7:  evidence code
#   column 13: organism/taxon

gaf_column_names <- c(
  "Database",
  "Database_object_ID",
  "Gene_symbol",
  "Qualifier",
  "GO_ID",
  "Database_reference",
  "Evidence_code",
  "With_or_from",
  "Aspect",
  "Database_object_name",
  "Synonym",
  "Database_object_type",
  "Taxon",
  "Date",
  "Assigned_by",
  "Annotation_extension",
  "Gene_product_form_ID"
)

source_term_from_filename <- function(filename) {
  basename(filename) |>
    stringr::str_remove("_3\\.txt$") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_squish()
}

read_amigo_annotation <- function(filename) {
  annotation <- suppressWarnings(
    utils::read.delim(
      filename,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "!",
      fill = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  if (ncol(annotation) < 13L) {
    stop(
      paste0(
        "The AmiGO file does not contain the expected annotation columns: ",
        basename(filename)
      ),
      call. = FALSE
    )
  }

  names(annotation) <- gaf_column_names[
    seq_len(
      min(
        ncol(annotation),
        length(gaf_column_names)
      )
    )
  ]

  annotation |>
    dplyr::filter(
      stringr::str_detect(
        as.character(Taxon),
        "(^|:)9606($|\\|)"
      ),
      Evidence_code %in% experimental_evidence_codes
    ) |>
    dplyr::transmute(
      geneSymbol = standardise_gene_symbols(Gene_symbol),
      GO_ID = as.character(GO_ID),
      evidence_code = as.character(Evidence_code),
      taxon = as.character(Taxon),
      source_term = source_term_from_filename(filename),
      source_file = basename(filename),
      source_database = "AmiGO/Gene Ontology"
    )
}

# Clean annotation-level records ------------------------------------------------

clean_amigo_annotations <- function(files) {
  purrr::map_dfr(
    files,
    read_amigo_annotation
  ) |>
    dplyr::filter(
      !is.na(geneSymbol),
      geneSymbol != "",
      geneSymbol != "GENE_SYMBOL",

      # Remove protein complexes or assemblies represented by labels such as
      # "rrm1-rrm2_human" or "bola2-glrx3_human".
      !stringr::str_detect(
        geneSymbol,
        "_HUMAN$"
      ),

      # Remove microRNA annotations because downstream network enrichment uses
      # standard gene-symbol assignments.
      !stringr::str_detect(
        geneSymbol,
        stringr::regex(
          "MIR",
          ignore_case = TRUE
        )
      )
    ) |>
    dplyr::distinct(
      geneSymbol,
      GO_ID,
      evidence_code,
      source_term,
      source_file,
      .keep_all = TRUE
    ) |>
    dplyr::arrange(
      geneSymbol,
      source_term,
      evidence_code
    )
}

copper_annotations <- clean_amigo_annotations(
  copper_amigo_files
)

broad_iron_annotations <- clean_amigo_annotations(
  broad_iron_amigo_files
)

core_iron_annotations <- clean_amigo_annotations(
  core_iron_amigo_files
)

# Construct the final unique gene sets ------------------------------------------

copper_genes <- copper_annotations |>
  dplyr::distinct(geneSymbol) |>
  dplyr::arrange(geneSymbol) |>
  dplyr::pull(geneSymbol)

broad_iron_genes <- broad_iron_annotations |>
  dplyr::distinct(geneSymbol) |>
  dplyr::arrange(geneSymbol) |>
  dplyr::pull(geneSymbol)

core_iron_genes <- core_iron_annotations |>
  dplyr::distinct(geneSymbol) |>
  dplyr::arrange(geneSymbol) |>
  dplyr::pull(geneSymbol)

# Validate expected gene-set properties -----------------------------------------

observed_gene_set_sizes <- c(
  copper = length(copper_genes),
  broad_iron = length(broad_iron_genes),
  core_iron = length(core_iron_genes)
)

expected_gene_set_sizes <- analysis_constants$expected_gene_set_sizes

incorrect_sizes <- names(observed_gene_set_sizes)[
  observed_gene_set_sizes !=
    expected_gene_set_sizes[names(observed_gene_set_sizes)]
]

if (length(incorrect_sizes) > 0L) {
  size_message <- paste(
    paste0(
      incorrect_sizes,
      ": observed ",
      observed_gene_set_sizes[incorrect_sizes],
      ", expected ",
      expected_gene_set_sizes[incorrect_sizes]
    ),
    collapse = "; "
  )

  stop(
    paste0(
      "The reconstructed gene-set sizes do not match the thesis: ",
      size_message
    ),
    call. = FALSE
  )
}

core_not_in_broad <- setdiff(
  core_iron_genes,
  broad_iron_genes
)

if (length(core_not_in_broad) > 0L) {
  stop(
    paste0(
      "The following core-iron genes were absent from the broad-iron set: ",
      paste(core_not_in_broad, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (anyDuplicated(copper_genes) > 0L ||
    anyDuplicated(broad_iron_genes) > 0L ||
    anyDuplicated(core_iron_genes) > 0L) {
  stop(
    "Duplicated symbols remain in at least one final gene set.",
    call. = FALSE
  )
}

# Optional comparison with the original final lists -----------------------------

original_final_files <- c(
  copper = file.path(
    gene_set_source_root,
    "Cu_AmiGO_gene_list_3.csv"
  ),
  broad_iron = file.path(
    gene_set_source_root,
    "Fe_AmiGO_gene_list_3.csv"
  ),
  core_iron = file.path(
    gene_set_source_root,
    "Fe_core_AmiGO_gene_list_3.csv"
  )
)

reconstructed_gene_sets <- list(
  copper = copper_genes,
  broad_iron = broad_iron_genes,
  core_iron = core_iron_genes
)

if (all(file.exists(original_final_files))) {
  original_gene_sets <- lapply(
    original_final_files,
    function(filename) {
      original <- readr::read_csv(
        filename,
        show_col_types = FALSE
      )

      if (!"geneSymbol" %in% names(original)) {
        stop(
          paste0(
            "The original list lacks a geneSymbol column: ",
            filename
          ),
          call. = FALSE
        )
      }

      original$geneSymbol |>
        standardise_gene_symbols() |>
        unique() |>
        sort()
    }
  )

  for (gene_set_name in names(reconstructed_gene_sets)) {
    if (!setequal(
      reconstructed_gene_sets[[gene_set_name]],
      original_gene_sets[[gene_set_name]]
    )) {
      stop(
        paste0(
          "The reconstructed ",
          gene_set_name,
          " list differs from the original final list."
        ),
        call. = FALSE
      )
    }
  }

  message(
    "All reconstructed gene sets match the original final lists exactly."
  )
}

# Collapse provenance to one row per gene ---------------------------------------

collapse_unique_values <- function(values) {
  paste(
    sort(
      unique(
        values[
          !is.na(values) &
            values != ""
        ]
      )
    ),
    collapse = "; "
  )
}

make_gene_curation_table <- function(
    annotations,
    gene_set_label) {

  annotations |>
    dplyr::group_by(geneSymbol) |>
    dplyr::summarise(
      Gene_set = gene_set_label,
      GO_terms = collapse_unique_values(source_term),
      GO_identifiers = collapse_unique_values(GO_ID),
      Evidence_codes = collapse_unique_values(evidence_code),
      Source_files = collapse_unique_values(source_file),
      .groups = "drop"
    ) |>
    dplyr::arrange(geneSymbol) |>
    dplyr::mutate(
      No. = dplyr::row_number(),
      .before = 1
    )
}

supplementary_table_s1 <- make_gene_curation_table(
  copper_annotations,
  "Copper"
)

supplementary_table_s2 <- make_gene_curation_table(
  broad_iron_annotations,
  "Broad iron"
)

supplementary_table_s3 <- make_gene_curation_table(
  core_iron_annotations,
  "Core iron"
)

# Source-term manifest ----------------------------------------------------------

source_manifest <- dplyr::bind_rows(
  tibble::tibble(
    Metal_category = "Copper",
    Source_file = basename(copper_amigo_files),
    Source_term = source_term_from_filename(
      copper_amigo_files
    ),
    Included_in_broad_set = TRUE,
    Included_in_core_set = NA
  ),
  tibble::tibble(
    Metal_category = "Iron",
    Source_file = basename(broad_iron_amigo_files),
    Source_term = source_term_from_filename(
      broad_iron_amigo_files
    ),
    Included_in_broad_set = TRUE,
    Included_in_core_set =
      basename(broad_iron_amigo_files) %in%
      core_iron_source_filenames
  )
) |>
  dplyr::arrange(
    Metal_category,
    Source_term
  )

# Save final gene lists ---------------------------------------------------------

dir.create(
  paths$gene_sets,
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  tibble::tibble(
    geneSymbol = copper_genes
  ),
  file.path(
    paths$gene_sets,
    "copper_gene_set.csv"
  )
)

readr::write_csv(
  tibble::tibble(
    geneSymbol = broad_iron_genes
  ),
  file.path(
    paths$gene_sets,
    "broad_iron_gene_set.csv"
  )
)

readr::write_csv(
  tibble::tibble(
    geneSymbol = core_iron_genes
  ),
  file.path(
    paths$gene_sets,
    "core_iron_gene_set.csv"
  )
)

metal_gene_sets <- list(
  copper = copper_genes,
  broad_iron = broad_iron_genes,
  core_iron = core_iron_genes
)

saveRDS(
  metal_gene_sets,
  file.path(
    paths$gene_sets,
    "metal_gene_sets.rds"
  )
)

# Shared downstream object used by candidate, functional and validation scripts.
saveRDS(
  metal_gene_sets,
  file.path(
    paths$derived_data,
    "curated_metal_gene_sets.rds"
  )
)

# Save supplementary source tables ---------------------------------------------

readr::write_csv(
  supplementary_table_s1,
  file.path(
    paths$supplementary_tables,
    "Supplementary_Table_S1_copper_gene_set_curation.csv"
  )
)

readr::write_csv(
  supplementary_table_s2,
  file.path(
    paths$supplementary_tables,
    "Supplementary_Table_S2_broad_iron_gene_set_curation.csv"
  )
)

readr::write_csv(
  supplementary_table_s3,
  file.path(
    paths$supplementary_tables,
    "Supplementary_Table_S3_core_iron_gene_set_curation.csv"
  )
)

readr::write_csv(
  source_manifest,
  file.path(
    paths$gene_sets,
    "AmiGO_source_term_manifest.csv"
  )
)

readr::write_csv(
  evidence_code_information,
  file.path(
    paths$gene_sets,
    "AmiGO_experimental_evidence_codes.csv"
  )
)

# Save curation summary ---------------------------------------------------------

curation_summary <- tibble::tibble(
  Gene_set = c(
    "Copper",
    "Broad iron",
    "Core iron"
  ),
  Number_of_source_terms = c(
    length(copper_amigo_files),
    length(broad_iron_amigo_files),
    length(core_iron_amigo_files)
  ),
  Number_of_unique_genes = c(
    length(copper_genes),
    length(broad_iron_genes),
    length(core_iron_genes)
  ),
  Human_taxon = "NCBITaxon:9606",
  Evidence_codes = paste(
    experimental_evidence_codes,
    collapse = "; "
  )
)

readr::write_csv(
  curation_summary,
  file.path(
    paths$gene_sets,
    "metal_gene_set_curation_summary.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    paths$logs,
    "02_metal_gene_set_curation_session_info.txt"
  )
)

message("Metal gene-set curation completed successfully.")
message("Copper genes: ", length(copper_genes))
message("Broad-iron genes: ", length(broad_iron_genes))
message("Core-iron genes: ", length(core_iron_genes))
message("Core iron is a subset of broad iron: TRUE")
