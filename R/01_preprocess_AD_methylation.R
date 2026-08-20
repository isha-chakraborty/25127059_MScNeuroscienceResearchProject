# ==============================================================================
# 01_preprocess_AD_methylation.R
# Preprocessing and normalisation of GSE125895 AD DLPFC methylation data
#
# Input:
#   Illumina HumanMethylation450 IDAT files and sample sheet for GSE125895
#
# Main steps:
#   1. Import raw IDAT files
#   2. Assess sample quality and methylated/unmethylated signal intensities
#   3. Compare methylation-predicted and recorded sex
#   4. Filter low-quality and potentially ambiguous probes
#   5. Inspect beta-value distributions and sample clustering
#   6. Perform BMIQ normalisation
#   7. Repeat quality-control assessment
#   8. Examine major methylation variation using SVD
#
# Outputs:
#   BMIQ-normalised beta-value matrix
#   Aligned sample metadata
#   Quality-control tables and plots
# ==============================================================================

# Load shared setup -------------------------------------------------------------

if (!exists("paths", inherits = TRUE)) {
  source(file.path("R", "00_setup.R"))
}

check_packages(
  c(
    "ChAMP",
    "minfi",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  ),
  analysis_name = "AD methylation preprocessing"
)

# Input and output directories --------------------------------------------------

# Preferred clean-repository location
raw_data_directory <- paths$data_raw

# Temporary compatibility with the current project organisation
legacy_raw_directory <- file.path(project_root, "data_raw")

contains_idat_files <- function(directory) {
  dir.exists(directory) &&
    length(
      list.files(
        directory,
        pattern = "_(Red|Grn)\\.idat$",
        ignore.case = TRUE
      )
    ) > 0L
}

if (!contains_idat_files(raw_data_directory) &&
    contains_idat_files(legacy_raw_directory)) {
  raw_data_directory <- legacy_raw_directory

  message(
    "Using the existing raw-data directory: ",
    raw_data_directory
  )
}

if (!contains_idat_files(raw_data_directory)) {
  stop(
    paste0(
      "No IDAT files were found.\n",
      "Place the GSE125895 IDAT files and sample sheet in:\n",
      paths$data_raw
    ),
    call. = FALSE
  )
}

preprocessing_output_directory <- file.path(
  paths$outputs,
  "AD_methylation_preprocessing"
)

pre_qc_directory <- file.path(
  preprocessing_output_directory,
  "QC_pre_normalisation"
)

bmiq_directory <- file.path(
  preprocessing_output_directory,
  "BMIQ_normalisation"
)

post_qc_directory <- file.path(
  preprocessing_output_directory,
  "QC_post_normalisation"
)

svd_directory <- file.path(
  preprocessing_output_directory,
  "SVD"
)

invisible(
  lapply(
    c(
      preprocessing_output_directory,
      pre_qc_directory,
      bmiq_directory,
      post_qc_directory,
      svd_directory
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# Locate and inspect the sample sheet -------------------------------------------

sample_sheet_files <- list.files(
  raw_data_directory,
  pattern = "sample.*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(sample_sheet_files) != 1L) {
  stop(
    paste0(
      "Expected one sample-sheet CSV in ",
      raw_data_directory,
      " but found ",
      length(sample_sheet_files),
      "."
    ),
    call. = FALSE
  )
}

sample_sheet <- readr::read_csv(
  sample_sheet_files,
  show_col_types = FALSE
)

required_metadata_columns <- c(
  "Sample_Name",
  "Sample_Group",
  "Sex"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(sample_sheet)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    paste0(
      "The sample sheet is missing: ",
      paste(missing_metadata_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

if ("Region" %in% names(sample_sheet)) {
  observed_regions <- unique(
    stringr::str_to_upper(
      stringr::str_trim(sample_sheet$Region)
    )
  )

  if (!all(observed_regions == "DLPFC")) {
    warning(
      paste0(
        "The sample sheet contains regions other than DLPFC: ",
        paste(observed_regions, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

readr::write_csv(
  sample_sheet,
  file.path(
    preprocessing_output_directory,
    "GSE125895_DLPFC_sample_sheet_used.csv"
  )
)

# Initial import for sample-level quality control -------------------------------

message("Importing raw methylation data for initial quality control.")

raw_load <- ChAMP::champ.load(
  directory = raw_data_directory,
  method = "minfi",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = FALSE,
  ProbeCutoff = 0,
  SampleCutoff = 0.10,
  detPcut = 0.01,
  filterBeads = TRUE,
  beadCutoff = 0.05,
  filterNoCG = TRUE,
  filterSNPs = TRUE,
  population = NULL,
  filterMultiHit = TRUE,
  filterXY = FALSE,
  force = FALSE,
  arraytype = "450K"
)

raw_rgset <- raw_load$rgSet

if (is.null(raw_rgset)) {
  stop(
    "ChAMP did not return the expected RGChannelSet object.",
    call. = FALSE
  )
}

raw_mset <- minfi::preprocessRaw(raw_rgset)

raw_rset <- minfi::ratioConvert(
  raw_mset,
  what = "both",
  keepCN = TRUE
)

raw_grset <- minfi::mapToGenome(raw_rset)

raw_metadata <- Biobase::pData(raw_rgset) |>
  as.data.frame()

# Sample-level signal quality ---------------------------------------------------

qc_metrics <- minfi::getQC(raw_mset) |>
  as.data.frame()

qc_metrics$Sample_Name <- rownames(qc_metrics)

qc_metrics <- qc_metrics |>
  dplyr::relocate(Sample_Name)

readr::write_csv(
  qc_metrics,
  file.path(
    pre_qc_directory,
    "minfi_sample_signal_QC_metrics.csv"
  )
)

grDevices::png(
  filename = file.path(
    pre_qc_directory,
    "minfi_sample_signal_QC.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

minfi::plotQC(
  minfi::getQC(raw_mset)
)

grDevices::dev.off()

# Control-probe performance -----------------------------------------------------

grDevices::png(
  filename = file.path(
    pre_qc_directory,
    "bisulfite_conversion_control_probes.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

minfi::controlStripPlot(
  raw_rgset,
  controls = "BISULFITE CONVERSION II"
)

grDevices::dev.off()

# Raw methylation distributions ------------------------------------------------

sample_groups <- NULL

if ("Sample_Group" %in% names(raw_metadata)) {
  sample_groups <- raw_metadata$Sample_Group
} else {
  metadata_match <- match(
    rownames(raw_metadata),
    sample_sheet$Sample_Name
  )

  sample_groups <- sample_sheet$Sample_Group[metadata_match]
}

grDevices::pdf(
  file = file.path(
    pre_qc_directory,
    "raw_beta_density_beanplot.pdf"
  ),
  width = 10,
  height = 7
)

minfi::densityBeanPlot(
  raw_mset,
  sampGroups = sample_groups
)

grDevices::dev.off()

grDevices::pdf(
  file = file.path(
    pre_qc_directory,
    "raw_beta_density_plot.pdf"
  ),
  width = 10,
  height = 7
)

minfi::densityPlot(
  raw_mset,
  sampGroups = sample_groups
)

grDevices::dev.off()

# Predicted and recorded sex ----------------------------------------------------

predicted_sex_object <- minfi::getSex(
  raw_grset,
  cutoff = -2
)

predicted_sex <- tibble::tibble(
  Sample_Name = rownames(predicted_sex_object),
  Predicted_sex = as.character(
    predicted_sex_object$predictedSex
  )
)

sex_comparison <- sample_sheet |>
  dplyr::select(
    Sample_Name,
    Recorded_sex = Sex
  ) |>
  dplyr::left_join(
    predicted_sex,
    by = "Sample_Name"
  ) |>
  dplyr::mutate(
    Recorded_sex = stringr::str_to_upper(
      stringr::str_sub(
        stringr::str_trim(Recorded_sex),
        1,
        1
      )
    ),
    Predicted_sex = stringr::str_to_upper(
      stringr::str_sub(
        stringr::str_trim(Predicted_sex),
        1,
        1
      )
    ),
    Sex_match = Recorded_sex == Predicted_sex
  )

readr::write_csv(
  sex_comparison,
  file.path(
    pre_qc_directory,
    "recorded_and_methylation_predicted_sex.csv"
  )
)

sex_mismatches <- sex_comparison |>
  dplyr::filter(
    !is.na(Sex_match),
    !Sex_match
  )

if (nrow(sex_mismatches) == 0L) {
  message("No recorded-versus-predicted sex mismatches were detected.")
} else {
  warning(
    paste0(
      nrow(sex_mismatches),
      " recorded-versus-predicted sex mismatch(es) were detected. ",
      "Review the exported comparison table."
    ),
    call. = FALSE
  )
}

# ChAMP filtering for downstream normalisation ---------------------------------

message("Applying the prespecified probe- and sample-quality filters.")

filtered_load <- ChAMP::champ.load(
  directory = raw_data_directory,
  method = "ChAMP",
  methValue = "B",
  autoimpute = TRUE,
  filterDetP = TRUE,
  ProbeCutoff = 0,
  SampleCutoff = 0.10,
  detPcut = 0.01,
  filterBeads = TRUE,
  beadCutoff = 0.05,
  filterNoCG = TRUE,
  filterSNPs = TRUE,
  population = NULL,
  filterMultiHit = TRUE,
  filterXY = TRUE,
  force = FALSE,
  arraytype = "450K"
)

filtered_beta <- filtered_load$beta

if (is.null(filtered_beta) || !is.matrix(filtered_beta)) {
  stop(
    "ChAMP did not return the expected filtered beta-value matrix.",
    call. = FALSE
  )
}

filtered_metadata <- filtered_load$pd |>
  as.data.frame()

# Align metadata and beta values explicitly -------------------------------------

if ("Sample_Name" %in% names(filtered_metadata)) {
  rownames(filtered_metadata) <- as.character(
    filtered_metadata$Sample_Name
  )
}

shared_samples <- intersect(
  colnames(filtered_beta),
  rownames(filtered_metadata)
)

if (length(shared_samples) == 0L) {
  stop(
    "No shared sample identifiers were found between beta values and metadata.",
    call. = FALSE
  )
}

filtered_beta <- filtered_beta[, shared_samples, drop = FALSE]

filtered_metadata <- filtered_metadata[
  shared_samples,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(filtered_beta),
    rownames(filtered_metadata)
  )
)

# Quality control before BMIQ normalisation -------------------------------------

message("Generating pre-normalisation ChAMP quality-control plots.")

ChAMP::champ.QC(
  resultsDir = pre_qc_directory,
  beta = filtered_beta,
  densityPlot = TRUE,
  mdsPlot = TRUE,
  dendrogram = TRUE,
  pheno = filtered_metadata$Sample_Group
)

# BMIQ normalisation ------------------------------------------------------------

message("Performing beta-mixture quantile normalisation.")

filtered_mset <- filtered_load$mset

if (is.null(filtered_mset)) {
  filtered_mset <- filtered_load$mSet
}

normalised_beta <- ChAMP::champ.norm(
  beta = filtered_beta,
  rgSet = filtered_load$rgSet,
  mset = filtered_mset,
  resultsDir = bmiq_directory,
  method = "BMIQ",
  plotBMIQ = TRUE,
  arraytype = "450K",
  cores = 1
)

if (is.null(normalised_beta) || !is.matrix(normalised_beta)) {
  stop(
    "BMIQ normalisation did not return the expected beta-value matrix.",
    call. = FALSE
  )
}

if (anyNA(normalised_beta)) {
  warning(
    "The normalised beta-value matrix contains missing values.",
    call. = FALSE
  )
}

if (any(normalised_beta < 0 | normalised_beta > 1, na.rm = TRUE)) {
  stop(
    "Normalised beta values outside the expected 0–1 range were detected.",
    call. = FALSE
  )
}

# Quality control after BMIQ normalisation --------------------------------------

message("Generating post-normalisation ChAMP quality-control plots.")

ChAMP::champ.QC(
  resultsDir = post_qc_directory,
  beta = normalised_beta,
  densityPlot = TRUE,
  mdsPlot = TRUE,
  dendrogram = TRUE,
  pheno = filtered_metadata$Sample_Group
)

# Singular-value decomposition -------------------------------------------------

# M-values are used for the SVD because they are more suitable than beta values
# for statistical assessment of methylation variation.

epsilon <- 1e-6

bounded_beta <- pmin(
  pmax(normalised_beta, epsilon),
  1 - epsilon
)

normalised_m_values <- log2(
  bounded_beta / (1 - bounded_beta)
)

svd_metadata <- filtered_metadata[
  colnames(normalised_m_values),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(normalised_m_values),
    rownames(svd_metadata)
  )
)

message("Examining major sources of methylation variation using SVD.")

grDevices::png(
  filename = file.path(
    svd_directory,
    "GSE125895_DLPFC_SVD_post_BMIQ.png"
  ),
  width = 3600,
  height = 3000,
  res = 300
)

graphics::par(
  mar = c(5, 12, 4, 2)
)

ChAMP::champ.SVD(
  beta = as.data.frame(normalised_m_values),
  pd = svd_metadata
)

grDevices::dev.off()

grDevices::pdf(
  file = file.path(
    svd_directory,
    "GSE125895_DLPFC_SVD_post_BMIQ.pdf"
  ),
  width = 12,
  height = 10
)

graphics::par(
  mar = c(5, 12, 4, 2)
)

ChAMP::champ.SVD(
  beta = as.data.frame(normalised_m_values),
  pd = svd_metadata
)

grDevices::dev.off()

# Save processed objects --------------------------------------------------------

normalised_beta_file <- file.path(
  paths$data_processed,
  "GSE125895_DLPFC_BMIQ_normalised_beta.rds"
)

normalised_m_file <- file.path(
  paths$data_processed,
  "GSE125895_DLPFC_BMIQ_normalised_M_values.rds"
)

metadata_file <- file.path(
  paths$data_processed,
  "GSE125895_DLPFC_processed_metadata.rds"
)

saveRDS(
  normalised_beta,
  normalised_beta_file,
  compress = FALSE
)

saveRDS(
  normalised_m_values,
  normalised_m_file,
  compress = FALSE
)

saveRDS(
  filtered_metadata,
  metadata_file
)

readr::write_csv(
  tibble::rownames_to_column(
    filtered_metadata,
    var = "Methylation_sample_ID"
  ),
  file.path(
    paths$data_processed,
    "GSE125895_DLPFC_processed_metadata.csv"
  )
)

# Preprocessing summary ---------------------------------------------------------

preprocessing_summary <- tibble::tibble(
  measure = c(
    "Samples in submitted sample sheet",
    "Samples after ChAMP filtering",
    "Probes after ChAMP filtering",
    "Samples after BMIQ normalisation",
    "Probes after BMIQ normalisation",
    "Recorded/predicted sex mismatches",
    "Minimum normalised beta value",
    "Maximum normalised beta value"
  ),
  value = c(
    nrow(sample_sheet),
    ncol(filtered_beta),
    nrow(filtered_beta),
    ncol(normalised_beta),
    nrow(normalised_beta),
    nrow(sex_mismatches),
    min(normalised_beta, na.rm = TRUE),
    max(normalised_beta, na.rm = TRUE)
  )
)

readr::write_csv(
  preprocessing_summary,
  file.path(
    preprocessing_output_directory,
    "GSE125895_DLPFC_preprocessing_summary.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    preprocessing_output_directory,
    "GSE125895_DLPFC_preprocessing_session_info.txt"
  )
)

message("AD DLPFC methylation preprocessing completed successfully.")
message("Normalised beta values: ", normalised_beta_file)
message("Normalised M-values: ", normalised_m_file)
message("Processed metadata: ", metadata_file)
