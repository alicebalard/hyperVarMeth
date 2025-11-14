## data preparation
## Alice Balard
## Last major update Nov 2025

#' Prepare input data for hvCpG analysis
#'
#' This function loads the metadata, dataset-level parameters, and CpG names
#' required to run the hvCpG detection algorithm.
#'
#' @param analysis Character string. Name of the analysis.
#'   If it contains `"MariasarraysREDUCED"`, a special directory structure is expected.
#' @param dataDir Character string. Path to the directory containing input data files.
#'   Should contain `sample_metadata.tsv`, `all_medsd_lambda.tsv`, and an HDF5 matrix file.
#' @param subsetMetadata Logical or data frame. If `FALSE` (default), the full metadata
#'   is used. Otherwise, a subset of the metadata can be provided to restrict
#'   the analysis to specific samples or datasets.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{metadata}{Data frame containing `sample` and `dataset` columns.}
#'     \item{medsd_lambdas}{Data frame containing dataset-level parameters in
#'       `dataset`, `median_sd`, and `lambda` columns.}
#'     \item{cpg_names_all}{Character vector of CpG names (columns in the HDF5 matrix).}
#'     \item{h5file}{Character string with the path to the HDF5 file containing the methylation matrix (samples x CpG names).}
#'   }
#'
#' @export
#'
#' @importFrom rhdf5 h5read
#' @importFrom utils read.table
prepData <- function(analysis, dataDir, subsetMetadata = FALSE) {
  if (grepl("MariasarraysREDUCED", analysis)) {
    x <- sub("MariasarraysREDUCED", "", analysis)
    basepath <- file.path("/home/alice/arraysh5_reducedMimicAtlas", x)
    metapath <- file.path(basepath, "all_metadata.tsv")
    if (!file.exists(metapath)) {
      stop("Run 03_prepDatasetsMaria/S02.prepare_REDUCED_arrays_mimicAtlas.py")
    }
    metadata <- read.table(metapath, sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(basepath, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    cpg_names_all <- h5read(file.path(basepath, "all_scaled_matrix.h5"), "cpg_names")
    h5file <- file.path(basepath, "all_scaled_matrix.h5")
  } else {
    metadata <- read.table(file.path(dataDir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
    medsd_lambdas <- read.table(file.path(dataDir, "all_medsd_lambda.tsv"), sep = "\t", header = TRUE)
    cpg_names_all <- h5read(file.path(dataDir, "all_matrix_noscale.h5"), "cpg_names")
    h5file <- file.path(dataDir, "all_matrix_noscale.h5")
  }
  ## Possible subset of samples or datasets
  if (!isFALSE(subsetMetadata)) {
    message("The algorithm runs on a subset of metadata")
    metadata = subsetMetadata
  }
  return(list(
    metadata = metadata,
    medsd_lambdas = medsd_lambdas,
    cpg_names_all = cpg_names_all,
    h5file = h5file
  ))
}
