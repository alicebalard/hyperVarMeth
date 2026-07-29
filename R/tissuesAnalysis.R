## hv probability comparison by tissue
## Alice Balard
## Last major update Nov 2025

#' Compute log posterior probability of variability per dataset (CpG-level)
#'
#' For a given CpG, computes the **log posterior probability** that it is
#' *hypervariable*, separately for each dataset.
#' Within each dataset, we take the average over all sampled individuals to ensure comparability
#' across datasets with differing numbers of samples.
#'
#' @param Mdf Numeric matrix (1 CpG x samples). Methylation values for a single CpG.
#' @param metadata Data frame containing `sample` and `dataset` columns.
#' @param dataset_groups Named list mapping dataset names to row indices in `metadata`.
#' @param ds_params Data frame with precomputed parameters per dataset:
#'   columns `sd0`, `sd1`, and rownames corresponding to dataset names.
#' @param p1 Numeric scalars: True positive rate.
#'
#' @return A named numeric vector giving the log posterior probability (`log(Phv|Dk)`)
#'   for each dataset where the CpG has valid data.
#' @export
#'
#' @importFrom stats dnorm
getLogPhv_oneCpG_byTissue <- function(Mdf, metadata, dataset_groups, ds_params, p1) {
  # Must be a numeric named vector
  if (!is.numeric(Mdf))
    stop("Mdf must be a numeric vector")

  if (is.null(names(Mdf)))
    stop("Mdf must have names = sample IDs")

  datasets <- unique(metadata$dataset)
  out <- setNames(rep(NA_real_, length(datasets)), datasets)

  ## Run for each dataset
  for (k in datasets) {

    ## Safe extraction of the correct individual values
    samp_ids <- metadata$sample[dataset_groups[[k]]]
    samp_ids <- intersect(samp_ids, names(Mdf))
    if (length(samp_ids) == 0)
      next

    vals <- Mdf[samp_ids]
    vals_ok <- vals[!is.na(vals)]

    # Skip dataset if insufficient usable values
    if (length(vals_ok) < 3 || all(is.na(vals_ok))) {
      next
    }

    ## Precompute mean
    mu <- mean(vals_ok)

    ## Call precomputed sds for both cases
    params =  ds_params[k, ]

    p_stable <- dnorm(vals_ok, mu, params$sd0)
    p_var <- dnorm(vals_ok, mu, params$sd1)

    denom <- p_var * p1 + p_stable * (1 - p1)
    valid <- denom > 0
    if (!any(valid)) next

    ## Compute log(Phv|Dk) as the sum of logs on all individuals
    post <- (p_var[valid] * p1) / denom[valid]
    post <- pmin(pmax(post, 1e-12), 1 - 1e-12) ## avoid 0 and 1

    out[k] <- mean(log(post)) ## take the mean out of all individuals in one dataset
  }
  return(out) # named vector of a logPr per dataset for this CpG
}

#' Compute log(Phv|Dk) across all CpGs in batches (HDF5-based)
#'
#' Loads methylation values from an HDF5 file in batches and computes,
#' for each CpG, the log posterior probability that it is hypervariable
#' in each dataset. Supports parallel execution across multiple cores.
#'
#' @param cpg_names_vec Character vector of CpG identifiers to process.
#' @param NCORES Integer; number of CPU cores to use in parallel processing.
#' @param p1 Numeric scalars: True positive rate.
#' @param prep List returned by [prepData()], containing metadata and HDF5 paths.
#' @param batch_size Integer; number of CpGs to process per batch (default 10000).
#' @param Nds Integer; minimum number of datasets required per CpG (default 3).
#'
#' @return A numeric matrix of dimensions `(CpG x dataset)`, containing
#'   log posterior probabilities (`log(Phv|Dk)`) for each CpG and dataset.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rhdf5 h5read
#' @importFrom stats setNames
#' @importFrom parallel mclapply
getALL_LogPhv_oneCpG_byTissue_batch <- function(
    cpg_names_vec, NCORES, p1, prep,
    batch_size = 1000L, Nds = 3L) {

  metadata       <- prep$metadata
  cpg_names_all  <- prep$cpg_names_all
  h5file         <- prep$h5file
  medsd_lambdas  <- prep$medsd_lambdas

  ## dataset params
  ds_params <- medsd_lambdas %>%
    dplyr::select(dataset, median_sd, lambda) %>%
    dplyr::mutate(sd0 = pmax(median_sd, 0.005),
                  sd1 = pmax(lambda * median_sd, 0.005)) %>%
    as.data.frame()
  rownames(ds_params) <- ds_params$dataset
  dataset_names <- rownames(ds_params)

  ## Build a list of row indices grouped by dataset
  dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)

  ## Read sample names once
  samples <- rhdf5::h5read(h5file, "samples")

  ## Map CpG names to indices
  cpg_indices <- match(cpg_names_vec, cpg_names_all)
  if (anyNA(cpg_indices)) {
    stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
  }

  n_cpgs <- length(cpg_indices)
  ## Preallocate result matrix (rows = requested CpGs, cols = datasets)
  result_mat <- matrix(NA_real_, nrow = n_cpgs, ncol = length(dataset_names),
                       dimnames = list(cpg_names_vec, dataset_names))

  batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))

  for (b in seq_along(batches)) {

    row_CpGs_batch <- batches[[b]]
    message(sprintf("Loading batch %d / %d (%d CpGs) at %s",
                    b, length(batches), length(row_CpGs_batch), Sys.time()))

    M_batch <- rhdf5::h5read(
      file = h5file,
      name = "matrix",
      index = list(row_CpGs_batch, NULL),  # rows = subset of CpGs, columns = all samples
      native = TRUE ## Important for portability between programming languages!
      ## E.g. python outputs in R majors would otherwise be read in col major in R!
    )

    ## Assign dimnames
    rownames(M_batch) <- cpg_names_all[row_CpGs_batch]
    colnames(M_batch) <- samples

    ## Check sample compatibility
    if (!setequal(colnames(M_batch), metadata$sample)) {
      stop("Mismatch between HDF5 sample names and metadata$sample")
    }

    ## Reorder only if needed
    if (!identical(colnames(M_batch), metadata$sample)) {
      M_batch <- M_batch[, metadata$sample, drop = FALSE]
    }

    ## Split rows among workers (safely)
    nr <- nrow(M_batch)
    if (nr == 0) next

    NCORES_use <- max(1L, min(as.integer(NCORES), nr))
    if (NCORES_use == 1L) {
      row_chunks <- list(seq_len(nr))
    } else {
      breaks_vec <- seq(0.5, nr + 0.5, length.out = NCORES_use + 1L)
      row_chunks <- split(seq_len(nr), cut(seq_len(nr), breaks = breaks_vec,
                                           labels = FALSE, include.lowest = TRUE))
    }

    ## process chunks in parallel (or sequential if NCORES_use == 1)
    chunk_results <- parallel::mclapply(row_chunks,
      function(chunk_rows) {
        out_list <- lapply(chunk_rows, function(i) {
          vec <- M_batch[i, ]
          names(vec) <- colnames(M_batch)

          ds_present <- unique(metadata$dataset[!is.na(vec)])
          if (length(ds_present) < Nds)
            return(rep(NA_real_, length(dataset_names)))

          res_vec <- getLogPhv_oneCpG_byTissue(
            Mdf = vec,
            metadata = metadata,
            dataset_groups = dataset_groups,
            ds_params = ds_params,
            p1 = p1
          )
          ## Ensure order matches dataset_names and return numeric vector
          res_vec[dataset_names]
        })
        do.call(rbind, out_list)   # always returns a matrix
      },
      mc.cores = NCORES_use
    )

    ## chunk_results is a list of matrices (each n_rows_in_chunk x n_datasets)
    batch_combined <- do.call(rbind, chunk_results)

    ## ensure matching rownames AFTER rbind
    rownames(batch_combined) <- rownames(M_batch)

    ## write into result_mat
    pos_rows <- match(rownames(batch_combined), rownames(result_mat))
    result_mat[pos_rows, ] <- batch_combined[, dataset_names, drop = FALSE]
  }
  return(result_mat)
}

#' Run and save tissue-level hvCpG log-probabilities
#'
#' Wrapper to compute and (optionally) save the CpG x dataset matrix of log(Phv|Dk).
#' Returns the matrix invisibly.
#'
#' @param analysis Character string. Name of the analysis.
#'   If it contains `"MariasarraysREDUCED"`, a special directory structure is expected.
#' @param cpg_names_vec Character vector of CpG identifiers to process.
#' @param resultDir Character string; output directory for result files.
#' @param NCORES Integer; number of CPU cores to use in parallel processing.
#' @param p1 Numeric scalars: True positive rate.
#' @param overwrite Logical; if `TRUE`, overwrite existing result file (default FALSE).
#' @param batch_size Integer; number of CpGs to process per batch (default 10000).
#' @param dataDir Character string. Path to the directory containing input data files.
#'   Should contain `sample_metadata.tsv`, `all_medsd_lambda.tsv`, and an HDF5 matrix file.
#' @param skipsave Logical; if `TRUE`, skip saving results to disk.
#' @param Nds Integer; minimum number of datasets required per CpG (default 3).
#' @param subsetMetadata Logical or data frame. If `FALSE` (default), the full metadata
#'   is used. Otherwise, a subset of the metadata can be provided to restrict
#'   the analysis to specific samples or datasets.
#'
#' @return Invisibly returns the result matrix `(CpG x dataset)` of `log(Phv|Dk)` values.
#'   Also saves an `.RData` file in `resultDir` unless `skipsave = TRUE`.
#' @export
#'
runAndSave_tissueAnalysis <- function(
    analysis, cpg_names_vec, resultDir, NCORES, p1, overwrite = FALSE, batch_size = 1000L,
    dataDir, skipsave = FALSE, Nds = 3L, subsetMetadata = FALSE) {
  t0 <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ", round(difftime(Sys.time(), t0, units = "secs")), " seconds")

  obj_name <- paste0("results_bytissue_", analysis, "_", length(cpg_names_vec), "CpGs_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)
  resultDir <- normalizePath(resultDir, mustWork = FALSE)
  if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

  file_name <- file.path(resultDir, paste0(obj_name, ".rds"))
  if (!overwrite && file.exists(file_name)) {
    message("!! File already exists: ", file_name)
    res <- readRDS(file_name)
    return(invisible(res))
  }

  result <- getALL_LogPhv_oneCpG_byTissue_batch(
    cpg_names_vec = cpg_names_vec, NCORES = NCORES, p1 = p1,
    prep = prep, batch_size = batch_size, Nds = Nds
  )

  if (!skipsave) {
    saveRDS(result, file = file_name)
    message("Saved results to: ", file_name)
  }

  invisible(result)
}

#' Compute tissue-level hvCpG components for one CpG
#'
#' For a single CpG, computes the tissue-specific posterior probability of being
#' hypervariable, together with the numerator and denominator used in the posterior
#' calculation.
#'
#' The function returns one value per dataset and is intended as a lower-level helper
#' for batch processing and debugging.
#'
#' @param Mdf Numeric named vector of methylation values for one CpG, with sample IDs as names.
#' @param metadata Data frame with at least `sample` and `dataset` columns.
#' @param dataset_groups Named list of row indices in `metadata`, grouped by dataset.
#' @param ds_params Data frame of dataset-level parameters, including columns `sd0` and `sd1`,
#'   with row names matching dataset names.
#' @param p1 Numeric scalar giving the prior probability of hypervariability.
#'
#' @return A list with three named numeric vectors, one element per dataset:
#' \describe{
#'   \item{logPhv}{Mean log posterior probability of hypervariability.}
#'   \item{numerator}{Mean numerator term used in the posterior calculation.}
#'   \item{denominator}{Mean denominator term used in the posterior calculation.}
#' }
#' @export
#'
#' @importFrom stats dnorm
getPhv_oneCpG_byTissue_components <- function(Mdf, metadata, dataset_groups, ds_params, p1) {
  if (!is.numeric(Mdf)) stop("Mdf must be a numeric vector")
  if (is.null(names(Mdf))) stop("Mdf must have names = sample IDs")

  datasets <- unique(metadata$dataset)
  out_post <- setNames(rep(NA_real_, length(datasets)), datasets)
  out_num  <- setNames(rep(NA_real_, length(datasets)), datasets)
  out_den  <- setNames(rep(NA_real_, length(datasets)), datasets)

  for (k in datasets) {
    samp_ids <- metadata$sample[dataset_groups[[k]]]
    samp_ids <- intersect(samp_ids, names(Mdf))
    if (length(samp_ids) == 0) next

    vals <- Mdf[samp_ids]
    vals_ok <- vals[!is.na(vals)]
    if (length(vals_ok) < 3 || all(is.na(vals_ok))) next

    mu <- mean(vals_ok)
    params <- ds_params[k, ]

    p_stable <- dnorm(vals_ok, mean = mu, sd = params$sd0)
    p_var    <- dnorm(vals_ok, mean = mu, sd = params$sd1)

    denom <- p_var * p1 + p_stable * (1 - p1)
    valid <- denom > 0
    if (!any(valid)) next

    numer <- p_var[valid] * p1
    post  <- numer / denom[valid]
    post  <- pmin(pmax(post, 1e-12), 1 - 1e-12)

    out_post[k] <- mean(log(post))
    out_num[k]  <- mean(numer)
    out_den[k]  <- mean(denom[valid])
  }

  list(
    logPhv = out_post,
    numerator = out_num,
    denominator = out_den
  )
}

#' Compute tissue-level hvCpG components in batches
#'
#' Loads CpG methylation values from an HDF5 file in batches and computes, for
#' each CpG and each dataset, the log posterior probability of hypervariability
#' together with the numerator and denominator terms used in the calculation.
#'
#' @param cpg_names_vec Character vector of CpG identifiers to process.
#' @param NCORES Integer number of CPU cores to use.
#' @param p1 Numeric scalar giving the prior probability of hypervariability.
#' @param prep List returned by [prepData()], containing metadata, CpG names, HDF5 path,
#'   and dataset parameters.
#' @param batch_size Integer number of CpGs to process per batch.
#' @param Nds Integer minimum number of datasets required for a CpG.
#'
#' @return A list with components:
#' \describe{
#'   \item{logPhv}{CpG-by-dataset matrix of mean log posterior probabilities.}
#'   \item{numerator}{CpG-by-dataset matrix of mean numerator terms.}
#'   \item{denominator}{CpG-by-dataset matrix of mean denominator terms.}
#'   \item{dataset_names}{Character vector of dataset names.}
#' }
#' @export
#'
#' @importFrom rhdf5 h5read
#' @importFrom parallel mclapply
#' @importFrom stats setNames
getALL_Phvv_components_byTissue_batch <- function(cpg_names_vec, NCORES, p1, prep,
                                                  batch_size = 1000L, Nds = 3L) {
  metadata <- prep$metadata
  cpg_names_all <- prep$cpg_names_all
  h5file <- prep$h5file
  medsd_lambdas <- prep$medsd_lambdas

  ds_params <- medsd_lambdas %>%
    dplyr::select(dataset, median_sd, lambda) %>%
    dplyr::mutate(
      sd0 = pmax(median_sd, 0.005),
      sd1 = pmax(lambda * median_sd, 0.005)
    ) %>%
    as.data.frame()
  rownames(ds_params) <- ds_params$dataset
  dataset_names <- rownames(ds_params)

  dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)
  samples <- rhdf5::h5read(h5file, "samples")

  cpg_indices <- match(cpg_names_vec, cpg_names_all)
  if (anyNA(cpg_indices)) {
    stop("Some CpG names not found in HDF5: ",
         paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
  }

  n_cpgs <- length(cpg_indices)

  result_logPhv <- matrix(NA_real_, nrow = n_cpgs, ncol = length(dataset_names),
                          dimnames = list(cpg_names_vec, dataset_names))
  result_num <- result_den <- result_logPhv

  batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))

  for (b in seq_along(batches)) {
    row_CpGs_batch <- batches[[b]]
    message(sprintf("Loading batch %d / %d (%d CpGs) at %s",
                    b, length(batches), length(row_CpGs_batch), Sys.time()))

    M_batch <- rhdf5::h5read(
      file = h5file,
      name = "matrix",
      index = list(row_CpGs_batch, NULL),
      native = TRUE
    )

    rownames(M_batch) <- cpg_names_all[row_CpGs_batch]
    colnames(M_batch) <- samples

    if (!setequal(colnames(M_batch), metadata$sample)) {
      stop("Mismatch between HDF5 sample names and metadata$sample")
    }
    if (!identical(colnames(M_batch), metadata$sample)) {
      M_batch <- M_batch[, metadata$sample, drop = FALSE]
    }

    nr <- nrow(M_batch)
    if (nr == 0) next

    NCORES_use <- max(1L, min(as.integer(NCORES), nr))
    if (NCORES_use == 1L) {
      row_chunks <- list(seq_len(nr))
    } else {
      breaks_vec <- seq(0.5, nr + 0.5, length.out = NCORES_use + 1L)
      row_chunks <- split(seq_len(nr), cut(seq_len(nr), breaks = breaks_vec,
                                           labels = FALSE, include.lowest = TRUE))
    }

    chunk_results <- parallel::mclapply(
      row_chunks,
      function(chunk_rows) {
        out_list <- lapply(chunk_rows, function(i) {
          vec <- M_batch[i, ]
          names(vec) <- colnames(M_batch)

          ds_present <- unique(metadata$dataset[!is.na(vec)])
          if (length(ds_present) < Nds) {
            return(list(
              logPhv = rep(NA_real_, length(dataset_names)),
              numerator = rep(NA_real_, length(dataset_names)),
              denominator = rep(NA_real_, length(dataset_names))
            ))
          }

          res <- getPhv_oneCpG_byTissue_components(
            Mdf = vec,
            metadata = metadata,
            dataset_groups = dataset_groups,
            ds_params = ds_params,
            p1 = p1
          )

          res
        })

        list(
          logPhv = do.call(rbind, lapply(out_list, `[[`, "logPhv")),
          numerator = do.call(rbind, lapply(out_list, `[[`, "numerator")),
          denominator = do.call(rbind, lapply(out_list, `[[`, "denominator"))
        )
      },
      mc.cores = NCORES_use
    )

    batch_logPhv <- do.call(rbind, lapply(chunk_results, `[[`, "logPhv"))
    batch_num    <- do.call(rbind, lapply(chunk_results, `[[`, "numerator"))
    batch_den    <- do.call(rbind, lapply(chunk_results, `[[`, "denominator"))

    rownames(batch_logPhv) <- rownames(M_batch)
    rownames(batch_num)    <- rownames(M_batch)
    rownames(batch_den)    <- rownames(M_batch)

    pos_rows <- match(rownames(batch_logPhv), rownames(result_logPhv))
    result_logPhv[pos_rows, ] <- batch_logPhv[, dataset_names, drop = FALSE]
    result_num[pos_rows, ]    <- batch_num[, dataset_names, drop = FALSE]
    result_den[pos_rows, ]    <- batch_den[, dataset_names, drop = FALSE]
  }

  list(
    logPhv = result_logPhv,
    numerator = result_num,
    denominator = result_den,
    dataset_names = dataset_names
  )
}

#' Second wrapper to compute and (optionally) save a list of logPhv = CpG x dataset matrix of
#' log(Phv|Dk), numerator = CpG x dataset matrix of numerator, denominator = CpG x
#' dataset matrix of denominator, dataset_names = dataset names
#'
#' @param analysis Character string. Name of the analysis.
#'   If it contains `"MariasarraysREDUCED"`, a special directory structure is expected.
#' @param cpg_names_vec Character vector of CpG identifiers to process.
#' @param resultDir Character string; output directory for result files.
#' @param NCORES Integer; number of CPU cores to use in parallel processing.
#' @param p1 Numeric scalars: True positive rate.
#' @param overwrite Logical; if `TRUE`, overwrite existing result file (default FALSE).
#' @param batch_size Integer; number of CpGs to process per batch (default 10000).
#' @param dataDir Character string. Path to the directory containing input data files.
#'   Should contain `sample_metadata.tsv`, `all_medsd_lambda.tsv`, and an HDF5 matrix file.
#' @param skipsave Logical; if `TRUE`, skip saving results to disk.
#' @param Nds Integer; minimum number of datasets required per CpG (default 3).
#' @param subsetMetadata Logical or data frame. If `FALSE` (default), the full metadata
#'   is used. Otherwise, a subset of the metadata can be provided to restrict
#'   the analysis to specific samples or datasets.
#'
#' @return Invisibly returns the result list `(CpG x dataset)` of `log(Phv|Dk)` ,
#' numerator and denominator values.
#'   Also saves an `.RData` file in `resultDir` unless `skipsave = TRUE`.
#' @export
#'
runAndSave_tissueAnalysis_components <- function(
    analysis, cpg_names_vec, resultDir, NCORES, p1,
    overwrite = FALSE, batch_size = 1000L,
    dataDir, skipsave = FALSE, Nds = 3L, subsetMetadata = FALSE) {

  t0 <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ",
          round(difftime(Sys.time(), t0, units = "secs")),
          " seconds")

  obj_name <- paste0("results_bytissue_", analysis, "_",
                     length(cpg_names_vec), "CpGs_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)

  resultDir <- normalizePath(resultDir, mustWork = FALSE)
  if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

  file_name <- file.path(resultDir, paste0(obj_name, "_components.rds"))

  if (!overwrite && file.exists(file_name)) {
    message("!! File already exists: ", file_name)
    res <- readRDS(file_name)
    return(invisible(res))
  }

  result <- getALL_Phvv_components_byTissue_batch(
    cpg_names_vec = cpg_names_vec,
    NCORES = NCORES,
    p1 = p1,
    prep = prep,
    batch_size = batch_size,
    Nds = Nds
  )

  if (!skipsave) {
    saveRDS(result, file = file_name)
    message("Saved results to: ", file_name)
  }

  invisible(result)
}
