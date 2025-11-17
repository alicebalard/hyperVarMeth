## hvCpG algorithm (batched HDF5 loading)
## Alice Balard
## Last major update Nov 2025

#' Compute log-likelihood for one CpG across all datasets
#'
#' Calculates the overall log-likelihood of observing methylation values
#' for a given CpG across multiple datasets, given mixture parameters and
#' hypervariable probabilities.
#'
#' @param Mdf Numeric matrix (1 CpG x samples). Methylation values for a single CpG.
#' @param metadata Data frame containing `sample` and `dataset` columns.
#' @param dataset_groups Named list mapping dataset names to row indices in `metadata`.
#' @param ds_params Data frame with precomputed parameters per dataset:
#'   columns `sd0`, `sd1`, and rownames corresponding to dataset names.
#' @param p0,p1 Numeric scalars: True negative and true positive rates.
#' @param alpha Numeric scalar (0-1): Probability of being a hypervariable site.
#'
#' @return Numeric scalar: the summed log-likelihood for the CpG across datasets.
#'
#' @export
#'
#' @importFrom stats dnorm
getLogLik_oneCpG_optimized_fast <- function(Mdf, metadata, dataset_groups, ds_params, p0, p1, alpha) {

  # samples <- metadata$sample # to rm? useless?
  datasets <- unique(metadata$dataset)

  # Precompute p0/p1 matrix and mixture probs
  p0p1_mat <- matrix(c(p0, 1 - p1, 1 - p0, p1), nrow = 2, byrow = TRUE)
  proba_hvCpG_vec <- c(1 - alpha, alpha)

  log_P_Mj <- 0

  for (k in datasets) {
    Mij_vals <- as.numeric(Mdf[, dataset_groups[[k]], drop = FALSE])
    if (length(Mij_vals) < 3 || all(is.na(Mij_vals))) next

    # Precompute mean and SDs
    mu_jk <- mean(Mij_vals, na.rm = TRUE)

    ## Cal precomputed sds for both cases
    params =  ds_params[k, ]

    # Vectorized density
    norm_probs <- matrix(0, nrow = length(Mij_vals), ncol = 2)
    norm_probs[,1] <- dnorm(Mij_vals, mu_jk, params$sd0)
    norm_probs[,2] <- dnorm(Mij_vals, mu_jk, params$sd1)

    # Compute zjk_probs safely
    zjk_probs <- array(0, dim = c(length(Mij_vals), 2, 2))
    for (zjk in 0:1) {
      zjk_probs[,, zjk+1] <- cbind(
        norm_probs[, zjk+1] * p0p1_mat[zjk+1,1],
        norm_probs[, zjk+1] * p0p1_mat[zjk+1,2]
      )
    }

    # Sum across latent states and mixture
    col_sums <- apply(zjk_probs, c(1,3), sum)
    dataset_loglik <- sum(log(rowSums(col_sums %*% proba_hvCpG_vec)))
    if (!is.finite(dataset_loglik)) dataset_loglik <- 0

    log_P_Mj <- log_P_Mj + dataset_loglik
  }

  return(log_P_Mj)
}

#' Estimate the optimal alpha value for one CpG
#'
#' Performs a coarse grid search followed by local refinement using
#' Brent optimization to find the alpha value (probability of hypervariability)
#' that maximises the log-likelihood for a single CpG.
#'
#' @param Mdf Numeric matrix (1 CpG x samples). Methylation values for a single CpG.
#' @param metadata Data frame containing `sample` and `dataset` columns.
#' @param dataset_groups Named list mapping dataset names to row indices in `metadata`.
#' @param ds_params Data frame with precomputed parameters per dataset:
#'   columns `sd0`, `sd1`, and rownames corresponding to dataset names.
#' @param p0,p1 Numeric scalars: True negative and true positive rates.
#'
#' @return Numeric scalar: estimated optimal alpha for the CpG.
#' @export
#'
#' @importFrom stats optim
runOptim1CpG_gridrefine <- function(Mdf, metadata, dataset_groups, ds_params, p0, p1) {
  # Step 1. Coarse grid search (0, 0.05, 0.1...)
  grid <- seq(0, 1, length.out = 21)
  logliks <- vapply(grid, function(a) {
    getLogLik_oneCpG_optimized_fast(Mdf, metadata, dataset_groups, ds_params, p0, p1, a)
  }, numeric(1))

  best_idx <- which.max(logliks)
  alpha_start <- grid[best_idx]

  # Step 2. Local refinement with Brent, only in neighborhood
  lower <- ifelse(best_idx == 1, 0, grid[best_idx - 1])
  upper <- ifelse(best_idx == length(grid), 1, grid[best_idx + 1])

  resOpt <- optim(
    par = alpha_start,
    fn = function(alpha) {
      getLogLik_oneCpG_optimized_fast(Mdf, metadata, dataset_groups, ds_params, p0, p1, alpha)
    },
    method = "Brent",
    lower = lower, upper = upper,
    control = list(fnscale = -1)
  )
  return(resOpt$par)
}

#' Optimize alpha values for multiple CpGs in parallel (HDF5 batched loading)
#'
#' Loads CpG methylation data from an HDF5 matrix in batches, and computes
#' the optimal alpha value for each CpG using parallel processing across cores.
#' Designed for large-scale hvCpG detection in datasets such as Atlas10X.
#'
#' @param cpg_names_vec Character vector of CpG identifiers to analyse.
#' @param NCORES Integer; number of parallel cores to use.
#' @param p0,p1 Numeric scalars: True negative and true positive rates.
#' @param prep List returned by [prepData()], containing metadata and HDF5 paths.
#' @param batch_size Integer; number of CpGs per HDF5 batch.
#' @param Nds Integer; minimum number of datasets required to compute a CpG (default 3).
#'
#' @return A numeric matrix with one column (`alpha`) and rownames equal to CpG IDs.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rhdf5 h5read
#' @importFrom stats setNames
#' @importFrom parallel mclapply
getAllOptimAlpha_parallel_batch_fast <- function(cpg_names_vec, NCORES, p0, p1, prep, batch_size, Nds) {
  metadata       <- prep$metadata
  cpg_names_all  <- prep$cpg_names_all
  h5file         <- prep$h5file
  medsd_lambdas  <- prep$medsd_lambdas

  ## Precompute dataset-level parameters
  ds_params <- medsd_lambdas %>%
    dplyr::select(dataset, median_sd, lambda) %>%
    dplyr::mutate(sd0 = pmax(median_sd, 1e-4),
                  sd1 = pmax(lambda * median_sd, 1e-4)) %>%
    as.data.frame()
  rownames(ds_params) <- ds_params$dataset

  ## Build a list of row indices grouped by dataset
  dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)

  # Read sample names once
  samples <- h5read(h5file, "samples")

  # Map CpG names to indices
  cpg_indices <- match(cpg_names_vec, cpg_names_all)
  if (anyNA(cpg_indices)) {
    stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
  }

  # Initialize result vector (guaranteed correct length)
  all_results_vec <- rep(NA_real_, length(cpg_indices))

  # Split into batches
  batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))

  for (b in seq_along(batches)) {
    message(sprintf(
      "Loading batch %d / %d (%d CpGs) at %s",
      b, length(batches), length(batches[[b]]), Sys.time()
    ))

    row_batches <- batches[[b]]

    # Skip empty batch
    if (length(row_batches) == 0) next

    # Load block of matrix (samples x CpGs)
    M_batch <- h5read(
      file = h5file,
      name = "matrix",
      index = list(row_batches, NULL)  # rows = subset of CpGs, columns = all samples
    )

    # Force to matrix safely
    M_batch <- as.matrix(M_batch)  # converts vector, matrix, or array into proper 2D matrix

    if (nrow(M_batch) != length(row_batches) || ncol(M_batch) != length(samples)) {
      stop(sprintf("Matrix shape mismatch: expected %d x %d, got %d x %d",
                   length(row_batches), length(samples), nrow(M_batch), ncol(M_batch)))
    }

    # Assign dimnames
    rownames(M_batch) <- cpg_names_all[row_batches]
    colnames(M_batch) <- samples

    # Reorder samples to match metadata
    M_batch <- M_batch[, metadata$sample, drop = FALSE]

    sample_to_dataset <- setNames(metadata$dataset, metadata$sample)

    # Split CpGs into chunks (not one per worker) -- SAFELY
    nrows <- nrow(M_batch)
    if (is.null(nrows) || nrows == 0) {
      # nothing to process in this batch
      next
    }

    if (is.na(NCORES) || NCORES < 1) {
      message("!!Invalid NCORES (", NCORES, ") - defaulting to 1.")
      NCORES <- 1
    }

    if (is.null(nrows) || nrows == 0) {
      next  # nothing to process
    }

    if (nrows == 1) {
      idx_split <- list(1L)
    } else {
      NCORES_use <- min(as.integer(NCORES), as.integer(nrows))
      if (is.na(NCORES_use) || NCORES_use < 1) {
        message("!! Invalid NCORES_use (", NCORES_use, ") - forcing to 1.")
        NCORES_use <- 1L
      }

      if (nrows <= 1L || NCORES_use <= 1L) {
        # Single-row or single-core: one group only
        idx_split <- list(seq_len(nrows))
      } else {
        # Safely compute breaks so cut() always has valid range
        breaks_vec <- seq(0.5, nrows + 0.5, length.out = NCORES_use + 1L)
        idx_split <- split(
          seq_len(nrows),
          cut(seq_len(nrows), breaks = breaks_vec, labels = FALSE, include.lowest = TRUE)
        )
      }
    }

    # Run in parallel over chunks
    chunk_results <- mclapply(idx_split, function(idx) {
      sapply(idx, function(i) {
        Mdf <- M_batch[i, , drop = FALSE]

        # Require at least Nds datasets with data
        datasets_present <- unique(sample_to_dataset[colnames(Mdf)[!is.na(Mdf)]])
        if (length(datasets_present) < Nds) return(NA_real_)

        tryCatch(
          runOptim1CpG_gridrefine(Mdf = Mdf, metadata = metadata,
                                  dataset_groups = dataset_groups,
                                  ds_params = ds_params, p0 = p0, p1 = p1),
          error = function(e) NA_real_
        )
      })
    }, mc.cores = NCORES)

    batch_results <- unlist(chunk_results, use.names = FALSE)

    # Store results in correct positions
    pos_in_all <- match(row_batches, cpg_indices)
    all_results_vec[pos_in_all] <- batch_results
  }

  # Build final matrix
  my_matrix <- matrix(all_results_vec, ncol = 1)
  rownames(my_matrix) <- cpg_names_vec
  colnames(my_matrix) <- "alpha"

  return(my_matrix)
}

#' Run the hvCpG algorithm and save results to file
#'
#' Top-level driver that orchestrates hvCpG detection:
#' prepares data, runs the batched parallel optimisation over CpGs,
#' and saves the resulting alpha estimates to disk.
#'
#' @param analysis Character string. Name of the analysis.
#'   If it contains `"MariasarraysREDUCED"`, a special directory structure is expected.
#' @param cpg_names_vec Character vector of CpG identifiers to process.
#' @param resultDir Character string; output directory for result files.
#' @param NCORES Integer; number of CPU cores to use in parallel processing.
#' @param p0,p1 Numeric scalars: True negative and true positive rates.
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
#' @return Invisibly returns the result matrix (CpG x alpha).
#'   The function also saves an `.RData` file to `resultDir` unless `skipsave = TRUE`.
#'
#' @export
runAndSave_fast <- function(
    analysis, cpg_names_vec, resultDir, NCORES, p0, p1,
    overwrite = FALSE, batch_size = 10000, dataDir,
    skipsave = FALSE, Nds = 3, subsetMetadata = FALSE
) {
  t <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ", round(Sys.time() - t), " seconds")

  obj_name <- paste0("results_", analysis, "_", length(cpg_names_vec),
                     "CpGs_", p0, "p0_", p1, "p1")
  obj_name <- gsub("[^[:alnum:]_]", "_", obj_name)

  resultDir <- normalizePath(resultDir, mustWork = FALSE)
  if (!dir.exists(resultDir)) {
    dir.create(resultDir, recursive = TRUE)
    message("New result directory ", resultDir, " created")
  }

  file_name <- file.path(resultDir, paste0(obj_name, ".rds"))
  if (!overwrite && file.exists(file_name)) {
    message("!! File already exists: ", file_name)
    return(invisible(NULL))
  }

  # Run batch + parallel processing
  result <- getAllOptimAlpha_parallel_batch_fast(
    cpg_names_vec = cpg_names_vec, NCORES = NCORES,
    p0 = p0, p1 = p1, prep = prep, batch_size = batch_size, Nds = Nds
  )

  if (!skipsave) {
    message("Saving to file: ", file_name)
    saveRDS(result, file = file_name)
    message("Result saved successfully.")
  }

  invisible(result)
}

