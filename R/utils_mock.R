#' Create a mock hvCpG dataset with realistic lambda and SD values
#'
#' Simulates a small HDF5 methylation matrix and per-dataset statistics
#' that mimic the logic used in the WGBS Atlas preprocessing pipeline:
#'  - simulate per-sample beta values with skew toward 0 and 1,
#'  - compute per-dataset median row SDs,
#'  - compute lambda = 95th percentile of row SDs / median SD.
#'
#' @param outdir Directory for output files (default "mock_data").
#' @param n_datasets Number of mock datasets (default 3).
#' @param n_samples_per_dataset Samples per dataset (default 4).
#' @param n_cpgs Number of CpGs (default 200).
#' @param seed Random seed (default 1234).
#' @param lowvar_dataset Name of dataset with reduced variability (default "DS1").
#' @param lowvar_multiplier SD scaling factor for that dataset
#' @return Invisibly, a list of paths and simulated CpG names.
#'
#' @export
#' @importFrom rhdf5 h5createFile h5write
#' @importFrom utils write.table
#' @importFrom stats rnorm rbinom runif sd quantile median
create_mock_hvCpG_data <- function(
    outdir = "inst/mock_data",
    n_datasets = 10,
    n_samples_per_dataset = 3,
    n_cpgs = 200,
    seed = 1234,
    lowvar_dataset = "DS1",
    lowvar_multiplier = 0.25
) {
  # 1. Select the true set of hypervariable CpGs
  # (10% of them)
  # 2. For each hvCpG and for each dataset except the "stable dataset",
  # decide whether it is HV (prob 0.65) or stable (prob 0.35).
  # 3. For the "fully stable dataset", force all hvCpGs to be stable.
  # 4. Use dataset-specific SD values:
  # HV CpGs: sd ~ 0.20–0.35
  # Stable CpGs: sd ~ 0.005–0.015

  set.seed(seed)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------
  # Sample metadata
  # -------------------------------
  samples <- paste0("S", seq_len(n_datasets * n_samples_per_dataset))
  datasets <- rep(paste0("DS", seq_len(n_datasets)), each = n_samples_per_dataset)
  metadata <- data.frame(sample = samples, dataset = datasets, stringsAsFactors = FALSE)
  metadata$dataset <- factor(metadata$dataset, levels = unique(metadata$dataset))
  meta_path <- file.path(outdir, "sample_metadata.tsv")
  write.table(metadata, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)

  # -------------------------------
  # Simulate CpG matrix (0-1 betas)
  # -------------------------------
  cpg_names <- paste0("cg", sprintf("%06d", seq_len(n_cpgs)))
  mat <- matrix(NA_real_, nrow = n_cpgs, ncol = length(samples),
                dimnames = list(cpg_names, samples))

  # Assign some CpGs to be near 0 or 1 in all samples using biased Bernoulli means
  n_low <- round(0.1 * n_cpgs)   # 10% near 0
  n_high <- round(0.1 * n_cpgs)  # 10% near 1
  low_idx <- 1:n_low
  high_idx <- (n_cpgs - n_high + 1):n_cpgs
  normal_idx <- setdiff(seq_len(n_cpgs), c(low_idx, high_idx))

  n_draws <- 10        # number of Bernoulli draws per CpG (i.e. coverage 10X)
  bias_prob <- c(0.1, 0.9)

  # -------------------------------
  # 1. Choose true hvCpGs (10%)
  # -------------------------------
  prop_hv <- 0.10
  n_hv <- round(n_cpgs * prop_hv)

  hv_cpg_idx <- sample(seq_len(n_cpgs), n_hv)

  # The rest are stable everywhere
  stable_cpg_idx <- setdiff(seq_len(n_cpgs), hv_cpg_idx)

  # -------------------------------
  # 2. Assign hv/stable per dataset
  # -------------------------------
  p_hv <- 0.65
  dataset_names <- paste0("DS", seq_len(n_datasets))

  # Define the dataset where hvCpGs must always be stable
  stable_dataset <- "DS1"

  hv_mask <- matrix(
    FALSE,
    nrow = n_cpgs,
    ncol = n_datasets,
    dimnames = list(cpg_names, dataset_names)
  )

  # For true hvCpGs: HV in 65% datasets except DS1
  for (ds in dataset_names) {
    if (ds == stable_dataset) {
      hv_mask[hv_cpg_idx, ds] <- FALSE  # force stable
    } else {
      hv_mask[hv_cpg_idx, ds] <- as.logical(
        rbinom(n_hv, 1, p_hv)
      )
    }
  }

  # All non-hvCpGs remain FALSE (stable everywhere)
  for (i in seq_len(n_datasets)) {
    ds <- paste0("DS", i)
    sample_idx <- which(datasets == ds)

    # CpG-specific SDs (determined by HV mask)
    sd_vec <- ifelse(
      hv_mask[, ds],
      runif(n_cpgs, 0.20, 0.35),
      runif(n_cpgs, 0.005, 0.015)
    )

    # Add one dataset with low variability
    sd_vec <- sd_vec * (if (ds == lowvar_dataset) lowvar_multiplier else 1)

    # Dataset-specific mean for normal CpGs (fixed across all samples in DS)
    dataset_mean <- runif(length(normal_idx), 0, 1)

    for (j in sample_idx) {

      # Low CpGs
      mat[low_idx, j] <-
        sapply(seq_len(n_low), function(x) mean(rbinom(n_draws, 1, 0.1)))

      # High CpGs
      mat[high_idx, j] <-
        sapply(seq_len(n_high), function(x) mean(rbinom(n_draws, 1, 0.9)))

      # Normal CpGs
      mat[normal_idx, j] <-
        pmin(pmax(
          rnorm(
            length(normal_idx),
            mean = dataset_mean,
            sd   = sd_vec[normal_idx]
          ),
          0), 1)
    }
  }

  # Randomly drop a few values (simulate low coverage)
  drop_mask <- matrix(runif(length(mat)) < 0.05, nrow = nrow(mat))
  mat[drop_mask] <- NA_real_

  # -------------------------------
  # Compute per-dataset medians & lambdas
  # -------------------------------
  medsd_lambdas <- lapply(unique(datasets), function(ds) {
    submat <- mat[, datasets == ds, drop = FALSE]
    if (ncol(submat) >= 3) {
      row_sds <- apply(submat, 1, sd, na.rm = TRUE)
      median_sd <- median(row_sds, na.rm = TRUE)
      p95 <- quantile(row_sds, 0.95, na.rm = TRUE)
      lambda <- as.numeric(p95 / median_sd)
    } else {
      median_sd <- NA_real_
      lambda <- NA_real_
    }
    data.frame(dataset = ds, median_sd = median_sd, lambda = lambda)
  })
  medsd_lambdas <- do.call(rbind, medsd_lambdas)
  lambda_path <- file.path(outdir, "all_medsd_lambda.tsv")
  write.table(medsd_lambdas, lambda_path, sep = "\t", quote = FALSE, row.names = FALSE)

  # -------------------------------
  # Save matrix as HDF5
  # -------------------------------
  h5_path <- file.path(outdir, "all_matrix_noscale.h5")
  if (file.exists(h5_path)) file.remove(h5_path)
  rhdf5::h5createFile(h5_path)
  rhdf5::h5write(mat, h5_path, "matrix")
  rhdf5::h5write(samples, h5_path, "samples")
  rhdf5::h5write(cpg_names, h5_path, "cpg_names")
  rhdf5::h5write(datasets, h5_path, "sample_groups")

  message(sprintf(
    "Created mock hvCpG dataset in '%s' with %d datasets, %d samples, %d CpGs.",
    normalizePath(outdir), n_datasets, length(samples), n_cpgs
  ))

  invisible(list(
    metadata = meta_path,
    medsd_lambdas = lambda_path,
    h5file = h5_path,
    cpg_names = cpg_names
  ))
}
