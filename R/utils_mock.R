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
#' @return Invisibly, a list of paths and simulated CpG names.
#'
#' @export
#' @importFrom rhdf5 h5createFile h5write
#' @importFrom utils write.table
#' @importFrom stats rnorm runif sd quantile median
create_mock_hvCpG_data <- function(
    outdir = "mock_data",
    n_datasets = 3,
    n_samples_per_dataset = 4,
    n_cpgs = 200,
    seed = 1234
) {
  set.seed(seed)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------
  # Sample metadata
  # -------------------------------
  samples <- paste0("S", seq_len(n_datasets * n_samples_per_dataset))
  datasets <- rep(paste0("DS", seq_len(n_datasets)), each = n_samples_per_dataset)
  metadata <- data.frame(sample = samples, dataset = datasets, stringsAsFactors = FALSE)
  meta_path <- file.path(outdir, "sample_metadata.tsv")
  write.table(metadata, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)

  # -------------------------------
  # Simulate CpG matrix (0–1 betas)
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

  n_draws <- 10        # number of Bernoulli draws per CpG (i.e. coverage)
  bias_prob <- c(0.1, 0.9)

  for (i in seq_len(n_datasets)) {
    sample_idx <- which(datasets == paste0("DS", i))

    # dataset-specific variability for normal CpGs
    if (i == 1) noise_sd <- runif(1, 0.02, 0.03) else noise_sd <- runif(1, 0.05, 0.25)

    for (j in sample_idx) {
      # Low CpGs: biased Bernoulli means near 0
      mat[low_idx, j] <- sapply(1:n_low, function(x) mean(rbinom(n_draws, 1, 0.1)))
      # High CpGs: biased Bernoulli means near 1
      mat[high_idx, j] <- sapply(1:n_high, function(x) mean(rbinom(n_draws, 1, 0.9)))
      # Normal CpGs: random around dataset_mean
      dataset_mean <- runif(length(normal_idx), 0, 1)
      mat[normal_idx, j] <- pmin(pmax(rnorm(length(normal_idx), mean = dataset_mean, sd = noise_sd), 0), 1)
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
    "✅ Created mock hvCpG dataset in '%s' with %d datasets, %d samples, %d CpGs.",
    normalizePath(outdir), n_datasets, length(samples), n_cpgs
  ))

  invisible(list(
    metadata = meta_path,
    medsd_lambdas = lambda_path,
    h5file = h5_path,
    cpg_names = cpg_names
  ))
}
