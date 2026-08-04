## hvCpG algorithm (batched HDF5 loading)
## Alice Balard

## Output columns produced per CpG (see posterior_hv_1CpG):
##   post_hv      posterior P(Z=1 | data); saturates near 0/1, good for a hard call
##   logBF        log Bayes factor (lb1 - lb0); continuous strength of evidence
##   logBF_per_ds logBF / n_ds; comparable ACROSS CpGs with different coverage (use for ranking)
##   n_hv_ds      expected number of datasets in which the CpG looks hv (0..K); most interpretable
##   n_ds         number of datasets that contributed
##   lambda_hat   crude effect size: sqrt(variance excess / baseline); "how much more variable"
RESULT_COLS <- c("post_hv", "logBF", "logBF_per_ds", "n_hv_ds", "n_ds", "lambda_hat")

#' Per-CpG hypervariability scores (hierarchical model)
#'
#' Computes a set of complementary per-CpG scores under the 3-level hierarchical
#' model. The headline posterior `post_hv` = P(Z=1 | data) saturates near 0/1
#' when many datasets/individuals point at the same yes/no question; the other
#' columns expose the underlying continuous evidence so you get a graded score.
#'
#' Returned values (a named length-6 numeric vector):
#' \describe{
#'   \item{post_hv}{Posterior P(Z=1 | data), in (0,1). Saturates; use for a hard call.}
#'   \item{logBF}{Log Bayes factor (lb1 - lb0): continuous strength of evidence for hv.}
#'   \item{logBF_per_ds}{logBF divided by the number of contributing datasets.
#'     Comparable ACROSS CpGs with different coverage — use this for ranking.}
#'   \item{n_hv_ds}{Expected number of datasets in which the CpG looks hv (0..K).
#'     Most interpretable score: "in how many datasets is it hypervariable".}
#'   \item{n_ds}{Number of datasets that contributed (>= minind individuals).}
#'   \item{lambda_hat}{Crude, unshrunken effect size = sqrt(variance excess over
#'     baseline sd0): "how much more variable than expected".}
#' }
#'
#' @param alpha0 Numeric in (0,1): global prior Pr(CpG is hv). Default 0.01.
#'   Only shifts `post_hv`; the other scores do not depend on it.
#' @inheritParams getLogLik_oneCpG_hierarchical
#' @return Named numeric vector of length 6 (see Details); all `NA` if no dataset
#'   met `minind`.
#' @importFrom stats dnorm
#' @export
posterior_hv_1CpG <- function(Mdf, metadata, dataset_groups, ds_params,
                              p0, p1, minind, alpha0 = 0.01) {
  na_out <- stats::setNames(rep(NA_real_, 6L),
                            c("post_hv","logBF","logBF_per_ds","n_hv_ds","n_ds","lambda_hat"))

  lb1 <- 0; lb0 <- 0; n_used <- 0L
  soft_k <- 0    # expected number of datasets that "look hv"
  ss <- 0        # accumulated scale-free variance excess (for lambda_hat)
  dfree <- 0     # accumulated degrees of freedom

  for (k in unique(metadata$dataset)) {
    v <- as.numeric(Mdf[, dataset_groups[[k]], drop = FALSE]); v <- v[is.finite(v)]
    if (length(v) < minind) next
    pr <- ds_params[k, ]
    if (!isTRUE(is.finite(pr$sd0) && is.finite(pr$sd1))) next  # dataset absent from medsd_lambdas
    mu <- mean(v)

    logd1 <- sum(dnorm(v, mu, pr$sd1, log = TRUE))   # log d1_k (Z_k = 1)
    logd0 <- sum(dnorm(v, mu, pr$sd0, log = TRUE))   # log d0_k (Z_k = 0)
    m <- max(logd1, logd0)

    lb1 <- lb1 + m + log(p1     * exp(logd1 - m) + (1 - p1) * exp(logd0 - m))
    lb0 <- lb0 + m + log((1-p0) * exp(logd1 - m) +  p0     * exp(logd0 - m))

    soft_k <- soft_k + 1 / (1 + exp(-(logd1 - logd0)))  # Pr(dataset k looks hv)
    ss     <- ss + sum((v - mu)^2) / pr$sd0^2           # scale-free variance excess
    dfree  <- dfree + (length(v) - 1L)

    n_used <- n_used + 1L
  }
  if (n_used == 0L) return(na_out)

  logBF <- lb1 - lb0
  lo    <- log(alpha0) - log1p(-alpha0) + logBF   # logit(post) = logit(alpha0) + logBF

  c(post_hv      = 1 / (1 + exp(-lo)),
    logBF        = logBF,
    logBF_per_ds = logBF / n_used,
    n_hv_ds      = soft_k,
    n_ds         = n_used,
    lambda_hat   = if (dfree > 0) sqrt(ss / dfree) else NA_real_)
}

#' Optimize alpha values for multiple CpGs in parallel (HDF5 batched loading)
#'
#' Loads CpG methylation data from an HDF5 matrix in batches, and computes
#' per-CpG hypervariability scores using parallel processing across cores.
#' Designed for large-scale hvCpG detection in datasets such as Atlas10X.
#'
#' @param cpg_names_vec Character vector of CpG identifiers to analyse.
#' @param NCORES Integer; number of parallel cores to use.
#' @param p0,p1 Numeric scalars: True negative and true positive rates.
#' @param prep List returned by [prepData()], containing metadata and HDF5 paths.
#' @param batch_size Integer; number of CpGs per HDF5 batch.
#' @param Nds Integer; minimum number of datasets required to compute a CpG (default 3).
#' @param minind Numeric scalar: Minimum number of individuals covered per dataset by CpG (default 3).
#' @param alpha0 Numeric in (0,1): global prior Pr(CpG is hv). Default 0.01.
#'
#' @return A numeric matrix with columns
#'   `c("post_hv","logBF","logBF_per_ds","n_hv_ds","n_ds","lambda_hat")`
#'   and rownames equal to CpG IDs. See [posterior_hv_1CpG()] for column meanings.
#'   For ranking CpGs use `logBF_per_ds` or `n_hv_ds`, not `post_hv` (which saturates).
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rhdf5 h5read
#' @importFrom stats setNames
#' @importFrom parallel mclapply
getAllOptimAlpha_parallel_batch_fast <- function(cpg_names_vec, NCORES, p0, p1, prep,
                                                 batch_size, Nds, minind, alpha0 = 0.01) {
  metadata       <- prep$metadata
  cpg_names_all  <- prep$cpg_names_all
  h5file         <- prep$h5file
  medsd_lambdas  <- prep$medsd_lambdas

  ncols <- length(RESULT_COLS)

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
  samples <- rhdf5::h5read(h5file, "samples")

  # Map CpG names to indices
  cpg_indices <- match(cpg_names_vec, cpg_names_all)
  if (anyNA(cpg_indices)) {
    stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
  }

  # Initialize result matrix (guaranteed correct shape)
  all_results <- matrix(NA_real_, nrow = length(cpg_indices), ncol = ncols,
                        dimnames = list(cpg_names_vec, RESULT_COLS))

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

    # Load block of matrix (some CpGs x all samples)
    M_batch <- rhdf5::h5read(
      file = h5file,
      name = "matrix",
      index = list(row_batches, NULL),  # rows = subset of CpGs, columns = all samples
      native = TRUE ## Important for portability between programming languages!
      ## E.g. python outputs in R majors would otherwise be read in col major in R!
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

    # SAFE sample reordering
    sample_idx <- match(metadata$sample, samples)

    if (anyNA(sample_idx)) {
      stop("Some metadata samples not found in HDF5 samples: ",
           paste(metadata$sample[is.na(sample_idx)], collapse=", "))
    }

    M_batch <- M_batch[, sample_idx, drop = FALSE]
    colnames(M_batch) <- metadata$sample

    sample_to_dataset <- metadata$dataset
    names(sample_to_dataset) <- metadata$sample

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

    # Run in parallel over chunks; each CpG returns a length-6 numeric vector
    chunk_results <- mclapply(idx_split, function(idx) {
      vapply(idx, function(i) {
        Mdf <- M_batch[i, , drop = FALSE]

        # Require at least Nds datasets with data
        datasets_present <- unique(sample_to_dataset[colnames(Mdf)[!is.na(Mdf)]])
        if (length(datasets_present) < Nds) return(rep(NA_real_, ncols))

        unname(tryCatch(
          posterior_hv_1CpG(Mdf = Mdf, metadata = metadata,
                            dataset_groups = dataset_groups, ds_params = ds_params,
                            p0 = p0, p1 = p1, minind = minind, alpha0 = alpha0),
          error = function(e) rep(NA_real_, ncols)
        ))
      }, numeric(ncols))   # vapply returns ncols x length(idx) matrix
    }, mc.cores = NCORES)

    # Each chunk is (ncols x n_cpg_in_chunk); bind columns then transpose to rows
    batch_mat <- t(do.call(cbind, chunk_results))   # (n_cpg_in_batch x ncols)

    # Store results in correct positions
    pos_in_all <- match(row_batches, cpg_indices)
    all_results[pos_in_all, ] <- batch_mat
  }

  return(all_results)
}

#' Run the hvCpG algorithm and save results to file
#'
#' Top-level driver that orchestrates hvCpG detection:
#' prepares data, runs the batched parallel scoring over CpGs,
#' and saves the resulting score matrix to disk.
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
#' @param minind Numeric scalar: Minimum number of individuals covered per dataset by CpG (default 3).
#' @param alpha0 Numeric in (0,1): global prior Pr(CpG is hv). Default 0.01.
#'
#' @return Invisibly returns the result matrix (CpG x score columns).
#'   The function also saves an `.rds` file to `resultDir` unless `skipsave = TRUE`.
#'
#' @export
runAndSave_fast <- function(
    analysis, cpg_names_vec, resultDir, NCORES, p0, p1,
    overwrite = FALSE, batch_size = 10000, dataDir,
    skipsave = FALSE, Nds = 3, subsetMetadata = FALSE, minind = 3, alpha0 = 0.01) {
  t <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ", round(Sys.time() - t), " seconds")

  obj_name <- paste0("results_", analysis, "_", length(cpg_names_vec),
                     "CpGs_", p0, "p0_", p1, "p1_", alpha0, "a0")
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
    p0 = p0, p1 = p1, prep = prep, batch_size = batch_size,
    Nds = Nds, minind = minind, alpha0 = alpha0
  )

  if (!skipsave) {
    message("Saving to file: ", file_name)
    saveRDS(result, file = file_name)
    message("Result saved successfully.")
  }

  invisible(result)
}

#' Compute CpG-level hv marginal log-likelihood (hierarchical, 3-level)
#'
#' Evaluates the marginal log-likelihood of one CpG's methylation data under a
#' three-level hierarchical model, for a given value of `alpha`:
#'
#' \itemize{
#'   \item \strong{Z} (CpG level): the CpG is hypervariable (Z=1) with prior
#'     probability `alpha`, or not (Z=0) with probability `1 - alpha`. There is a
#'     single Z shared by all datasets.
#'   \item \strong{Z_k} (dataset level): given Z, each dataset independently is in
#'     its hv state with probability `p1` (if Z=1) or `1 - p0` (if Z=0). Datasets
#'     are conditionally independent \emph{given Z}, so their evidence is
#'     multiplied only inside each Z-branch.
#'   \item \strong{M_ik} (individual level): given Z_k, each individual's value is
#'     drawn from a Normal with SD `sd0` (stable) or `sd1 = lambda * sd0` (hv).
#' }
#'
#' The CpG-level indicator Z is marginalised \strong{once}, at the top:
#' \deqn{P(M) = \alpha \prod_k [p_1 d^1_k + (1-p_1) d^0_k]
#'            + (1-\alpha) \prod_k [(1-p_0) d^1_k + p_0 d^0_k]}
#' where \eqn{d^1_k}, \eqn{d^0_k} are the products over individuals in dataset k of
#' the hv / stable Normal densities. Because `alpha` multiplies the whole product
#' over datasets (rather than entering per individual or per dataset), a large
#' dataset contributes a single confident dataset-level factor rather than one
#' vote per sample: this removes the sample-size (N) dominance of the older
#' summed-per-individual likelihood while keeping each dataset's full within-set
#' power. Note the per-CpG MLE of `alpha` is typically at the boundary (0 or 1);
#' informative in-between values come from estimating `alpha` across CpGs
#' (e.g. empirical Bayes), not from a single CpG.
#'
#'   Not used by the default pipeline (which reports scores via
#'   `posterior_hv_1CpG`); provided for empirical-Bayes estimation of the global
#'   prior `alpha0` by maximising the summed log-likelihood across CpGs.
#'
#' @param Mdf Numeric matrix (1 CpG x samples): methylation values for one CpG.
#' @param metadata Data frame with a `dataset` column (and `sample`).
#' @param dataset_groups Named list mapping dataset names to column indices of `Mdf`.
#' @param ds_params Data frame (rownames = dataset names) with columns `sd0`, `sd1`.
#' @param p0 Numeric scalar: Pr(Z_k = 0 | Z = 0), the true-negative rate at the
#'   dataset level (e.g. 0.95).
#' @param p1 Numeric scalar: Pr(Z_k = 1 | Z = 1), the true-positive rate at the
#'   dataset level (e.g. 0.65).
#' @param alpha Numeric scalar in (0, 1): prior probability that the CpG is hv
#'   (the value being optimised).
#' @param minind Integer: minimum non-missing individuals a dataset must have to
#'   contribute (datasets below this are skipped).
#'
#' @return Numeric scalar: the marginal log-likelihood \eqn{\log P(M)} for this
#'   CpG at the given `alpha`. Returns `NA_real_` if no dataset met `minind`.
#'
#' @importFrom stats dnorm
#' @export
getLogLik_oneCpG_hierarchical <- function(Mdf, metadata, dataset_groups, ds_params,
                                          p0, p1, alpha, minind) {
  datasets <- unique(metadata$dataset)
  alpha <- min(max(alpha, 1e-9), 1 - 1e-9)

  # accumulate the TWO branch log-products across datasets (log d1_k, log d0_k combined)
  logbranch1 <- 0   # sum_k log[ p1*d1_k + (1-p1)*d0_k ]   (given Z=1)
  logbranch0 <- 0   # sum_k log[ (1-p0)*d1_k + p0*d0_k ]   (given Z=0)
  n_used <- 0L

  for (k in datasets) {
    v <- as.numeric(Mdf[, dataset_groups[[k]], drop = FALSE]); v <- v[is.finite(v)]
    if (length(v) < minind) next
    pr <- ds_params[k, ]
    if (!isTRUE(is.finite(pr$sd0) && is.finite(pr$sd1))) next  # dataset absent from medsd_lambdas
    mu <- mean(v)

    logd1_k <- sum(dnorm(v, mu, pr$sd1, log = TRUE))   # log d1_k  (Z_k = 1)
    logd0_k <- sum(dnorm(v, mu, pr$sd0, log = TRUE))   # log d0_k  (Z_k = 0)

    # log[ p1*d1_k + (1-p1)*d0_k ]  via log-sum-exp
    m1 <- max(logd1_k, logd0_k)
    logbranch1 <- logbranch1 + m1 +
      log(p1 * exp(logd1_k - m1) + (1 - p1) * exp(logd0_k - m1))

    # log[ (1-p0)*d1_k + p0*d0_k ]
    logbranch0 <- logbranch0 + m1 +
      log((1 - p0) * exp(logd1_k - m1) + p0 * exp(logd0_k - m1))

    n_used <- n_used + 1L
  }
  if (n_used == 0L) return(NA_real_)

  # marginalise the CpG-level Z ONCE, at the top:  log[ alpha*B1 + (1-alpha)*B0 ]
  M <- max(logbranch1, logbranch0)
  M + log(alpha * exp(logbranch1 - M) + (1 - alpha) * exp(logbranch0 - M))
}
