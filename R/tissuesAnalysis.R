## hv probability comparison by tissue
## Alice Balard
## Last major update Nov 2025

#' Compute log posterior probability of variability per dataset (CpG-level)
#'
#' For a given CpG, computes the **log posterior probability** that it belongs to
#' the *variable* (hypervariable) state vs the *stable* state, separately for each dataset.
#' Within each dataset, only 3 random individuals are sampled to ensure comparability
#' across datasets with differing numbers of samples.
#'
#' @param Mdf Numeric matrix (samples x 1 CpG). Methylation values for a single CpG.
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

    datasets <- unique(metadata$dataset)

    ## Run for each dataset, and save in an list
    vector_datasets_logphv <- vector()
    for (k in datasets) {
        ## Safe extraction of the correct individual values
        sample_names <- metadata$sample[dataset_groups[[k]]]
        Mij_vals <- as.numeric(Mdf[sample_names, , drop = FALSE])

        if (length(Mij_vals) < 3 || all(is.na(Mij_vals))) next

        ## Precompute mean and SDs
        mu_jk <- mean(Mij_vals, na.rm = TRUE)

        ## Call precomputed sds for both cases
        params =  ds_params[k, ]

        ## NB: select only 3 individuals per dataset! Makes it comparable across datasets with different N
        Mij_vals_3 <- sample(Mij_vals, 3, replace = F)

        ## Vectorized density
        norm_probs <- matrix(0, nrow = length(Mij_vals_3), ncol = 2)
        norm_probs_stable <- dnorm(Mij_vals_3, mu_jk, params$sd0)
        norm_probs_variable <- dnorm(Mij_vals_3, mu_jk, params$sd1)

        ## Compute log(Phv|Dk) as the sum of logs on all individuals
        dataset_logphv <- sum(log((norm_probs_variable*p1)/
                                    (norm_probs_variable*p1 + norm_probs_stable*(1-p1))))
        vector_datasets_logphv[[k]] <- dataset_logphv
    }
    return(vector_datasets_logphv) # named vector of a logPr per dataset for this CpG
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
getALL_LogPhv_oneCpG_byTissue_batch <- function(cpg_names_vec, NCORES, p1, prep, batch_size = 1000, Nds) {
    metadata       <- prep$metadata
    cpg_names_all  <- prep$cpg_names_all
    h5file         <- prep$h5file
    medsd_lambdas  <- prep$medsd_lambdas

    ## Precompute dataset-level parameters
    ds_params = medsd_lambdas %>%
        dplyr::select(dataset, median_sd, lambda) %>%
        dplyr::mutate(sd0 = pmax(median_sd, 1e-4),
                      sd1 = pmax(lambda * median_sd, 1e-4)) %>%
        as.data.frame()
    rownames(ds_params) = ds_params$dataset

    ## Build a list of row indices grouped by dataset
    dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)

    ## Read sample names once
    samples <- h5read(h5file, "samples")

    ## Map CpG names to indices
    cpg_indices <- match(cpg_names_vec, cpg_names_all)
    if (anyNA(cpg_indices)) {
        stop("Some CpG names not found in HDF5: ", paste(cpg_names_vec[is.na(cpg_indices)], collapse = ", "))
    }

    ## Container for results
    all_results <- vector("list", length(cpg_indices))

    ## Split into batches
    batches <- split(cpg_indices, ceiling(seq_along(cpg_indices) / batch_size))

    for (b in seq_along(batches)) {
        message(sprintf(
            "Loading batch %d / %d (%d CpGs) at %s",
            b, length(batches), length(batches[[b]]), Sys.time()
        ))

     ## Load block of matrix (samples x CpGs) with rhdf5::h5read direct slice
        col_batches <- batches[[b]]
        M_batch <- h5read(
            file   = h5file,
            name   = "matrix",
            index  = list(NULL, col_batches)  # NULL = all rows, subset columns
        )

        ## Assign dimnames
        rownames(M_batch) <- samples
        colnames(M_batch) <- cpg_names_all[col_batches]

        ## Reorder rows to match metadata
        M_batch <- M_batch[metadata$sample, , drop = FALSE]

        sample_to_dataset <- setNames(metadata$dataset, metadata$sample)

        ## Split CpGs into chunks (not one per worker)
        idx_split <- split(seq_len(ncol(M_batch)), cut(seq_len(ncol(M_batch)), NCORES, labels = FALSE))

        ## Run in parallel over chunks
        chunk_results <- mclapply(idx_split, function(idx) {
            sapply(idx, function(i) {
                Mdf <- M_batch[, i, drop = FALSE]

                ## Require at least Nds datasets with data
                datasets_present <- unique(sample_to_dataset[names(Mdf[!is.na(Mdf), ])])
                if (length(datasets_present) < Nds) return(NA_real_)

                res <- tryCatch(
                    getLogPhv_oneCpG_byTissue(Mdf = Mdf, metadata = metadata, dataset_groups = dataset_groups,
                                              ds_params = ds_params, p1 = p1),
                    error = function(e) NA_real_
                )
                return(res) # named vector for one CpG
            })
        }, mc.cores = NCORES)

        ## Convert per-chunk results into a clean list of named vectors
        batch_results <- do.call(c, chunk_results)

        ## Make sure each CpG has a named vector of same length (datasets)
        dataset_names <- rownames(ds_params)
        batch_matrix <- matrix(NA_real_, nrow = length(batch_results), ncol = length(dataset_names),
                               dimnames = list(names(batch_results), dataset_names))

        ## Fill in each CpG's log-probabilities
        for (i in seq_along(batch_results)) {
            v <- batch_results[[i]]
            if (is.null(v)) next
            batch_matrix[i, names(v)] <- unlist(v)
        }

        ## Append to global results (as list of matrices)
        all_results[[b]] <- batch_matrix
    }

    ## Combine batches by rows (CpGs)
    my_matrix <- do.call(rbind, all_results)

    return(my_matrix)
}

#' Run and save tissue-level hvCpG log-probabilities
#'
#' Top-level driver that orchestrates the computation of tissue-specific
#' log posterior probabilities of hypervariability (`log(Phv|Dk)`) for each CpG,
#' saving the result as an `.RData` object.
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
    analysis, cpg_names_vec, resultDir, NCORES, p1,
    overwrite = FALSE, batch_size = 1000, dataDir,
    skipsave = FALSE, Nds = 3, subsetMetadata = FALSE
) {
  t <- Sys.time()
  prep <- prepData(analysis, dataDir, subsetMetadata)
  message("Preparing the data took ", round(Sys.time() - t), " seconds")

  obj_name <- paste0("results_bytissue_", analysis, "_",
                     length(cpg_names_vec), "CpGs_", p1, "p1")
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
  result <- getALL_LogPhv_oneCpG_byTissue_batch(
    cpg_names_vec = cpg_names_vec,
    NCORES = NCORES,
    p1 = p1,
    prep = prep,
    batch_size = batch_size,
    Nds = Nds
  )

  if (!skipsave) {
    message("Saving to file: ", file_name)
    saveRDS(result, file = file_name)
    message("Result saved successfully.")
  }

  invisible(result)
}
