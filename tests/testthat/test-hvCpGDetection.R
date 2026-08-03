library(testthat)
library(hyperVarMeth)  # load your package namespace (Ctrl-shift-L to load_all during dvp)

#################################
######### Perform tests #########
#################################

test_that("prepData loads mock data correctly", {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = 10, n_samples_per_dataset = 3)
  prep <- prepData(analysis = "mock", dataDir = dirname(mock$metadata))

  expect_type(prep, "list")
  expect_true(all(c("metadata", "medsd_lambdas", "cpg_names_all", "h5file") %in% names(prep)))
  expect_true(file.exists(prep$h5file))
})

# helper: shared mock + params, avoids repetition
.mock_setup <- function(n_datasets = 10) {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = n_datasets, n_samples_per_dataset = 3)
  prep <- prepData("mock", dirname(mock$metadata))
  ds_groups <- split(seq_len(nrow(prep$metadata)), prep$metadata$dataset)
  ds_params <- data.frame(
    sd0 = rep(0.1, n_datasets),
    sd1 = rep(0.2, n_datasets),
    row.names = unique(prep$metadata$dataset)
  )
  list(prep = prep, ds_groups = ds_groups, ds_params = ds_params)
}

test_that("getLogLik_oneCpG_hierarchical returns a finite marginal log-likelihood", {
  s <- .mock_setup()
  Mdf <- matrix(runif(30), nrow = 1)
  ll <- getLogLik_oneCpG_hierarchical(
    Mdf, s$prep$metadata, s$ds_groups, s$ds_params,
    p0 = 0.95, p1 = 0.65, alpha = 0.5, minind = 3)

  expect_type(ll, "double")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
})

test_that("getLogLik_oneCpG_hierarchical returns NA when no dataset meets minind", {
  s <- .mock_setup()
  Mdf <- matrix(runif(30), nrow = 1)
  ll <- getLogLik_oneCpG_hierarchical(
    Mdf, s$prep$metadata, s$ds_groups, s$ds_params,
    p0 = 0.95, p1 = 0.65, alpha = 0.5, minind = 999)   # impossible to meet
  expect_true(is.na(ll))
})

test_that("posterior_hv_1CpG returns a probability in (0,1)", {
  s <- .mock_setup()
  Mdf <- matrix(runif(30), nrow = 1)
  post <- posterior_hv_1CpG(
    Mdf, s$prep$metadata, s$ds_groups, s$ds_params,
    p0 = 0.95, p1 = 0.65, minind = 3, alpha0 = 0.05)

  expect_type(post, "double")
  expect_length(post, 1)
  expect_true(is.finite(post))
  expect_gte(post, 0)
  expect_lte(post, 1)
})

test_that("posterior_hv_1CpG: hv-looking data scores higher than stable-looking data", {
  s <- .mock_setup()
  n <- nrow(s$prep$metadata)
  # tight values -> should look stable (low posterior)
  tight <- matrix(rnorm(n, 0, 0.02), nrow = 1)
  # wide values -> should look hv (higher posterior)
  wide  <- matrix(rnorm(n, 0, 0.30), nrow = 1)

  p_tight <- posterior_hv_1CpG(tight, s$prep$metadata, s$ds_groups, s$ds_params,
                               p0 = 0.95, p1 = 0.65, minind = 3, alpha0 = 0.05)
  p_wide  <- posterior_hv_1CpG(wide,  s$prep$metadata, s$ds_groups, s$ds_params,
                               p0 = 0.95, p1 = 0.65, minind = 3, alpha0 = 0.05)

  expect_gt(p_wide, p_tight)
})

test_that("posterior_hv_1CpG increases with alpha0 (monotonic in the prior)", {
  s <- .mock_setup()
  Mdf <- matrix(runif(30), nrow = 1)
  p_lo <- posterior_hv_1CpG(Mdf, s$prep$metadata, s$ds_groups, s$ds_params,
                            p0 = 0.95, p1 = 0.65, minind = 3, alpha0 = 0.01)
  p_hi <- posterior_hv_1CpG(Mdf, s$prep$metadata, s$ds_groups, s$ds_params,
                            p0 = 0.95, p1 = 0.65, minind = 3, alpha0 = 0.20)
  expect_gte(p_hi, p_lo)
})

test_that("getAllOptimAlpha_parallel_batch_fast returns a matrix with correct dimensions", {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = 10, n_samples_per_dataset = 3)
  prep <- prepData("mock", dirname(mock$metadata))

  result <- getAllOptimAlpha_parallel_batch_fast(
    cpg_names_vec = prep$cpg_names_all[1:10],
    NCORES = 2,
    p0 = 0.95, p1 = 0.65,
    prep = prep,
    batch_size = 5,
    Nds = 3,
    minind = 3,
    alpha0 = 0.05
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(10, 1))
  expect_equal(colnames(result), "post_hv")           # was "alpha"
  # values are posteriors: NA or within [0,1]
  vals <- result[, 1]
  expect_true(all(is.na(vals) | (vals >= 0 & vals <= 1)))
})

test_that("runAndSave_fast runs and saves results", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData("mock", dirname(mock$metadata))

  tmp_dir <- tempdir()
  result <- runAndSave_fast(
    dataDir = system.file("mock_data", package = "hyperVarMeth"),
    analysis = "mock",
    cpg_names_vec = prep$cpg_names_all,
    resultDir = tmp_dir,
    NCORES = 1,
    p0 = 0.95,
    p1 = 0.65,
    alpha0 = 0.05,
    skipsave = FALSE
  )

  saved <- grep("results_mock", list.files(tmp_dir), value = TRUE)
  expect_true(length(saved) >= 1)
  expect_true(file.exists(file.path(tmp_dir, saved[1])))
})
