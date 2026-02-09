library(testthat)
library(hyperVarMeth)  # load your package namespace (Ctrl - shift - L to load_all during dvp)

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

test_that("getLogLik_oneCpG_optimized_fast computes log-likelihood", {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = 10, n_samples_per_dataset = 3)
  prep <- prepData("mock", dirname(mock$metadata))

  ds_groups <- split(seq_len(nrow(prep$metadata)), prep$metadata$dataset)
  ds_params <- data.frame(
    sd0 = rep(0.1, 10),
    sd1 = rep(0.2, 10),
    row.names = unique(prep$metadata$dataset)
  )

  Mdf <- matrix(runif(30), nrow = 1)
  loglik <- getLogLik_oneCpG_optimized_fast(Mdf, prep$metadata, ds_groups, ds_params,
                                            p0 = 0.9, p1 = 0.9, alpha = 0.5, minind = 3)

  expect_type(loglik, "double")
  expect_true(is.finite(loglik))
})

test_that("runOptim1CpG_gridrefine returns alpha between 0 and 1", {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = 10, n_samples_per_dataset = 3)
  prep <- prepData("mock", dirname(mock$metadata))

  ds_groups <- split(seq_len(nrow(prep$metadata)), prep$metadata$dataset)
  ds_params <- data.frame(
    sd0 = rep(0.1, 10),
    sd1 = rep(0.2, 10),
    row.names = unique(prep$metadata$dataset)
  )

  Mdf <- matrix(runif(30), nrow = 1)
  alpha <- runOptim1CpG_gridrefine(Mdf, prep$metadata, ds_groups, ds_params, p0 = 0.9, p1 = 0.9, minind = 3)

  expect_true(alpha >= 0 && alpha <= 1)
})

test_that("getAllOptimAlpha_parallel_batch_fast returns a matrix with correct dimensions", {
  mock <- hyperVarMeth::create_mock_hvCpG_data(
    n_datasets = 10, n_samples_per_dataset = 3)
  prep <- prepData("mock", dirname(mock$metadata))

  result <- getAllOptimAlpha_parallel_batch_fast(
    cpg_names_vec = prep$cpg_names_all[1:10],
    NCORES = 2,
    p0 = 0.9, p1 = 0.9,
    prep = prep,
    batch_size = 5,
    Nds = 3
  )

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(length(prep$cpg_names_all[1:10]), 1))
  expect_equal(colnames(result), "alpha")
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
    p0 = 0.9,
    p1 = 0.9,
    skipsave = FALSE
  )

  expect_true(file.exists(file.path(tmp_dir, grep("results_mock", list.files(tmp_dir), value = TRUE))))
})
