library(testthat)
library(hyperVarMeth)

#################################
######### Perform tests #########
#################################

test_that("prepData loads mock data correctly for tissue analysis", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData(analysis = "mock", dataDir = system.file("mock_data", package = "hyperVarMeth"))
  cpg_names_all <- prep$cpg_names_all
  metadata <- prep$metadata

  samples <- rhdf5::h5read(prep$h5file, "samples")
  Mdf <- rhdf5::h5read(
    file = prep$h5file,
    name = "matrix",
    index = list(1:nrow(cpg_names_all), NULL),
    native = TRUE
  )

  expect_type(prep, "list")
  expect_true(all(c("metadata", "medsd_lambdas", "cpg_names_all", "h5file") %in% names(prep)))
  expect_true(class(metadata) == "data.frame")
  expect_true(file.exists(prep$h5file))
})

test_that("runAndSave_tissueAnalysis runs and optionally saves", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData("mock", dirname(mock$metadata))

  tmp_dir <- tempdir()

  result <- runAndSave_tissueAnalysis(
    analysis = "mock",
    cpg_names_vec = prep$cpg_names_all,
    resultDir = tmp_dir,
    NCORES = 1,
    p1 = 0.9,
    overwrite = TRUE,
    batch_size = 100,
    dataDir = system.file("mock_data", package = "hyperVarMeth"),
    skipsave = TRUE
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(prep$cpg_names_all))
  expect_equal(ncol(result), length(unique(prep$metadata$dataset)))
})

test_that("getPhv_oneCpG_byTissue_components returns posterior and components", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData(analysis = "mock", dataDir = system.file("mock_data", package = "hyperVarMeth"))
  metadata <- prep$metadata
  dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)
  medsd_lambdas <- prep$medsd_lambdas

  ds_params <- medsd_lambdas[, c("dataset", "median_sd", "lambda")]
  ds_params$sd0 <- pmax(ds_params$median_sd, 0.005)
  ds_params$sd1 <- pmax(ds_params$lambda * ds_params$median_sd, 0.005)
  rownames(ds_params) <- ds_params$dataset

  cpg <- prep$cpg_names_all[1]
  samples <- rhdf5::h5read(prep$h5file, "samples")
  M_batch <- rhdf5::h5read(
    file = prep$h5file,
    name = "matrix",
    index = list(1, NULL),
    native = TRUE
  )
  names(M_batch) <- samples

  res <- getPhv_oneCpG_byTissue_components(
    Mdf = M_batch,
    metadata = metadata,
    dataset_groups = dataset_groups,
    ds_params = ds_params,
    p1 = 0.9
  )

  expect_type(res, "list")
  expect_true(all(c("logPhv", "numerator", "denominator") %in% names(res)))
  expect_equal(length(res$logPhv), length(unique(metadata$dataset)))
  expect_equal(length(res$numerator), length(unique(metadata$dataset)))
  expect_equal(length(res$denominator), length(unique(metadata$dataset)))
})

test_that("getALL_Phvv_components_byTissue_batch returns list of matrices", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData(analysis = "mock", dataDir = system.file("mock_data", package = "hyperVarMeth"))

  tmp_dir <- tempdir()

  result <- getALL_Phvv_components_byTissue_batch(
    cpg_names_vec = prep$cpg_names_all[1:5],
    NCORES = 1,
    p1 = 0.9,
    prep = prep,
    batch_size = 5,
    Nds = 3
  )

  expect_type(result, "list")
  expect_true(all(c("logPhv", "numerator", "denominator", "dataset_names") %in% names(result)))
  expect_true(is.matrix(result$logPhv))
  expect_true(is.matrix(result$numerator))
  expect_true(is.matrix(result$denominator))
  expect_equal(nrow(result$logPhv), 5)
  expect_equal(ncol(result$logPhv), length(unique(prep$metadata$dataset)))
})

test_that("runAndSave_tissueAnalysis_components runs and saves component output", {
  tmp_dir <- tempdir()
  prep <- prepData("mock", system.file("mock_data", package = "hyperVarMeth"))

  result <- runAndSave_tissueAnalysis_components(
    analysis = "mock",
    cpg_names_vec = prep$cpg_names_all[1:5],
    resultDir = tmp_dir,
    NCORES = 1,
    p1 = 0.9,
    overwrite = TRUE,
    batch_size = 5,
    dataDir = system.file("mock_data", package = "hyperVarMeth"),
    skipsave = TRUE
  )

  expect_type(result, "list")
  expect_true(all(c("logPhv", "numerator", "denominator", "dataset_names") %in% names(result)))
  expect_true(is.matrix(result$logPhv))
  expect_true(is.matrix(result$numerator))
  expect_true(is.matrix(result$denominator))
})
