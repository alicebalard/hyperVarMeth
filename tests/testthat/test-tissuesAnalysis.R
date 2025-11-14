library(testthat)
library(hyperVarMeth)

#################################
######### Perform tests #########
#################################

test_that("prepData loads mock data correctly for tissue analysis", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData(analysis = "mock", dataDir = dirname(mock$metadata))

  expect_type(prep, "list")
  expect_true(all(c("metadata", "medsd_lambdas", "cpg_names_all", "h5file") %in% names(prep)))
  expect_equal(nrow(prep$metadata), 12)  # 3 datasets × 4 samples
  expect_true(file.exists(prep$h5file))
})

test_that("getLogPhv_oneCpG_byTissue computes log-probabilities per dataset", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData("mock", dirname(mock$metadata))

  # Build dataset groups
  ds_groups <- split(seq_len(nrow(prep$metadata)), prep$metadata$dataset)

  # Minimal fake parameters
  ds_params <- data.frame(
    sd0 = rep(0.1, 3),
    sd1 = rep(0.3, 3),
    row.names = unique(prep$metadata$dataset)
  )

  # Create single CpG vector (named)
  Mdf <- matrix(runif(12), nrow = 1)
  vec <- as.numeric(Mdf[1, ])
  names(vec) <- prep$metadata$sample

  result <- getLogPhv_oneCpG_byTissue(
    Mdf = vec,
    metadata = prep$metadata,
    dataset_groups = ds_groups,
    ds_params = ds_params,
    p1 = 0.9,
    sample_n = 3
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(any(!is.na(result)))
})

test_that("getALL_LogPhv_oneCpG_byTissue_batch returns CpG x dataset matrix", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData("mock", dirname(mock$metadata))

  # Only test 10 CpGs for speed
  res <- getALL_LogPhv_oneCpG_byTissue_batch(
    cpg_names_vec = prep$cpg_names_all[1:10],
    NCORES = 1,
    p1 = 0.9,
    prep = prep,
    batch_size = 5,
    Nds = 3,
    sample_n = 3
  )

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 10)
  expect_equal(ncol(res), length(unique(prep$metadata$dataset)))
  expect_true(all(rownames(res) == prep$cpg_names_all[1:10]))
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
    batch_size = 20,
    dataDir = dirname(mock$metadata),
    skipsave = TRUE,
    Nds = 3,
    sample_n = 3
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(prep$cpg_names_all))
  expect_equal(ncol(result), length(unique(prep$metadata$dataset)))
})
