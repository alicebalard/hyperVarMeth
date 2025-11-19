library(testthat)
library(hyperVarMeth)

#################################
######### Perform tests #########
#################################

test_that("prepData loads mock data correctly for tissue analysis", {
  mock <- hyperVarMeth::create_mock_hvCpG_data()
  prep <- prepData(analysis = "mock", dataDir = system.file("mock_data", package = "hyperVarMeth"))
  cpg_names_all  <- prep$cpg_names_all
  metadata <- prep$metadata

  samples <- rhdf5::h5read(prep$h5file, "samples")
  Mdf <- rhdf5::h5read(
    file = prep$h5file,
    name = "matrix",
    index = list(1:nrow(cpg_names_all), NULL),
    native = TRUE ## IMPORTANT
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
