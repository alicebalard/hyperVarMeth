####################################
######### Create mock data #########
####################################

metadata <- data.frame(
  sample  = paste0("S", 1:6),
  dataset = rep(c("Tissue1", "Tissue2"), each = 3)
)

# CpG matrix: 6 samples × 2 CpGs
M_matrix <- matrix(
  c(0, 0, 0, 0.5, 0.1, 0.8,
    0.05, 0.1, 0.1, 0.6, 0.85, 0.9),
  nrow = 6, ncol = 2,
  dimnames = list(metadata$sample, c("cg001", "cg002"))
)

ds_params <- data.frame(
  dataset = c("Tissue1", "Tissue2"),
  median_sd = c(0.05, 0.1),
  lambda = c(2, 2)
)
rownames(ds_params) <- ds_params$dataset
ds_params <- ds_params %>%
  dplyr::mutate(sd0 = median_sd, sd1 = lambda * median_sd)
dataset_groups <- split(seq_len(nrow(metadata)), metadata$dataset)

#################################
######### Perform tests #########
#################################

test_that("getLogPhv_oneCpG_byTissue returns named non-empty vector", {

  Mdf <- M_matrix[, "cg001", drop = FALSE]
  res <- getLogPhv_oneCpG_byTissue(
    Mdf = Mdf,
    metadata = metadata,
    dataset_groups = dataset_groups,
    ds_params = ds_params,
    p1 = 0.65
  )

  expect_true(all(names(res) %in% metadata$dataset))
  expect_false(all(is.na(res)))  # should produce real numbers
})
