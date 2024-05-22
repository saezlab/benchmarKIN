test_that("test load_perturbData", {
  mat <- load_perturbData()
  expect_equal(dim(mat)[1], 65686)
  expect_equal(dim(mat)[2], 225)
  expect_equal(is.data.frame(mat), T)
})

test_that("test load_meta", {
  meta <- load_meta()
  expect_equal(is.character(meta$id), T)
  expect_equal(is.data.frame(meta), T)
})
