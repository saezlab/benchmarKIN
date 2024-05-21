test_that("test run_zscore", {
  set.seed(123)

  # Create random network and matrix
  net <- data.frame(source = rep(c("A", "B", "C"), each = 5),
                    target = sample(rep(c("A", "B", "C", "D", "E"), each = 3)),
                    mor = 1) %>%
    dplyr::distinct()
  mat <- data.frame(exp1 = runif(5, min = -2, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  # Activity estimation
  res <- run_zscore(mat = mat, network = net, minsize = 2)

  # Test function
  expect_equal(nrow(res), 3)
  expect_equal(round(res[1,1], 2), 0.38)
})
