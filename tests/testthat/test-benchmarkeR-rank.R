test_that("test run_rank", {
  # Create random meta and matrix
  set.seed(321)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2"),
                     target = c("E", "A"),
                     sign = c(1, -1))

  # run benchmark
  res <- run_rank(act = mat, meta = meta)
  expect_equal(ncol(res), 6)
  expect_equal(mean(res$scaled_rank), 0.5)
})
