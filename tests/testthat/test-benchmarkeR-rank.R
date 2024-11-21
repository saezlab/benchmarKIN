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

test_that("test run_phit", {
  # Create random meta and matrix
  set.seed(321)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2"),
                     target = c("E", "A"),
                     sign = c(1, -1))

  res_1 <- run_phit(act = mat, meta = meta, k = 1)
  res_2 <- run_phit(act = mat, meta = meta, k = 2)

  res_bool <- res_2 >= res_1
  expect_equal(res_bool, TRUE)
  expect_equal(res_2, 0.5)
})
