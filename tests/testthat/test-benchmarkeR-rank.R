test_that("test run_rank", {
  # Create random meta and matrix
  set.seed(321)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2),
                    exp3 = runif(5, min = -2, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2", "exp3"),
                     target = c("E", "A", "E"),
                     sign = c(1, -1, 1))

  # run benchmark
  res <- run_rank(act = mat, meta = meta, average = F)
  expect_equal(ncol(res), 6)

  res_avg <- run_rank(act = mat, meta = meta, average = T)
  expect_equal(ncol(res_avg), 4)
  expect_equal(mean(res_avg$scaled_rank) < mean(res$scaled_rank), T)
})

test_that("test run_phit", {
  # Create random meta and matrix
  set.seed(321)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2),
                    exp3 = runif(5, min = -2, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2", "exp3"),
                     target = c("B", "C", "C"),
                     sign = c(1, -1, 1))

  res_1 <- run_phit(act = mat, meta = meta, k = 1)
  res_2 <- run_phit(act = mat, meta = meta, k = 2)
  res_3 <- run_phit(act = mat, meta = meta, k = 1, average = F)

  res_bool <- res_2 >= res_1
  expect_equal(res_bool, TRUE)
  expect_equal(res_1 < res_3, T)
})
