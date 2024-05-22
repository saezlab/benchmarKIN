test_that("test run_perturbBench", {
  # Create random meta and matrix
  set.seed(123)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2"),
                     target = c("E", "A"),
                     sign = c(1, -1))

  # run benchmark
  res <- run_perturbBench(act = mat, meta = meta, method_id = "test")
  expect_equal(mean(res$auroc), 1)
  expect_equal(unique(res$method), "test")
})

test_that("test prepareBench", {
  # Create random meta and matrix
  set.seed(123)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  meta <- data.frame(id = c("exp1", "exp2"),
                     target = c("A", "C"),
                     sign = c(1, -1))

  #Prepare benchmark
  out_list <- prepareBench(act = mat, meta = meta, method_id = "test")

  expect_equal(names(out_list$act), "test")
  expect_equal(length(out_list), 2)
  expect_equal(nrow(out_list$act$test), nrow(out_list$obs))
})


test_that("test scale_scores", {
  # Create random matrix
  set.seed(123)

  mat <- data.frame(exp1 = runif(5, min = -2, max = 2),
                    exp2 = runif(5, min = -3, max = 2))
  rownames(mat) <- c("A", "B", "C", "D", "E")

  # Scaling
  scaled_mat <- scale_scores(mat = mat, scaling = "sd")

  expect_equal(dim(mat), dim(scaled_mat))
})

test_that("test remove_bg", {
  # Create random network and matrix
  set.seed(123)

  mat <- data.frame(experiment = c("exp1", "exp2"),
                    A = runif(2, min = -2, max = 2),
                    B = runif(2, min = -3, max = 2),
                    C = runif(2, min = -2, max = 2),
                    D = runif(2, min = -3, max = 2))

  meta <- data.frame(sample = c("exp1", "exp2"),
                     perturb = c("A", "F"))

  # Remove experiment from background
  mat_filtered <- remove_bg(mat = mat, meta = meta)
  expect_equal(nrow(mat_filtered) < nrow(mat), T)
})
