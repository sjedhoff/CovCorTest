test_that("get_hypothesis: wrong input", {
  V <- matrix(c(1,2,3,4,5,6), nrow = 2)
  v0 <- c(1,2,3)
  expect_error(get_hypothesis(v0, V))
  expect_error(get_hypothesis(v0, v0))
  expect_error(get_hypothesis(V, V))
  v0_a <- c("a", "b", "c")
  V <- matrix(c(1,2,3,4,5,6), nrow = 3)
  expect_error(get_hypothesis(v0_a, V))
  expect_error(get_hypothesis(c(1,2,3,4), V))
  V <- matrix(c(1,2,3,4), nrow = 2)
  expect_error(get_hypothesis(v0, V))
})

test_that("get_hypothesis: right answers", {
  V <- matrix(c(1,0,1,2,0,4), nrow = 3)
  v0 <- c(1,0,1)
  expect_equal(get_hypothesis(v0, V)$hypothesis_matrix,
               matrix(c(0,1,0), nrow = 1))
  expect_equal(as.numeric(get_hypothesis(v0, V)$hypothesis_vector), 0)
})

test_that("get_hypothesis: get_extended_matrix", {
  expect_error(get_extended_matrix(matrix(c(1,2,3,4), nrow = 2)))
  V <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 3)  # 3x2 matrix
  expect_equal(get_extended_matrix(V), matrix(c(1,0,0,1,1,0,0,0,1), nrow = 3))
})
