test_that("Jacobian wrong function", {
  X <- matrix(rnorm(25), 5, 5)
  d <- 4
  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))
  expect_no_error(Jacobian(X = X, a = a, d = d, p = p,
                           fun = "subdiagonal_mean_ratio_fct"))
  expect_error(Jacobian(X = X, a = a, d = d, p = p,
                        fun = "wrong_name"))
})

test_that("dvech square matrix", {
  d <- 4
  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))

  expect_no_error(
    dvech(
      X = matrix(rnorm(25), 5, 5),
      a = a,
      d = d,
      p = pu,
      inc_diag = FALSE
    )
  )

  expect_error(
    dvech(
      X = matrix(rnorm(50), 10, 5),
      a = a,
      d = d,
      p = pu,
      inc_diag = FALSE
    )
  )
})
test_that("vechp no square matrix", {
  expect_error(
    vechp(X = matrix(rnorm(10), 2, 5))
  )
})

test_that("weighted direct sum one group", {
  expect_no_error(WDirect.sumL(matrix(5, 1, 1), 2))
})

test_that("CheckBandwidth", {
  d <- 6

  expect_equal(CheckBandwidth(1, d), 1L)
  expect_equal(CheckBandwidth(d - 2, d), 4L)

  expect_error(CheckBandwidth(NA, d))
  expect_error(CheckBandwidth(c(1, 2), d))
  expect_error(CheckBandwidth(1.5, d))
  expect_error(CheckBandwidth("1", d))
  expect_error(CheckBandwidth(Inf, d))
  expect_error(CheckBandwidth(0, d))
  expect_error(CheckBandwidth(d - 1, d))
  expect_error(CheckBandwidth(1, 2))
})
