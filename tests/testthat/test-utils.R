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


test_that("CheckRepetitions accepts positive integers", {
  expect_identical(CheckRepetitions(500), 500L)
  expect_identical(CheckRepetitions(1000L), 1000L)
})

test_that("CheckRepetitions warns for small positive integers", {
  expect_warning(
    expect_identical(CheckRepetitions(100), 100L),
    "Fewer than 500 repetitions"
  )
})

test_that("CheckRepetitions rejects invalid values", {
  expect_error(CheckRepetitions(0), "single positive integer")
  expect_error(CheckRepetitions(-1), "single positive integer")
  expect_error(CheckRepetitions(10.5), "single positive integer")
  expect_error(CheckRepetitions(NA_real_), "single positive integer")
  expect_error(CheckRepetitions(Inf), "single positive integer")
  expect_error(CheckRepetitions(c(500, 1000)), "single positive integer")
  expect_error(CheckRepetitions("1000"), "single positive integer")
  expect_error(CheckRepetitions(NULL), "single positive integer")
  expect_error(CheckRepetitions(500, NA))
  expect_error(CheckRepetitions(500, 0))
  expect_error(CheckRepetitions(500, c(1, 2)))
})

test_that("zero resampling p-values use the correct resolution", {
  expect_identical(
    format_resampling_pvalue(0, 1000),
    "p < 1e-03"
  )

  expect_identical(
    format_resampling_pvalue(0, 3000),
    "p < 3.333e-04"
  )
})

test_that("nonzero resampling p-values are printed directly", {
  expect_identical(
    format_resampling_pvalue(0.025, 1000),
    "p = 0.02500"
  )
})

test_that("p-value formatting validates its inputs", {
  expect_error(
    format_resampling_pvalue(-0.1, 1000),
    "between 0 and 1"
  )

  expect_error(
    format_resampling_pvalue(1.1, 1000),
    "between 0 and 1"
  )

})
expect_error(
  CheckRepetitions(.Machine$integer.max + 1),
  "single positive integer"
)


test_that("get_extended_matrix wrong input", {
expect_error(get_extended_matrix(V, tol = NA))
expect_error(get_extended_matrix(V, tol = -1))
expect_error(get_extended_matrix(V, tol = c(1e-8, 1e-6)))
expect_error(get_extended_matrix(V, tol = "small"))
})
