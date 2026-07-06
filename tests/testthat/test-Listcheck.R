test_that("Listcheck accepts numeric matrix without nv", {
  X <- matrix(rnorm(25), 5, 5)

  out <- expect_no_error(Listcheck(X, nv = NULL))

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
  expect_equal(dim(out$X), c(ncol(X), nrow(X)))
})

test_that("Listcheck accepts numeric data frame without nv", {
  X <- data.frame(
    a = rnorm(5),
    b = rnorm(5),
    c = rnorm(5)
  )

  out <- expect_no_error(Listcheck(X, nv = NULL))

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
  expect_equal(dim(out$X), c(ncol(X), nrow(X)))
})

test_that("Listcheck rejects non-numeric data frame", {
  X <- data.frame(
    a = rnorm(5),
    b = letters[1:5]
  )

  expect_error(
    Listcheck(X, nv = NULL),
    "numeric matrix or numeric data frame"
  )
})

test_that("Listcheck rejects non-numeric matrix", {
  X <- matrix(letters[1:9], 3, 3)

  expect_error(
    Listcheck(X, nv = NULL),
    "numeric matrix or numeric data frame"
  )
})

test_that("Listcheck rejects invalid X type", {
  X <- rnorm(10)

  expect_error(
    Listcheck(X, nv = NULL),
    "X has to be a numeric matrix, numeric data frame, or a list of those"
  )
})

test_that("Listcheck accepts list of numeric matrices", {
  X <- list(
    matrix(rnorm(20), 5, 4),
    matrix(rnorm(24), 6, 4)
  )

  expect_warning(
    out <- Listcheck(X, nv = NULL),
    "will proceed with nv"
  )

  expect_true(is.list(out$X))
  expect_equal(length(out$X), 2)
  expect_equal(out$nv, c(5, 6))
  expect_equal(vapply(out$X, nrow, integer(1)), c(4, 4))
  expect_equal(vapply(out$X, ncol, integer(1)), c(5, 6))
})

test_that("Listcheck accepts list of numeric data frames", {
  X <- list(
    data.frame(a = rnorm(4), b = rnorm(4), c = rnorm(4)),
    data.frame(a = rnorm(5), b = rnorm(5), c = rnorm(5))
  )

  expect_warning(
    out <- Listcheck(X, nv = NULL),
    "will proceed with nv"
  )

  expect_true(is.list(out$X))
  expect_equal(length(out$X), 2)
  expect_equal(out$nv, c(4, 5))
  expect_equal(vapply(out$X, nrow, integer(1)), c(3, 3))
})

test_that("Listcheck accepts mixed list of matrix and data frame", {
  X <- list(
    matrix(rnorm(12), 4, 3),
    data.frame(a = rnorm(5), b = rnorm(5), c = rnorm(5))
  )

  expect_warning(
    out <- Listcheck(X, nv = NULL),
    "will proceed with nv"
  )

  expect_true(is.list(out$X))
  expect_equal(length(out$X), 2)
  expect_equal(out$nv, c(4, 5))
  expect_equal(vapply(out$X, nrow, integer(1)), c(3, 3))
})

test_that("Listcheck rejects list with non-matrix/data-frame element", {
  X <- list(
    matrix(rnorm(9), 3, 3),
    rnorm(5)
  )

  expect_error(
    Listcheck(X, nv = NULL),
    "All list elements have to be matrices or data frames"
  )
})

test_that("Listcheck rejects list containing non-numeric data frame", {
  X <- list(
    matrix(rnorm(9), 3, 3),
    data.frame(a = rnorm(3), b = letters[1:3])
  )

  expect_error(
    Listcheck(X, nv = NULL),
    "numeric matrices or numeric data frames"
  )
})

test_that("Listcheck rejects empty list", {
  expect_error(
    Listcheck(list(), nv = NULL),
    "empty list"
  )
})

test_that("Listcheck accepts matrix with matching single nv", {
  X <- matrix(rnorm(25), 5, 5)

  expect_no_error(
    out <- Listcheck(X, nv = 5)
  )

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
  expect_equal(dim(out$X), c(5, 5))
})

test_that("Listcheck warns for matrix with non-matching single nv", {
  X <- matrix(rnorm(25), 5, 5)

  expect_warning(
    out <- Listcheck(X, nv = 4),
    "do not align"
  )

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
})

test_that("Listcheck accepts one-element list as one group", {
  X <- list(matrix(rnorm(25), 5, 5))

  expect_no_error(
    out <- Listcheck(X, nv = NULL)
  )

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
  expect_equal(dim(out$X), c(5, 5))
})

test_that("Listcheck warns for one-element list with non-matching single nv", {
  X <- list(matrix(rnorm(25), 5, 5))

  expect_warning(
    out <- Listcheck(X, nv = 4),
    "do not align"
  )

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
})

test_that("Listcheck splits matrix according to multiple nv", {
  X <- matrix(rnorm(50), 10, 5)
  nv <- c(4, 6)

  expect_no_error(
    out <- Listcheck(X, nv = nv)
  )

  expect_true(is.list(out$X))
  expect_equal(length(out$X), 2)
  expect_equal(out$nv, nv)
  expect_equal(dim(out$X[[1]]), c(5, 4))
  expect_equal(dim(out$X[[2]]), c(5, 6))
})

test_that("Listcheck errors when matrix rows do not match sum of nv", {
  X <- matrix(rnorm(50), 10, 5)
  nv <- c(4, 5)

  expect_error(
    Listcheck(X, nv = nv),
    "sum of group sizes"
  )
})

test_that("Listcheck accepts list with matching multiple nv", {
  X <- list(
    matrix(rnorm(20), 4, 5),
    matrix(rnorm(30), 6, 5)
  )

  expect_no_error(
    out <- Listcheck(X, nv = c(4, 6))
  )

  expect_true(is.list(out$X))
  expect_equal(out$nv, c(4, 6))
  expect_equal(dim(out$X[[1]]), c(5, 4))
  expect_equal(dim(out$X[[2]]), c(5, 6))
})

test_that("Listcheck warns and replaces wrong nv for multi-element list", {
  X <- list(
    matrix(rnorm(20), 4, 5),
    matrix(rnorm(30), 6, 5)
  )

  expect_warning(
    out <- Listcheck(X, nv = c(5, 5)),
    "will proceed with nv"
  )

  expect_equal(out$nv, c(4, 6))
})

test_that("Listcheck warns for one-element list with multiple nv", {
  X <- list(matrix(rnorm(25), 5, 5))

  expect_warning(
    out <- Listcheck(X, nv = c(2, 3)),
    "X contains only one group"
  )

  expect_true(is.matrix(out$X))
  expect_null(out$nv)
})

test_that("Listcheck rejects invalid nv", {
  X <- matrix(rnorm(25), 5, 5)

  expect_error(Listcheck(X, nv = 0), "positive integer")
  expect_error(Listcheck(X, nv = -1), "positive integer")
  expect_error(Listcheck(X, nv = c(2, NA)), "positive integer")
  expect_error(Listcheck(X, nv = c(2.5, 2.5)), "positive integer")
  expect_error(Listcheck(X, nv = "5"), "positive integer")
})

test_that("Listcheck removes all-NA variable rows in one group", {
  X <- matrix(rnorm(25), 5, 5)
  X[, 2] <- NA

  expect_warning(
    out <- Listcheck(X, nv = NULL),
    "row\\(s\\) with only NA values were removed"
  )

  expect_equal(nrow(out$X), 4)
  expect_false(any(apply(out$X, 1, function(row) all(is.na(row)))))
})

test_that("Listcheck removes subjects with missing values in one group", {
  X <- matrix(rnorm(25), 5, 5)
  X[1, 1] <- NA

  expect_warning(
    out <- Listcheck(X, nv = NULL),
    "subject\\(s\\) is/are removed due to missing values"
  )

  expect_equal(ncol(out$X), 4)
  expect_false(anyNA(out$X))
})

test_that("Listcheck removes all-NA variable rows across multiple groups", {
  X <- list(
    matrix(rnorm(20), 4, 5),
    matrix(rnorm(30), 6, 5)
  )

  X[[1]][, 2] <- NA

  expect_warning(
    out <- Listcheck(X, nv = c(4, 6)),
    "row\\(s\\) with only NA values were removed"
  )

  expect_equal(vapply(out$X, nrow, integer(1)), c(4, 4))
})

test_that("Listcheck removes subjects with missing values across multiple groups", {
  X <- list(
    matrix(rnorm(20), 4, 5),
    matrix(rnorm(30), 6, 5)
  )

  X[[1]][1, 1] <- NA
  X[[2]][2, 3] <- NA

  expect_warning(
    out <- Listcheck(X, nv = c(4, 6)),
    "subject\\(s\\) is/are removed due to missing values"
  )

  expect_equal(out$nv, c(3, 5))
  expect_false(any(vapply(out$X, anyNA, logical(1))))
})

test_that("Listcheck errors when groups have different numbers of variables", {
  X <- list(
    matrix(rnorm(20), 4, 5),
    matrix(rnorm(24), 6, 4)
  )

  expect_error(
    Listcheck(X, nv = c(4, 6)),
    "same number of variables"
  )
})

test_that("Listcheck errors when one group has only one subject", {
  X <- list(
    matrix(rnorm(5), 1, 5),
    matrix(rnorm(25), 5, 5)
  )

  expect_error(
    Listcheck(X, nv = c(1, 5)),
    "only one subject"
  )
})

test_that("Listcheck errors when only one variable in multiple groups", {
  X <- list(
    matrix(rnorm(4), 4, 1),
    matrix(rnorm(6), 6, 1)
  )

  expect_error(
    Listcheck(X, nv = c(4, 6)),
    "only one variable"
  )
})

test_that("Listcheck errors when only one subject in one group", {
  X <- matrix(rnorm(5), 1, 5)

  expect_error(
    Listcheck(X, nv = NULL),
    "only one subject"
  )
})

test_that("Listcheck errors when only one variable in one group", {
  X <- matrix(rnorm(5), 5, 1)

  expect_error(
    Listcheck(X, nv = NULL),
    "only one variable"
  )
})

