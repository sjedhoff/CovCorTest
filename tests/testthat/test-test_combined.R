# Loading the dataset
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]
d <- 6
p <- d * (d + 1) / 2

X_list <- list(
  EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "AD", vars],
  EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "MCI", vars],
  EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "SCC", vars],
  EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD", vars],
  EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "MCI", vars],
  EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "SCC", vars]
)

X_list_mat <- lapply(X_list, as.matrix)
X_matrix <- do.call(rbind, X_list)
X <- X_list[[1]]
X_mat <- as.matrix(X)
nv <- c(12, 27, 20, 24, 30, 47)

## Combined Test for Covariance and Correlation ---------------------------

test_that("test_combined pvalues", {
  set.seed(31415)
  res <- test_combined(
    X = X_list[1:2],
    nv = nv[1:2]
  )

  expect_equal(res$pvalue_Variances, 0.418, tolerance = 0.001)
  expect_equal(res$pvalue_Correlations, 0.016, tolerance = 0.001)
  expect_equal(res$pvalue_Total, 0.016, tolerance = 0.001)
})


test_that("test_combined returns expected object structure", {
  set.seed(31415)
  res <- test_combined(
    X = X_list[1:2],
    nv = nv[1:2]
  )

  expect_true(is.list(res))
  expect_named(
    res,
    c("pvalue_Variances", "pvalue_Correlations", "pvalue_Total",
      "Teststatistic", "repetitions", "nv"),
    ignore.order = TRUE
  )
  expect_true(is.numeric(res$pvalue_Variances))
  expect_true(is.numeric(res$pvalue_Correlations))
  expect_true(is.numeric(res$pvalue_Total))
  expect_length(res$pvalue_Variances, 1)
  expect_length(res$pvalue_Correlations, 1)
  expect_length(res$pvalue_Total, 1)
})


test_that("test_combined gives same pvalues for equivalent input formats", {
  set.seed(31415)
  res_list <- test_combined(
    X = X_list[1:2],
    nv = nv[1:2]
  )

  set.seed(31415)
  res_list_mat <- test_combined(
    X = X_list_mat[1:2],
    nv = nv[1:2]
  )

  set.seed(31415)
  res_matrix <- test_combined(
    X = do.call(rbind, X_list[1:2]),
    nv = nv[1:2]
  )

  set.seed(31415)
  expect_warning(
    res_list_inferred_nv <- test_combined(
      X = X_list[1:2],
      nv = NULL
    ),
    "will proceed with nv"
  )

  expect_equal(res_list$pvalue_Variances, res_list_mat$pvalue_Variances)
  expect_equal(res_list$pvalue_Correlations, res_list_mat$pvalue_Correlations)
  expect_equal(res_list$pvalue_Total, res_list_mat$pvalue_Total)

  expect_equal(res_list$pvalue_Variances, res_matrix$pvalue_Variances)
  expect_equal(res_list$pvalue_Correlations, res_matrix$pvalue_Correlations)
  expect_equal(res_list$pvalue_Total, res_matrix$pvalue_Total)

  expect_equal(res_list$pvalue_Variances, res_list_inferred_nv$pvalue_Variances)
  expect_equal(res_list$pvalue_Correlations, res_list_inferred_nv$pvalue_Correlations)
  expect_equal(res_list$pvalue_Total, res_list_inferred_nv$pvalue_Total)
})


test_that("test_combined wrong number of groups", {
  expect_error(
    test_combined(
      X = X_list[1],
      nv = nv[1]
    )
  )

  expect_error(
    test_combined(
      X = X_list[1:3],
      nv = nv[1:3]
    )
  )
})


test_that("test_combined wrong inputs", {
  expect_error(
    test_combined(
      X = X_list[1:2],
      nv = nv[1:2],
      hypothesis = "equal"
    )
  )

  expect_error(
    test_combined(
      X = X_list[1:2],
      nv = nv[1:2],
      method = "BT"
    )
  )
})


test_that("test_combined rejects invalid nv values", {
  expect_error(test_combined(X = X_list[1:2], nv = c(0, 27)), "positive")
  expect_error(test_combined(X = X_list[1:2], nv = c(-12, 27)), "positive")
  expect_error(test_combined(X = X_list[1:2], nv = c(12, NA)), "positive")
  expect_error(test_combined(X = X_list[1:2], nv = c(12.5, 27)), "positive")
  expect_error(test_combined(X = X_list[1:2], nv = c("12", "27")), "positive")
})


test_that("test_combined errors when nv does not match matrix rows", {
  expect_error(
    test_combined(
      X = do.call(rbind, X_list[1:2]),
      nv = c(12, 26)
    ),
    "sum of group sizes|do not align"
  )
})


test_that("test_combined warns and replaces wrong nv for list input", {
  expect_warning(
    res <- test_combined(
      X = X_list[1:2],
      nv = c(11, 28)
    ),
    "will proceed with nv"
  )

  expect_true(is.numeric(res$pvalue_Total))
})


test_that("test_combined rejects non-numeric input", {
  expect_error(
    test_combined(
      X = list(
        data.frame(a = rnorm(12), b = letters[1:12]),
        X_list[[2]]
      ),
      nv = nv[1:2]
    ),
    "numeric"
  )

  expect_error(
    test_combined(
      X = list(
        matrix(letters[1:72], nrow = 12, ncol = 6),
        as.matrix(X_list[[2]])
      ),
      nv = nv[1:2]
    ),
    "numeric"
  )
})


test_that("test_combined rejects non-matrix/data-frame list elements", {
  expect_error(
    test_combined(
      X = list(X_list[[1]], rnorm(27)),
      nv = nv[1:2]
    ),
    "matrices or data frames"
  )
})


test_that("test_combined d=1", {
  expect_warning(expect_error(
    test_combined(
      X = list(X_list[[1]][1, ], X_list[[2]][1, ]),
      nv = nv[1:2]
    )
  ))

  expect_warning(expect_error(
    test_combined(
      X = list(
        X_list[[1]][1, , drop = FALSE],
        X_list[[2]][1, , drop = FALSE]
      ),
      nv = nv[1:2]
    )
  ))
})


test_that("test_combined errors when groups have different numbers of variables", {
  expect_error(
    test_combined(
      X = list(
        matrix(rnorm(12 * 6), nrow = 12, ncol = 6),
        matrix(rnorm(27 * 5), nrow = 27, ncol = 5)
      ),
      nv = nv[1:2]
    ),
    "same number of variables|dimensions"
  )
})


test_that("test_combined errors when one group has only one subject", {
  expect_error(
    test_combined(
      X = list(
        matrix(rnorm(6), nrow = 1, ncol = 6),
        matrix(rnorm(27 * 6), nrow = 27, ncol = 6)
      ),
      nv = c(1, 27)
    ),
    "only one subject"
  )
})


test_that("CombTest object", {
  test <- CombTest()
  expect_null(print(test))
})
