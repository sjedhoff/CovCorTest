# Loading the dataset
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]
d <- 6
p <- d * (d + 1) / 2


X_list <- list(EEGwide[EEGwide$sex == "M" &
                         EEGwide$diagnosis == "AD", vars],
               EEGwide[EEGwide$sex == "M" &
                         EEGwide$diagnosis == "MCI", vars],
               EEGwide[EEGwide$sex == "M" &
                         EEGwide$diagnosis == "SCC", vars],
               EEGwide[EEGwide$sex == "W" &
                         EEGwide$diagnosis == "AD", vars],
               EEGwide[EEGwide$sex == "W" &
                         EEGwide$diagnosis == "MCI", vars],
               EEGwide[EEGwide$sex == "W" &
                         EEGwide$diagnosis == "SCC", vars])
X_list_mat <- lapply(X_list, as.matrix)
X_matrix <- do.call(rbind, X_list)
X <- X_list[[1]]
X_mat <- as.matrix(X)
nv <- c(12, 27, 20, 24, 30, 47)


## Multiple Groups
test_that("test_covariance multi groups Teststatistic", {
  # Equal-Trace
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000
    )$Teststatistic,
    4.9047045
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    4.9047045
  )
  # Equal
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000
    )$Teststatistic,
    2.9384292
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    2.9384292
  )
  # Equal-Diagonals
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "MC",
      repetitions = 1000
    )$Teststatistic,
    2.7304562
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    2.7304562
  )
})

test_that("test_covariance multi groups pvalues", {
  # Equal-Trace
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000
    )$pvalue,
    0.019
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0.046
  )
  # Equal
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000
    )$pvalue,
    0.042
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0.065
  )
  # Equal-Diagonals
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "MC",
      repetitions = 1000
    )$pvalue,
    0.074
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0.085
  )
})

test_that("test_covariance multi groups wrong hypothesis", {
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equals",
      method = "MC",
      repetitions = 1000
    )
  )

})

test_that("test_covariance multi groups dimensions do not fit", {
  expect_warning(expect_error(
    test_covariance(
      X = lapply(X_list, t),
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000

    )
  ))
  expect_warning(
    test_covariance(
      X = X_list,
      nv = nv[-1],
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000

    )
  )
})

test_that("test_covariance multi groups different input formats", {
  res_list <- test_covariance(
    X = X_list,
    nv = nv,
    hypothesis = "equal-trace",
    method = "MC",
    repetitions = 1000
  )

  res_list_mat <- test_covariance(
    X = X_list_mat,
    nv = nv,
    hypothesis = "equal-trace",
    method = "MC",
    repetitions = 1000
  )

  expect_warning(
    res_list_inferred_nv <- test_covariance(
      X = X_list,
      nv = NULL,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000
    ),
    "will proceed with nv"
  )

  res_matrix <- test_covariance(
    X = X_matrix,
    nv = nv,
    hypothesis = "equal-trace",
    method = "MC",
    repetitions = 1000
  )

  expect_equal(res_list$Teststatistic, 4.9047045)
  expect_equal(res_list_mat$Teststatistic, 4.9047045)
  expect_equal(res_list_inferred_nv$Teststatistic, 4.9047045)
  expect_equal(res_matrix$Teststatistic, 4.9047045)

  expect_equal(res_list$Teststatistic, res_list_mat$Teststatistic)
  expect_equal(res_list$Teststatistic, res_list_inferred_nv$Teststatistic)
  expect_equal(res_list$Teststatistic, res_matrix$Teststatistic)

  expect_error(
    test_covariance(
      X = X_matrix,
      nv = nv[-1],
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000
    )
  )
})

## Single Group
test_that("test_covariance single group teststatistics", {
  # Equal
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000
    )$Teststatistic,
    1.59396724
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    1.59396724
  )

  # Given Trace
  expect_warning(
    res_given_trace_bt <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "BT",
      repetitions = 1000
    )
  )
  expect_equal(res_given_trace_bt$Teststatistic, 7.1450555)

  expect_warning(
    res_given_trace_mc <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_equal(res_given_trace_mc$Teststatistic, 7.1450555)

  # Given Matrix
  expect_warning(
    res_given_matrix_bt <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "BT",
      repetitions = 1000
    )
  )
  expect_equal(res_given_matrix_bt$Teststatistic, 3.1849761)

  expect_warning(
    res_given_matrix_mc <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_equal(res_given_matrix_mc$Teststatistic, 3.1849761)

  # Uncorrelated
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000
    )$Teststatistic,
    4.8878026
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    4.8878026
  )
})

test_that("test_covariance single group pvalue", {
  # Equal
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000
    )$pvalue,
    0.221
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0.265
  )

  # Given Trace
  set.seed(31415)
  expect_warning(
    res_given_trace_bt <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "BT",
      repetitions = 1000
    )
  )
  expect_equal(res_given_trace_bt$pvalue, 0.028)

  set.seed(31415)
  expect_warning(
    res_given_trace_mc <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_equal(res_given_trace_mc$pvalue, 0.004)

  # Given Matrix
  set.seed(31415)
  expect_warning(
    res_given_matrix_bt <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "BT",
      repetitions = 1000
    )
  )
  expect_equal(res_given_matrix_bt$pvalue, 0.083)

  set.seed(31415)
  expect_warning(
    res_given_matrix_mc <- test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000
    )
  )
  expect_equal(res_given_matrix_mc$pvalue, 0.032)

  # Uncorrelated
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000
    )$pvalue,
    0.017
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0.047
  )
})

test_that("test_covariance single group given trace/matrix", {
  # Given Trace
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      A = c(1, 2, 3),
      repetitions = 1000

    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      A = "a",
      repetitions = 1000

    )
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      A = 18,
      repetitions = 1000

    )$pvalue,
    0.981
  )

  # Given Matrix
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      A = 1,
      repetitions = 1000

    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      A = matrix(rnorm(10), 2, 5),
      repetitions = 1000

    )
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      A = var(X_list[[1]]) + seq(0, 0.7, length.out = 36),
      repetitions = 1000

    )$pvalue,
    0.954
  )
})

test_that("test_covariance single group wrong hypothesis", {
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000

    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal-diagonals",
      method = "MC",
      repetitions = 1000

    )
  )
  expect_warning(test_covariance(
    X = X_list[[1]],
    nv = NULL,
    hypothesis = "equal",
    A = 1
  ))
  expect_error(test_covariance(
    X = X_list,
    nv = nv,
    hypothesis = "uncorrelated"
  ))
})


test_that("test_covariance single group different input formats", {
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000

    )
  )
  expect_warning(
    test_covariance(
      X = X_list[[1]],
      nv = 17,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000

    )
  )
  expect_warning(test_covariance(
    X = list(X_list[[1]]),
    nv = 15,
    hypothesis = "equal"
  ))
  expect_warning(test_covariance(
    X = list(X_list[[1]]),
    nv = nv,
    hypothesis = "equal"
  ))
})


test_that("test_covariance wrong method", {
  set.seed(31415)
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "abc",
      repetitions = 1000

    )
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "mc",
      repetitions = 1000

    )$Teststatistic,
    1.59396723756539
  )
})

test_that("test_covariance AM", {
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )
})


## Structure
test_that("test_covariance_structure teststatistics", {
  expect_equal(
    test_covariance_structure(X, structure = "autoregressive",
                              method = "MC")$Teststatistic,
    2.14534594320388
  )
  expect_equal(
    test_covariance_structure(X, structure = "ar", method = "BT")$Teststatistic,
    2.14534594320388
  )
  expect_equal(
    test_covariance_structure(X, structure = "FO-autoregressive",
                              method = "MC")$Teststatistic,
    1.63857996457449
  )
  expect_equal(
    test_covariance_structure(X, structure = "FO-ar",
                              method = "BT")$Teststatistic,
    1.63857996457449
  )
  expect_equal(
    test_covariance_structure(X, structure = "diagonal",
                              method = "MC")$Teststatistic,
    4.88780263620412
  )
  expect_equal(
    test_covariance_structure(X, structure = "diag",
                              method = "BT")$Teststatistic,
    4.88780263620412
  )
  expect_equal(
    test_covariance_structure(X, structure = "sphericity",
                              method = "MC")$Teststatistic,
    3.10441188918666
  )
  expect_equal(
    test_covariance_structure(X, structure = "spher",
                              method = "BT")$Teststatistic,
    3.10441188918666
  )
  expect_equal(
    test_covariance_structure(X, structure = "compoundsymmetry",
                              method = "MC")$Teststatistic,
    1.65692260204835
  )
  expect_equal(
    test_covariance_structure(X, structure = "cs", method = "BT")$Teststatistic,
    1.65692260204835
  )
  expect_equal(
    test_covariance_structure(X, structure = "toeplitz",
                              method = "MC")$Teststatistic,
    1.63921322978212
  )
  expect_equal(
    test_covariance_structure(X, structure = "toep",
                              method = "BT")$Teststatistic,
    1.63921322978212
  )
  expect_equal(
    test_covariance_structure(X, structure = "constant-offdiagonal",
                              method = "MC")$Teststatistic,
    1.98912114172805
  )
  expect_equal(
    test_covariance_structure(X, structure = "const-offdiag",
                              method = "BT")$Teststatistic,
    1.98912114172805
  )

  expect_equal(
    test_covariance_structure(X, structure = "standard-toeplitz",
                              method = "MC")$Teststatistic,
    2.14686664888321
  )
  expect_equal(
    test_covariance_structure(X, structure = "std-toep",
                              method = "BT")$Teststatistic,
    2.14686664888321
  )

  expect_equal(
    test_covariance_structure(X, structure = "banded",
                              bandwidth = 4,
                              method = "MC")$Teststatistic,
    5.09281161962977
  )
  expect_equal(
    test_covariance_structure(X, structure = "band",
                              bandwidth = 4,
                              method = "BT")$Teststatistic,
    5.09281161962977
  )

  expect_equal(
    test_covariance_structure(X, structure = "banded-toeplitz",
                              bandwidth = 3,
                              method = "MC")$Teststatistic,
    1.89115799397207
  )
  expect_equal(
    test_covariance_structure(X, structure = "band-toep",
                              bandwidth = 3,
                              method = "BT")$Teststatistic,
    1.89115799397207
  )

})

test_that("test_covariance_structure bandwidth", {
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded-toeplitz",
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      bandwidth = 0,
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      bandwidth = d - 1,
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      bandwidth = 1.5,
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      bandwidth = "1",
      method = "MC"
    )
  )
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "banded",
      bandwidth = Inf,
      method = "MC"
    )
  )
})

test_that("test_covariance_structure pvalue", {
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "autoregressive",
      method = "MC",

    )$pvalue,
    0.131
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "ar",
      method = "BT",

    )$pvalue,
    0.142
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "FO-autoregressive",
      method = "MC",

    )$pvalue,
    0.176
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "FO-ar",
      method = "BT",

    )$pvalue,
    0.214
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "diagonal",
      method = "MC",

    )$pvalue,
    0.017
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "diag",
      method = "BT",

    )$pvalue,
    0.053
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "sphericity",
      method = "MC",

    )$pvalue,
    0.041
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "spher",
      method = "BT",

    )$pvalue,
    0.077
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "compoundsymmetry",
      method = "MC",

    )$pvalue,
    0.175
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "cs",
      method = "BT",

    )$pvalue,
    0.227
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "toeplitz",
      method = "MC",

    )$pvalue,
    0.195
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "toep",
      method = "BT",

    )$pvalue,
    0.232
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "constant-offdiagonal",
      method = "MC"
    )$pvalue,
    0.126
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "const-offdiag",
      method = "BT"
    )$pvalue,
    0.174
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "standard-toeplitz",
      method = "MC"
    )$pvalue,
    0.109
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "std-toep",
      method = "BT"
    )$pvalue,
    0.165
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "banded",
      bandwidth = 4,
      method = "MC"
    )$pvalue,
    0.012
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "band",
      bandwidth = 4,
      method = "BT"
    )$pvalue,
    0.063
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "banded-toeplitz",
      bandwidth = 3,
      method = "MC"
    )$pvalue,
    0.173
  )

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X,
      structure = "band-toep",
      bandwidth = 3,
      method = "BT"
    )$pvalue,
    0.196
  )

})

test_that("test_covariance_structure AM", {
  expect_equal(
    test_covariance_structure(
      X,
      structure = "cs",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance_structure(
      X,
      structure = "cs",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance_structure(
      X,
      structure = "constant-offdiagonal",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance_structure(
      X,
      structure = "constant-offdiagonal",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance_structure(
      X,
      structure = "standard-toeplitz",
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance_structure(
      X,
      structure = "standard-toeplitz",
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance_structure(
      X,
      structure = "banded",
      bandwidth = 4,
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance_structure(
      X,
      structure = "banded",
      bandwidth = 4,
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )

  expect_equal(
    test_covariance_structure(
      X,
      structure = "banded-toeplitz",
      bandwidth = 3,
      AM = 0,
      method = "MC"
    )$Teststatistic,
    test_covariance_structure(
      X,
      structure = "banded-toeplitz",
      bandwidth = 3,
      AM = 1,
      method = "MC"
    )$Teststatistic,
    tolerance = 1e-10
  )
})

test_that("test_covariance_structure wrong method", {
  set.seed(31415)
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "cs",
      method = "abc",
      repetitions = 1000

    )
  )
  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X = X,
      structure = "cs",
      method = "mc",
      repetitions = 1000
    )$Teststatistic,
    1.65692260204835
  )
})

test_that("test_covariance_structure input list", {
  set.seed(31415)
  expect_warning(
    res_list <- test_covariance_structure(
      X = X_list,
      structure = "cs",
      method = "mc",
      repetitions = 1000
    )
  )
  expect_equal(res_list$pvalue, 0.175)

  set.seed(31415)
  expect_equal(
    test_covariance_structure(
      X = list(X),
      structure = "cs",
      method = "mc",
      repetitions = 1000
    )$pvalue,
    0.175
  )
})


test_that("test_covariance_structure wrong dimension / hypothesis", {
  expect_error(test_covariance_structure(
    X = X_list[[1]][1, 1:12, drop = FALSE],
    structure = "cs",
    method = "BT"
  ))
  expect_error(test_covariance_structure(
    X = X_list[[1]],
    structure = "a",
    method = "BT"
  ))
})

test_that("test_covariance_structure constant-offdiagonal dimension too small", {
  expect_error(test_covariance_structure(
  X[,1:2],
  structure = "constant-offdiagonal",
  AM = 0,
  method = "MC")
)
})

## Base
C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
            nrow = 1,
            ncol = 21)
Xi <- 2

test_that("test_covariance no hypothesis, C or Xi missing", {
  expect_error(test_covariance(
    X = X_list, nv = nv, C = C
  ))
})

test_that("test_covariance pvalue,statistic", {
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000

    )$pvalue,
    0.038
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000

    )$Teststatistic,
    6.32190311190589
  )
  set.seed(31415)
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000
      ,
      hypothesis = NULL
    )$pvalue,
    0.038
  )
})

test_that("test_covariance dimensions", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
              nrow = 1,
              ncol = 21)
  Xi <- 2
  set.seed(31415)
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20],
      Xi = Xi,
      method = "BT",
      repetitions = 1000

    )
  )
  set.seed(31415)
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20, drop = FALSE],
      Xi = Xi,
      method = "BT",
      repetitions = 1000

    )
  )
})

test_that("test_covariance wrong method", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
              nrow = 1,
              ncol = 21)
  Xi <- 2
  set.seed(31415)
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20],
      Xi = Xi,
      method = "abc",
      repetitions = 1000

    )
  )
})

## Printing
test_that("print covtest", {
  X <- CovTest()
  expect_null(print(X))
})


## Missing values
expect_missing_warning <- function(expr) {
  warnings <- character()

  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(
    any(grepl("removed", warnings)),
    info = paste(warnings, collapse = "\n")
  )

  invisible(value)
}

expect_missing_error <- function(expr) {
  warnings <- character()

  expect_error(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )

  expect_true(
    any(grepl("removed", warnings)),
    info = paste(warnings, collapse = "\n")
  )
}

test_that("missing values one group", {
  X_miss <- matrix(
    c(1, NA, 3, NA, NA, NA, 2, 4, NA, NA, 6, 7, 1, 2, 3, 4, 5, 6),
    nrow = 3,
    byrow = TRUE
  )

  expect_missing_error(
    test_covariance(
      X = X_miss,
      nv = NULL,
      hypothesis = "uncorrelated"
    )
  )

  X_miss <- matrix(
    c(NA, NA, NA, NA, 1, NA, 2, 4, 4, NA, 6, 7),
    nrow = 3,
    byrow = TRUE
  )

  res <- expect_missing_warning(
    test_covariance(
      X = X_miss,
      nv = NULL,
      hypothesis = "uncorrelated"
    )
  )

  expect_true(inherits(res, "CovTest"))
})

test_that("missing values multiple groups", {
  X_miss <- list(
    matrix1 = t(matrix(
      c(1, NA, 3, NA, 2, NA, 2, 4, 1, NA, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = 3,
      byrow = TRUE
    )),
    matrix2 = t(matrix(
      c(5, NA, 1, NA, 9, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    )),
    matrix3 = t(matrix(
      c(2, NA, 2, NA, 4, NA, 2, 4, 1, NA, 6, 7),
      nrow = 3,
      byrow = TRUE
    )),
    matrix4 = t(matrix(
      c(7, NA, 6, NA, 7, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    ))
  )

  res <- expect_missing_warning(
    test_covariance(
      X = X_miss,
      nv = c(7, 4, 4, 4),
      hypothesis = "equal"
    )
  )

  expect_true(inherits(res, "CovTest"))

  X_miss <- list(
    matrix1 = t(matrix(
      c(1, NA, 3, NA, 2, NA, 2, 4, 1, NA, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = 3,
      byrow = TRUE
    )),
    matrix2 = t(matrix(
      c(NA, NA, NA, NA, 9, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    )),
    matrix3 = t(matrix(
      c(2, NA, 2, NA, 4, NA, 2, 4, 1, NA, 6, 7),
      nrow = 3,
      byrow = TRUE
    )),
    matrix4 = t(matrix(
      c(7, NA, 6, NA, 7, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    ))
  )

  res <- expect_missing_warning(
    test_covariance(
      X = X_miss,
      nv = c(7, 4, 4, 4),
      hypothesis = "equal"
    )
  )

  expect_true(inherits(res, "CovTest"))
})

## wrong dimensions: only one subject, only one variable
test_that("dim = 1 one group", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  expect_error(test_covariance(X, hypothesis = "equal"))
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1)
  expect_error(test_covariance(X, hypothesis = "equal"))
})

test_that("dim = 1 multiple groups", {
  X <- list(t(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)), t(matrix(c(1, 2, 3), nrow = 3)))
  expect_error(test_covariance(X, nv = c(2, 1), hypothesis = "equal"))
  X <- list(t(matrix(c(1, 2, 3, 4, 5, 6), nrow = 1)), t(matrix(c(1, 2, 3), nrow = 1)))
  expect_error(test_covariance(X, nv = c(6, 3), hypothesis = "equal"))

})

## Additional regression tests for input preprocessing and validation

test_that("test_covariance returns expected object structure", {
  res <- test_covariance(
    X = X_list,
    nv = nv,
    hypothesis = "equal",
    method = "MC",
    repetitions = 1000
  )

  expect_true(is.list(res))
  expect_true(
    all(c(
      "Teststatistic",
      "pvalue",
      "C",
      "Xi",
      "AM",
      "C_original",
      "Xi_original"
    ) %in% names(res))
  )
  expect_type(res$Teststatistic, "double")
  expect_type(res$pvalue, "double")
  expect_length(res$Teststatistic, 1)
  expect_length(res$pvalue, 1)
})

test_that("test_covariance accepts lower-case method names", {
  expect_no_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "mc",
      repetitions = 1000
    )
  )

  expect_no_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "bt",
      repetitions = 1000
    )
  )
})

test_that("test_covariance rejects non-numeric data", {
  expect_error(
    test_covariance(
      X = data.frame(a = 1:5, b = letters[1:5]),
      hypothesis = "equal"
    ),
    "numeric"
  )

  expect_error(
    test_covariance(
      X = matrix(letters[1:25], 5, 5),
      hypothesis = "equal"
    ),
    "numeric"
  )

  expect_error(
    test_covariance(
      X = list(
        matrix(rnorm(25), 5, 5),
        data.frame(a = 1:5, b = letters[1:5])
      ),
      hypothesis = "equal"
    ),
    "numeric"
  )
})

test_that("test_covariance rejects invalid nv values", {
  expect_error(test_covariance(X = X_matrix, nv = 0, hypothesis = "equal"), "positive")
  expect_error(test_covariance(X = X_matrix, nv = -1, hypothesis = "equal"), "positive")
  expect_error(test_covariance(X = X_matrix, nv = c(12, NA), hypothesis = "equal"), "positive")
  expect_error(test_covariance(X = X_matrix, nv = c(12.5, 27), hypothesis = "equal"), "positive")
  expect_error(test_covariance(X = X_matrix, nv = "12", hypothesis = "equal"), "positive")
})



test_that("test_covariance works with default method and repetitions", {
  set.seed(31415)
  expect_no_error(
    res <- test_covariance(
      X = X_list[[1]],
      hypothesis = "equal"
    )
  )

  expect_true(is.numeric(res$Teststatistic))
  expect_true(is.numeric(res$pvalue))
  expect_length(res$Teststatistic, 1)
  expect_length(res$pvalue, 1)
})

test_that("test_covariance checks Xi dimensions", {
  C <- matrix(c(1, rep(0, 20)), nrow = 1)

  expect_error(
    test_covariance(
      X = X,
      C = C,
      Xi = c(1, 2),
      method = "BT",
      repetitions = 1000
    )
  )

  expect_no_error(
    test_covariance(
      X = X,
      C = C,
      Xi = 2,
      method = "BT",
      repetitions = 1000
    )
  )
})

test_that("test_covariance wrong AM", {
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "mc",
      AM = -1,
      repetitions = 1000
    )
  )

  expect_error(
    test_covariance_structure(
      X = X_list[[1]],
      structure = "diagonal",
      AM = 3,
      method = "bt",
      repetitions = 1000
    )
  )
})
