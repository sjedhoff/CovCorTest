# Loading the dataset
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]
d <- 6
p <- d * (d + 1) / 2


X_list <- list(t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "AD", vars]),
               t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "MCI", vars]),
               t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "SCC", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "AD", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "MCI", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "SCC", vars]))
X_matrix <- matrix(unlist(X_list), nrow = 6)
X <- X_list[[1]]
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
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.9047045
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "BT",
      repetitions = 1000,
      seed = NULL
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
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    2.9384292
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000,
      seed = NULL
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
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    2.7304562
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "BT",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    2.7304562
  )
})

test_that("test_covariance multi groups pvalues", {
  # Equal-Trace
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.027
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.046
  )
  # Equal
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.035
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.065
  )
  # Equal-Diagonals
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.064
  )
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-diagonals",
      method = "BT",
      repetitions = 1000,
      seed = 31415
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
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equals",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_warning(expect_error(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = "abc"
    )
  ))

})

test_that("test_covariance multi groups dimensions do not fit", {
  expect_warning(expect_error(
    test_covariance(
      X = lapply(X_list, t),
      nv = nv,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )
  ))
  expect_warning(
    test_covariance(
      X = X_list,
      nv = nv[-1],
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )
  )
})

test_that("test_covariance multi groups different input formats", {
  expect_equal(
    test_covariance(
      X = X_list,
      nv = nv,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.9047045
  )
  expect_warning(expect_equal(
    test_covariance(
      X = X_list,
      nv = NULL,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.9047045
  ))
  expect_equal(
    test_covariance(
      X = X_matrix,
      nv = nv,
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.9047045
  )
  expect_error(
    test_covariance(
      X = X_matrix,
      nv = nv[-1],
      hypothesis = "equal-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
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
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    1.59396724
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    1.59396724
  )
  # Given Trace
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "BT",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    7.1450555
  ))
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    7.1450555
  ))

  # Given Matrix
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "BT",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    3.1849761
  ))
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    3.1849761
  ))

  # Uncorrelated
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.8878026
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT",
      repetitions = 1000,
      seed = NULL
    )$Teststatistic,
    4.8878026
  )
})

test_that("test_covariance single group pvalue", {
  # Equal
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.212
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.265
  )

  # Given Trace
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.028
  ))
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.004
  ))

  # Given Matrix
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.083
  ))
  expect_warning(expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.032
  ))

  # Uncorrelated
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.017
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT",
      repetitions = 1000,
      seed = 31415
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
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      A = "a",
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-trace",
      method = "MC",
      A = 18,
      repetitions = 1000,
      seed = 31415
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
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      A = matrix(rnorm(10), 2, 5),
      repetitions = 1000,
      seed = NULL
    )
  )

  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "given-matrix",
      method = "MC",
      A = var(t(X_list[[1]])) + seq(0, 0.7, length.out = 36),
      repetitions = 1000,
      seed = 31415
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
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal-diagonals",
      method = "MC",
      repetitions = 1000,
      seed = NULL
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
      repetitions = 1000,
      seed = NULL
    )
  )
  expect_warning(
    test_covariance(
      X = X_list[[1]],
      nv = 17,
      hypothesis = "equal",
      method = "MC",
      repetitions = 1000,
      seed = NULL
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
  expect_error(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "abc",
      repetitions = 1000,
      seed = 31415
    )
  )
  expect_equal(
    test_covariance(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "equal",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$Teststatistic,
    1.59396723756539
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
})

test_that("test_covariance_structure pvalue", {
  expect_equal(
    test_covariance_structure(
      X,
      structure = "autoregressive",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.129
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "ar",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.142
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "FO-autoregressive",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.197
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "FO-ar",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.214
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "diagonal",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.017
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "diag",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.053
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "sphericity",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.033
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "spher",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.077
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "compoundsymmetry",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.177
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "cs",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.227
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "toeplitz",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.18
  )
  expect_equal(
    test_covariance_structure(
      X,
      structure = "toep",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.232
  )
})

test_that("test_covariance_structure wrong method", {
  expect_error(
    test_covariance_structure(
      X = X,
      structure = "cs",
      method = "abc",
      repetitions = 1000,
      seed = 31415
    )
  )
  expect_equal(
    test_covariance_structure(
      X = X,
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$Teststatistic,
    1.65692260204835
  )
})

test_that("test_covariance_structure input list", {
  expect_warning(expect_equal(
    test_covariance_structure(
      X = X_list,
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.177
  ))
  expect_equal(
    test_covariance_structure(
      X = list(X),
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.177
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


## Base
C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
            nrow = 1,
            ncol = 21)
Xi <- 2
test_that("test_covariance pvalue,statistic", {
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.038
  )
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )$Teststatistic,
    6.32190311190589
  )
  expect_equal(
    test_covariance(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      seed = 31415,
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
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20],
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )
  )
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20, drop = FALSE],
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      seed = 31415
    )
  )
})

test_that("test_covariance wrong method", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
              nrow = 1,
              ncol = 21)
  Xi <- 2
  expect_error(
    test_covariance(
      X = X,
      nv = NULL,
      C = C[1, 1:20],
      Xi = Xi,
      method = "abc",
      repetitions = 1000,
      seed = 31415
    )
  )
})

## Printing
test_that("print covtest", {
  X <- CovTest()
  expect_null(print(X))
})


## Missing values
test_that("missing values one group", {
  matrix <- matrix(
    c(1, NA, 3, NA, NA, NA, 2, 4, NA, NA, 6, 7, 1, 2, 3, 4, 5, 6),
    nrow = 3,
    byrow = TRUE
  )
  expect_error(expect_warning(
    test_covariance(
      X = matrix,
      nv = NULL,
      hypothesis = "uncorrelated"
    )
  ))
  matrix <- matrix(c(NA, NA, NA, NA, 1, NA, 2, 4, 4, NA, 6, 7),
                   nrow = 3,
                   byrow = TRUE)
  expect_warning(expect_warning(
    test_covariance(
      X = matrix,
      nv = NULL,
      hypothesis = "uncorrelated"
    )
  ))
})

test_that("missing values mutliple groups", {
  X <- list(
    matrix1 = matrix(
      c(1, NA, 3, NA, 2, NA, 2, 4, 1, NA, 6, 7, 1 , 2 , 3, 4, 5, 6, 7, 8, 9),
      nrow = 3,
      byrow = TRUE
    ),
    matrix2 = matrix(
      c(5, NA, 1, NA, 9, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    ),
    matrix3 = matrix(
      c(2, NA, 2, NA, 4, NA, 2, 4, 1, NA, 6, 7),
      nrow = 3,
      byrow = TRUE
    ),
    matrix4 = matrix(
      c(7, NA, 6, NA, 7, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    )
  )
  expect_warning(test_covariance(
    X = X,
    nv = c(7, 4, 4, 4),
    hypothesis = "equal"
  ))
  X <- list(
    matrix1 = matrix(
      c(1, NA, 3, NA, 2, NA, 2, 4, 1, NA, 6, 7, 1 , 2 , 3, 4, 5, 6, 7, 8, 9),
      nrow = 3,
      byrow = TRUE
    ),
    matrix2 = matrix(
      c(NA, NA, NA, NA, 9, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    ),
    matrix3 = matrix(
      c(2, NA, 2, NA, 4, NA, 2, 4, 1, NA, 6, 7),
      nrow = 3,
      byrow = TRUE
    ),
    matrix4 = matrix(
      c(7, NA, 6, NA, 7, NA, 4, NA, 3, 5, 8, NA),
      nrow = 3,
      byrow = TRUE
    )
  )
  expect_warning(expect_warning(
    test_covariance(
      X = X,
      nv = c(7, 4, 4, 4),
      hypothesis = "equal"
    )
  ))
})

## wrong dimensions: only one subject, only one variable
test_that("dim = 1 one group", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1)
  expect_error(test_covariance(X, hypothesis = "equal"))
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
  expect_error(test_covariance(X, hypothesis = "equal"))
})

test_that("dim = 1 multiple groups", {
  X <- list(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3), matrix(c(1, 2, 3), nrow = 3))
  expect_error(test_covariance(X, nv = c(2, 1), hypothesis = "equal"))
  X <- list(matrix(c(1, 2, 3, 4, 5, 6), nrow = 1), matrix(c(1, 2, 3), nrow = 1))
  expect_error(test_covariance(X, nv = c(6, 3), hypothesis = "equal"))

})
