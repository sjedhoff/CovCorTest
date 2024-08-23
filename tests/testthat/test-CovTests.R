# Loading the dataset
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]
d <- 6
p <- d*(d+1)/2


X_list <- list(t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "AD",vars]),
          t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "MCI",vars]),
          t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "SCC",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "MCI",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "SCC",vars]))
X_matrix <- matrix(unlist(X_list), nrow = 6)
X <- X_list[[1]]
nv <- c(12,27,20,24,30,47)


## Multiple Groups
test_that("TestCovariance_simple multi groups Teststatistic", {
  # Equal-Trace
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 4.9047045)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 4.9047045)
  # Equal
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 2.9384292)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 2.9384292)
  # Equal-Diagonals
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 2.7304562)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 2.7304562)
})

test_that("TestCovariance_simple multi groups pvalues", {
  # Equal-Trace
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.027)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.046)
  # Equal
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.035)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.065)
  # Equal-Diagonals
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.064)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.085)
})

test_that("TestCovariance_simple multi groups wrong hypothesis", {
  expect_error(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "uncorrelated", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "given-trace", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "given-matrix", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equals", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_warning(expect_error(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = "abc")))

})

test_that("TestCovariance_simple multi groups dimensions do not fit",{
  expect_warning(expect_error(TestCovariance_simple(X = lapply(X_list, t), nv = nv, hypothesis = "equal", method = "MC",
                                      repetitions = 1000, seed = NULL)))
  expect_warning(TestCovariance_simple(X = X_list, nv = nv[-1], hypothesis = "equal", method = "MC",
                                       repetitions = 1000, seed = NULL))
})

test_that("TestCovariance_simple multi groups different input formats",{
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 4.9047045)
  expect_warning(expect_equal(TestCovariance_simple(X = X_list, nv = NULL, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 4.9047045))
  expect_equal(TestCovariance_simple(X = X_matrix, nv = nv, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 4.9047045)
  expect_error(TestCovariance_simple(X = X_matrix, nv = nv[-1], hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL))
})

## Single Group
test_that("TestCovariance_simple single group teststatistics",{
  # Equal
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 1.59396724)
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 1.59396724)
  # Given Trace
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 7.1450555))
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 7.1450555))

  # Given Matrix
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "BT",
                                                    repetitions = 1000, seed = NULL)$Teststatistic, 3.1849761))
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "MC",
                                                    repetitions = 1000, seed = NULL)$Teststatistic, 3.1849761))
})

test_that("TestCovariance_simple single group pvalue",{
  # Equal
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.212)
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.265)

  # Given Trace
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "BT",
                                                    repetitions = 1000, seed = 31415)$pvalue, 0.028))
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                                    repetitions = 1000, seed = 31415)$pvalue, 0.004))

  # Given Matrix
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "BT",
                                                    repetitions = 1000, seed = 31415)$pvalue, 0.083))
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "MC",
                                                    repetitions = 1000, seed = 31415)$pvalue, 0.032))
})

test_that("TestCovariance_simple single group given trace/matrix",{
  # Given Trace
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                     A = c(1,2,3),repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                     A = "a",repetitions = 1000, seed = NULL))
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                     A = 18, repetitions = 1000, seed = 31415)$pvalue, 0.981)

  # Given Matrix
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "MC",
                                     A = 1, repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "MC",
                                     A = matrix(rnorm(10),2,5), repetitions = 1000, seed = NULL))

  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-matrix", method = "MC",
                        A = var(t(X_list[[1]])) + seq(0,0.7,length.out = 36), repetitions = 1000, seed = 31415)$pvalue, 0.954)
})

test_that("TestCovariance_simple single group wrong hypothesis", {
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal-diagonals", method = "MC",
                                     repetitions = 1000, seed = NULL))
})


test_that("TestCovariance_simple single group different input formats",{
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = nv, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_warning(TestCovariance_simple(X = X_list[[1]], nv = 17, hypothesis = "equal", method = "MC",
                                       repetitions = 1000, seed = NULL))
})

test_that("TestCovariance_simple wrong method",{
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "abc",
                                     repetitions = 1000, seed = 31415))
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "mc",
                                     repetitions = 1000, seed = 31415)$Teststatistic,
               1.59396723756539)
})

## Structure
test_that("TestCovariance_structure teststatistics",{
  expect_equal(TestCovariance_structure(X, structure = "autoregressive", method = "MC")$Teststatistic,
                        2.13297722090024)
  expect_equal(TestCovariance_structure(X, structure = "ar", method = "BT")$Teststatistic,
                        2.13297722090024)
  expect_equal(TestCovariance_structure(X, structure = "FO-autoregressive", method = "MC")$Teststatistic,
                        1.63857996457449)
  expect_equal(TestCovariance_structure(X, structure = "FO-ar", method = "BT")$Teststatistic,
                        1.63857996457449)
  expect_equal(TestCovariance_structure(X, structure = "diagonal", method = "MC")$Teststatistic,
                        4.88780263620412)
  expect_equal(TestCovariance_structure(X, structure = "diag", method = "BT")$Teststatistic,
                        4.88780263620412)
  expect_equal(TestCovariance_structure(X, structure = "sphericity", method = "MC")$Teststatistic,
                        3.10441188918666)
  expect_equal(TestCovariance_structure(X, structure = "spher", method = "BT")$Teststatistic,
                        3.10441188918666)
  expect_equal(TestCovariance_structure(X, structure = "compoundsymmetry", method = "MC")$Teststatistic,
                        1.65692260204835)
  expect_equal(TestCovariance_structure(X, structure = "cs", method = "BT")$Teststatistic,
                        1.65692260204835)
  expect_equal(TestCovariance_structure(X, structure = "toeplitz", method = "MC")$Teststatistic,
                        1.63921322978212)
  expect_equal(TestCovariance_structure(X, structure = "toep", method = "BT")$Teststatistic,
                        1.63921322978212)
})

test_that("TestCovariance_structure teststatistics",{
  expect_equal(TestCovariance_structure(X, structure = "autoregressive", method = "MC", seed = 31415)$pvalue,
               0.103)
  expect_equal(TestCovariance_structure(X, structure = "ar", method = "BT", seed = 31415)$pvalue,
               0.158)
  expect_equal(TestCovariance_structure(X, structure = "FO-autoregressive", method = "MC", seed = 31415)$pvalue,
               0.197)
  expect_equal(TestCovariance_structure(X, structure = "FO-ar", method = "BT", seed = 31415)$pvalue,
               0.224)
  expect_equal(TestCovariance_structure(X, structure = "diagonal", method = "MC", seed = 31415)$pvalue,
               0.017)
  expect_equal(TestCovariance_structure(X, structure = "diag", method = "BT", seed = 31415)$pvalue,
               0.053)
  expect_equal(TestCovariance_structure(X, structure = "sphericity", method = "MC", seed = 31415)$pvalue,
               0.033)
  expect_equal(TestCovariance_structure(X, structure = "spher", method = "BT", seed = 31415)$pvalue,
               0.077)
  expect_equal(TestCovariance_structure(X, structure = "compoundsymmetry", method = "MC", seed = 31415)$pvalue,
               0.177)
  expect_equal(TestCovariance_structure(X, structure = "cs", method = "BT", seed = 31415)$pvalue,
               0.227)
  expect_equal(TestCovariance_structure(X, structure = "toeplitz", method = "MC", seed = 31415)$pvalue,
               0.18)
  expect_equal(TestCovariance_structure(X, structure = "toep", method = "BT", seed = 31415)$pvalue,
               0.232)
})

test_that("TestCovariance_structure wrong method",{
  expect_error(TestCovariance_structure(X = X, structure = "cs", method = "abc",
                                     repetitions = 1000, seed = 31415))
  expect_equal(TestCovariance_structure(X = X, structure = "cs", method = "mc",
                                     repetitions = 1000, seed = 31415)$Teststatistic,
               1.65692260204835)
})

test_that("TestCovariance_structure input list",{
  expect_warning(expect_equal(TestCovariance_structure(X = X_list, structure = "cs", method = "mc",
                                        repetitions = 1000, seed = 31415)$pvalue, 0.177))
  expect_equal(TestCovariance_structure(X = list(X), structure = "cs", method = "mc",
                                        repetitions = 1000, seed = 31415)$pvalue, 0.177)
})

