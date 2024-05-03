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
nv <- c(12,27,20,24,30,47)

# x <- TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "MC",
#                      repetitions = 1000, seed = 31415)
#
# TestCovariance_simple(X_list[[1]], hypothesis = "given-trace", A = 3)
# "equal", "equal-trace", "equal-diagonals", "given-trace", "given-matrix" and "uncorrelated"
# Testing Teststatistc

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
                                     repetitions = 1000, seed = 31415)$pvalue, 0.973)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.954)
  # Equal
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.965)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.935)
  # Equal-Diagonals
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "MC",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.936)
  expect_equal(TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-diagonals", method = "BT",
                                     repetitions = 1000, seed = 31415)$pvalue, 0.915)
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
   expect_error(TestCovariance_simple(X = lapply(X_list, t), nv = nv, hypothesis = "equal", method = "MC",
                                      repetitions = 1000, seed = NULL))
})

test_that("TestCovariance_simple single group teststatistics",{
  # Equal
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 3.5184964)
  expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 3.5184964)
  # Given Trace
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "BT",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 7.1450555))
  expect_warning(expect_equal(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "given-trace", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 7.1450555))

})

test_that("TestCovariance_simple single groups wrong hypothesis", {
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal-trace", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_error(TestCovariance_simple(X = X_list[[1]], nv = NULL, hypothesis = "equal-diagonals", method = "MC",
                                     repetitions = 1000, seed = NULL))
})

