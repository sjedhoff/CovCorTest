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



## TestCorrelation_simple
## Multiple Groups
test_that("TestCorrelation_simple multipe groups test statistics",{
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "BT")$Teststatistic, 1.93014066)
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "MC")$Teststatistic, 1.93014066)
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "TAY")$Teststatistic, 1.93014066)

})

test_that("TestCorrelation_simple multiple groups pvalue",{
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "BT",
                                      seed = 31415)$pvalue, 0.03)
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "MC",
                                      seed = 31415)$pvalue, 0.023)
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "TAY",
                                      seed = 31415)$pvalue, 0.094)
})

test_that("TestCorrelation_simple multi groups wrong hypothesis", {
  expect_error(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "uncorrelated", method = "MC"))
  expect_error(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equalcorrelated", method = "MC"))

})

test_that("TestCorrelation_simple multi groups dimensions do not fit", {
  expect_warning(expect_error(TestCorrelation_simple(X = lapply(X_list, t), nv = nv, hypothesis = "equal-correlated", method = "MC",
                        repetitions = 1000, seed = NULL)))
  expect_warning(TestCorrelation_simple(X = X_list, nv = nv[-1], hypothesis = "equal-correlated", method = "MC",
                                       repetitions = 1000, seed = NULL))
})

test_that("TestCorrelation_simple multi groups different input formats",{
  expect_equal(TestCorrelation_simple(X = X_list, nv = nv, hypothesis = "equal-correlated", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 1.93014066)
  expect_warning(expect_equal(TestCorrelation_simple(X = X_list, nv = NULL, hypothesis = "equal-correlated", method = "MC",
                                                    repetitions = 1000, seed = NULL)$Teststatistic, 1.93014066))
  expect_equal(TestCorrelation_simple(X = X_matrix, nv = nv, hypothesis = "equal-correlated", method = "MC",
                                     repetitions = 1000, seed = NULL)$Teststatistic, 1.93014066)
  expect_error(TestCorrelation_simple(X = X_matrix, nv = nv[-1], hypothesis = "equal-correlated", method = "MC",
                                     repetitions = 1000, seed = NULL))
})


## Single group
test_that("TestCorrelation_simple one group test statistics",{
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "BT")$Teststatistic, 68.906826)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "MC")$Teststatistic, 68.906826)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "TAY")$Teststatistic, 68.906826)

  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "BT")$Teststatistic, 5.2610535)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "MC")$Teststatistic, 5.2610535)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "TAY")$Teststatistic, 5.2610535)
})


test_that("TestCorrelation_simple one group pvalue",{
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "BT", seed = 31415)$pvalue, 0)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "MC", seed = 31415)$pvalue, 0)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "uncorrelated", method = "TAY")$pvalue, 0)

  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "BT", seed = 31415)$pvalue, 0.016)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "MC", seed = 31415)$pvalue, 0.006)
  expect_equal(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equal-correlated", method = "TAY", seed = 31415)$pvalue, 0.065)
})


test_that("TestCorrelation_simple one group wrong hypothesis", {
  expect_error(TestCorrelation_simple(X = X, nv = NULL, hypothesis = "equaltity", method = "BT"))
})


test_that("TestCorrelation_simple one grouo different input formats", {
  expect_error(TestCorrelation_simple(X = X_list[[1]], nv = nv, hypothesis = "uncorrelated", method = "MC",
                                     repetitions = 1000, seed = NULL))
  expect_warning(TestCorrelation_simple(X = X_list[[1]], nv = 17, hypothesis = "uncorrelated", method = "MC",
                                       repetitions = 1000, seed = NULL))
})


test_that("TestCorrelation_simple wrong method",{
  expect_error(TestCorrelation_simple(X = X_list[[1]], nv = NULL, hypothesis = "uncorrelated", method = "abc",
                                     repetitions = 1000, seed = 31415))
  expect_equal(TestCorrelation_simple(X = X_list[[1]], nv = NULL, hypothesis = "uncorrelated", method = "mc",
                                     repetitions = 1000, seed = 31415)$Teststatistic,
               68.9068261789833)
})
