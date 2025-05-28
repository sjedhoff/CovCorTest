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


##Combined Test for Covariance and Correlation
test_that("test_combined pvalues", {
  expect_equal(
    test_combined(
      X = X_list[1:2],
      nv = nv[1:2],
      seed = 31415
    )$pvalue_Variances,
    0.418
  )
  expect_equal(
    test_combined(
      X = X_list[1:2],
      nv = nv[1:2],
      seed = 31415
    )$pvalue_Correlations,
    0.016
  )
  expect_equal(
    test_combined(
      X = X_list[1:2],
      nv = nv[1:2],
      seed = 31415
    )$pvalue_Total,
    0.016
  )

})



test_that("test_combined wrong number of groups", {
  expect_error(test_combined(
    X = X_list[1],
    nv = nv[1],
    seed = 31415
  ))
  expect_error(test_combined(
    X = X_list[1:3],
    nv = nv[1:3],
    seed = 31415
  ))
})

test_that("test_combined wrong imputs", {
  expect_error(test_combined(
    X = X_list[1:2],
    nv = nv[1:2],
    seed = 31415,
    hypothesis = "equal"
  ))
  expect_error(test_combined(
    X = X_list[1:2],
    nv = nv[1:2],
    seed = 31415,
    method = "BT"
  ))
})

test_that("test_combined d=1", {
  expect_error(test_combined(
    X = list( X_list[[1]][1,], X_list[[2]][1,]),
    nv = nv[1:2],
    seed = 31415
  ))
})
