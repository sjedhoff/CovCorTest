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



## test_correlation
## Multiple Groups
test_that("test_correlation multipe groups test statistics", {
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "BT"
    )$Teststatistic,
    1.93014066
  )
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "MC"
    )$Teststatistic,
    1.93014066
  )
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "TAY"
    )$Teststatistic,
    1.93014066
  )

})

test_that("test_correlation multiple groups pvalue", {
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "BT"
    )$pvalue,
    0.03
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "MC"
    )$pvalue,
    0.023
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "TAY"
    )$pvalue,
    0.094
  )
})

test_that("test_correlation multi groups wrong hypothesis", {
  expect_error(test_correlation(
    X = X_list,
    nv = nv,
    hypothesis = "uncorrelated",
    method = "MC"
  ))
  expect_error(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equalcorrelated",
      method = "MC"
    )
  )

})

test_that("test_correlation multi groups dimensions do not fit", {
  expect_warning(expect_error(
    test_correlation(
      X = lapply(X_list, t),
      nv = nv,
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )
  ))
  expect_warning(
    test_correlation(
      X = X_list,
      nv = nv[-1],
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )
  )
})

test_that("test_correlation multi groups different input formats", {
  expect_equal(
    test_correlation(
      X = X_list,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )$Teststatistic,
    1.93014066
  )
  expect_warning(expect_equal(
    test_correlation(
      X = X_list,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )$Teststatistic,
    1.93014066
  ))
  expect_equal(
    test_correlation(
      X = X_matrix,
      nv = nv,
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )$Teststatistic,
    1.93014066
  )
  expect_error(
    test_correlation(
      X = X_matrix,
      nv = nv[-1],
      hypothesis = "equal-correlated",
      method = "MC",
      repetitions = 1000

    )
  )
})


## Single group
test_that("test_correlation one group test statistics", {
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT"
    )$Teststatistic,
    68.906826
  )
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC"
    )$Teststatistic,
    68.906826
  )
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "TAY"
    )$Teststatistic,
    68.906826
  )

  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "BT"
    )$Teststatistic,
    5.2610535
  )
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "MC"
    )$Teststatistic,
    5.2610535
  )
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "TAY"
    )$Teststatistic,
    5.2610535
  )
})


test_that("test_correlation one group pvalue", {
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "BT"
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "MC"
    )$pvalue,
    0
  )
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "TAY"
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "BT"
    )$pvalue,
    0.016
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "MC"
    )$pvalue,
    0.006
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      hypothesis = "equal-correlated",
      method = "TAY"
    )$pvalue,
    0.065
  )
})


test_that("test_correlation one group wrong hypothesis", {
  expect_error(test_correlation(
    X = X,
    nv = NULL,
    hypothesis = "equaltity",
    method = "BT"
  ))
})


test_that("test_correlation one group different input formats", {
  expect_error(
    test_correlation(
      X = X_list[[1]],
      nv = nv,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000

    )
  )
  expect_warning(
    test_correlation(
      X = X_list[[1]],
      nv = 17,
      hypothesis = "uncorrelated",
      method = "MC",
      repetitions = 1000

    )
  )
})


test_that("test_correlation wrong method", {
  set.seed(31415)
  expect_error(
    test_correlation(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "abc",
      repetitions = 1000
    )
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X_list[[1]],
      nv = NULL,
      hypothesis = "uncorrelated",
      method = "mc",
      repetitions = 1000
    )$Teststatistic,
    68.9068261789833
  )
})

test_that("test_correlation d=1", {
  expect_error(test_correlation(
    X = X_list[[1]][1, 1:12, drop = FALSE],
    nv = NULL,
    hypothesis = "uncorrelated"
  ))
})

test_that("test_correlation no hypothesis, C or Xi missing", {
  expect_error(test_correlation(
    X = X_list, nv = nv, C = c(1,2,3)
  ))
})

## Structure
test_that("test_correlation_structure teststatic", {
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "BT"
    )$Teststatistic,
    3.95677381790978
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "TAY"
    )$Teststatistic,
    3.95677381790978
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "MC"
    )$Teststatistic,
    3.95677381790978
  )

  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "BT"
    )$Teststatistic,
    68.9068261789833
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "TAY"
    )$Teststatistic,
    68.9068261789833
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "MC"
    )$Teststatistic,
    68.9068261789833
  )

  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "BT"
    )$Teststatistic,
    5.26105347405293
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "TAY"
    )$Teststatistic,
    5.26105347405293
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "MC"
    )$Teststatistic,
    5.26105347405293
  )

  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "BT"
    )$Teststatistic,
    4.88304727273834
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "TAY"
    )$Teststatistic,
    4.88304727273834
  )
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "MC"
    )$Teststatistic,
    4.88304727273834
  )

})

test_that("test_correlation_structure pvalue", {
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "BT"
    )$pvalue,
    0.021
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "TAY"
    )$pvalue,
    0.104
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "MC"
    )$pvalue,
    0.02
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[2]],
      structure = "diag",
      method = "BT"
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "TAY"
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "MC"
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "BT"
    )$pvalue,
    0.007
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "TAY"
    )$pvalue,
    0.065
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "MC"
    )$pvalue,
    0.006
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "BT"
    )$pvalue,
    0.009
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "TAY"
    )$pvalue,
    0.124
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "MC"
    )$pvalue,
    0.006
  )

})


test_that("test_correlation_structure wrong method/hypothesis", {
  set.seed(31415)
  expect_error(
    test_correlation_structure(
      X = X,
      structure = "hcs",
      method = "abc",
      repetitions = 1000
    )
  )
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = X,
      structure = "hcs",
      method = "mc",
      repetitions = 1000
    )$Teststatistic,
    5.26105347405293
  )
  expect_error(test_correlation_structure(X = X, structure = "a"))
})

test_that("test_correlation_structure input list", {
  set.seed(31415)
  expect_warning(expect_equal(
    test_correlation_structure(
      X = X_list,
      structure = "hcs",
      method = "mc",
      repetitions = 1000
    )$pvalue,
    0.006
  ))
  set.seed(31415)
  expect_equal(
    test_correlation_structure(
      X = list(X),
      structure = "hcs",
      method = "mc",
      repetitions = 1000
    )$pvalue,
    0.006
  )
})

test_that("test_correlation_structure d=1", {
expect_error(test_correlation_structure(X = X_list[[1]][1, 1:12, drop = FALSE],
                                         structure = "Har"))
})

## Base
test_that("test_correlation pvalue,statistic", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
              nrow = 1,
              ncol = 15)
  Xi <- 2
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000
    )$pvalue,
    0
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000
    )$Teststatistic,
    51.212638310026
  )
  set.seed(31415)
  expect_equal(
    test_correlation(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "BT",
      repetitions = 1000,
      hypothesis = NULL
    )$pvalue,
    0
  )
})

test_that("TestCovariance_base dimensions", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1),
              nrow = 1,
              ncol = 21)
  Xi <- 2
  set.seed(31415)
  expect_error(
    test_correlation(
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
    test_correlation(
      X = X,
      nv = NULL,
      C = C[1, 1:20, drop = FALSE],
      Xi = Xi,
      method = "BT",
      repetitions = 1000
    )
  )
})

test_that("test_correlation wrong method", {
  C <- matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
              nrow = 1,
              ncol = 15)
  Xi <- 2
  set.seed(31415)
  expect_error(
    test_correlation(
      X = X,
      nv = NULL,
      C = C,
      Xi = Xi,
      method = "abc",
      repetitions = 1000
    )
  )
})
