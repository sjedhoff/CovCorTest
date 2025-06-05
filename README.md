
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CovCorTest

## Statistical Tests for Covariance and Correlation Matrices and their Structures

<!-- badges: start -->
[![R-CMD-check](https://github.com/sjedhoff/CovCorTest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjedhoff/CovCorTest/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/sjedhoff/CovCorTest/graph/badge.svg)](https://app.codecov.io/gh/sjedhoff/CovCorTest)
<!-- badges: end -->


A compilation of tests for hypotheses regarding covariance and
correlation matrices for one or more groups. The hypothesis can be
specified through a corresponding hypothesis matrix and a vector or by
choosing one of the basic hypotheses, while for the structure test, only
the latter works. Thereby Monte-Carlo and Bootstrap-techniques are used,
and the respective method must be chosen, and the functions provide
p-values and mostly also estimators of calculated covariance matrices of
test statistics.

## Installation

The official version of CovCorTest can be installed using the R Console:

``` r
install.packages("CovCorTest")
```

You can install the development version of CovCorTest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sjedhoff/CovCorTest")
```

## Structure of the package

The package is structures in tests regarding the covariance matrix and
the correlation matrix and their structures. A combined test for both is implemented as well.
For each of the matrices, covariance and correlation, two test functions are defined. 
The best approach is to start with the simple functions:
`test_covariance` and `test_correlation` respectively allow to test for a selection
of different predefined hypotheses for the corresponding matrices. These function
take the dataset, the group sizes (when testing for multiple groups) and
the hypothesis, which should be tested. Since the hypothesis can be
chosen using a character string like “equal”, no further knowledge about
the matrices used to test the hypotheses is needed. 
For more advanced users, alternatively to using the `hypothesis` argument,
a specific hypothesis matrix `C` and a corresponding vector `Xi` can be passed
along to the function. This can be used to test all forms of hypotheses, but in-depth
knowledge is necessary.

The structures of the covariance and correlation matrices can be tested
using `test_covariance_structure`and `test_correlation_structure`
respectively. Instead of a hypothesis, a structure can be selected using
a string, which will then be tested.

The `combined_test` functions delivers a possibility to test for equality of the covariance
and correlation matrix of two groups.

## Example

We are using the `EEGwide` dataset from the `MANOVA.RM` package as an example.
For this, we are just focusing on two groups and the numerical variables.
``` r
library(CovCorTest)
data("EEGwide", package = "MANOVA.RM")
vars <- colnames(EEGwide)[1:6]

data <- list(t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "AD", vars]),
               t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "MCI", vars]))
```

For the two groups, we can check for equality of the covariance matrices
```r
test_covariance(X = data, nv = c(12,27), hypothesis = "equal")

```
The nv argument is for passing along group sizes. We can also leave it empty
and a warning message shows.


We could also test, if the two groups are equal-correlated
```r
test_correlation(X = data, hypothesis = "equal-correlated")
```

With the combined test, we can test for the covariance and the correlation matrices
```r
test_combined(X = data, nv = c(12, 27))
```

The test for the structure of the covariance and correlation matrices are just
for one matrix, i.e. just one group. Different structures can be tested:
```r
test_covariance_structure(X = data[[1]], structure = "diag")
```

## Literature

- Sattler, P. & Dobler, D. (2025). Testing for patterns and structures
  in covariance and correlation matrices. <em>arXiv preprint</em>
  <https://arxiv.org/abs/2310.11799>
- Sattler, P., Bathke, A.C. & Pauly, M. (2022). Testing hypotheses about covariance 
  matrices in general MANOVA designs. <em>Journal of
  Statistical Planning and Inference</em> 219, 134-146
  <https://doi.org/10.1007/s11749-023-00906-6>
- Sattler, P. & Pauly, M. (2024). Testing hypotheses about correlation
  matrices in general MANOVA designs. <em>TEST</em> 33, 496–516
  <https://doi.org/10.1007/s11749-023-00906-6>
