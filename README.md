
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CovCorTest

## Statistical Tests for Covariance and Correlation Matrices and their Structures

<!-- badges: start -->

[![R-CMD-check](https://github.com/sjedhoff/CovCorTest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjedhoff/CovCorTest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/sjedhoff/CovCorTest/graph/badge.svg)](https://app.codecov.io/gh/sjedhoff/CovCorTest)
[![CRAN
status](https://www.r-pkg.org/badges/version/CovCorTest)](https://CRAN.R-project.org/package=CovCorTest)
<!-- badges: end -->

A compilation of tests for hypotheses regarding covariance and
correlation matrices for one or more groups. The hypothesis can be
specified by choosing one of the predefined hypotheses or, for more
advanced applications, through a corresponding hypothesis matrix and a
vector of null values.

The package implements Monte Carlo, bootstrap, and Taylor approximation
procedures, depending on the selected test. The functions provide
p-values, test statistics, and additional information on the tested
hypotheses.

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

## Citation

If you use `CovCorTest` in a scientific publication, please cite:

Sattler, P. and Jedhoff, S. (2025). Testing Hypotheses regarding
Covariance and Correlation Matrices with the R Package CovCorTest.
*arXiv:2507.03406*. <https://doi.org/10.48550/arXiv.2507.03406>

The corresponding citation and BibTeX entry can also be obtained
directly in R:

``` r
citation("CovCorTest")
toBibtex(citation("CovCorTest"))
```

Depending on the methods used in the analysis, please additionally cite
the corresponding methodological paper listed in the [Methodological
references](#methodological-references) section below.

## Structure of the package

The package provides tests for covariance matrices, correlation
matrices, and their structures. A combined test for equality of both
covariance and correlation matrices of two groups is implemented as
well.

For covariance and correlation matrices, the main functions are
`test_covariance` and `test_correlation`. These functions allow users to
test a selection of predefined hypotheses by specifying the data, the
group sizes when required, and the hypothesis of interest. Since
predefined hypotheses can be selected using character strings such as
`"equal"` or `"equal-correlated"`, no direct specification of hypothesis
matrices is needed for standard applications.

For more advanced users, a custom hypothesis matrix `C` and a
corresponding vector `Xi` can be supplied instead of using a predefined
hypothesis. This allows for more general hypotheses, but requires
knowledge of the underlying matrix representation.

By default, the package uses alternative hypothesis matrices via the
argument `AM = 1`. These alternative matrices reduce the number of rows
of the hypothesis matrix while leaving the resulting test statistic
unchanged. The original hypothesis matrix and vector are also stored in
the output as `C_original` and `Xi_original`.

The structures of covariance and correlation matrices can be tested
using `test_covariance_structure` and `test_correlation_structure`,
respectively. Instead of a hypothesis, a structure is selected using a
character string.

For covariance matrices, implemented structures include autoregressive,
first-order autoregressive, diagonal, sphericity, compound symmetry,
Toeplitz, constant off-diagonal, standard Toeplitz, banded, and banded
Toeplitz structures.

For correlation matrices, implemented structures include heterogeneous
autoregressive, heterogeneous compound symmetry, heterogeneous Toeplitz,
diagonal, banded, and banded Toeplitz structures.

For banded and banded Toeplitz structures, the additional argument
`bandwidth` specifies the number of off-diagonals that are allowed to be
non-zero.

The function `test_combined` provides a combined test for equality of
covariance and correlation matrices of two groups.

## Example

We are using the `EEGwide` dataset from the `MANOVA.RM` package as an
example. For this, we focus on two groups and the first six numerical
variables.

``` r
library(CovCorTest)
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]

data <- list(
  EEGwide[EEGwide$sex == "M" &
            EEGwide$diagnosis == "AD", vars],
  EEGwide[EEGwide$sex == "M" &
            EEGwide$diagnosis == "MCI", vars]
)
```

The observations are stored in the rows and the variables in the
columns. For the two groups, we can test for equality of the covariance
matrices.

``` r
test_covariance(X = data, hypothesis = "equal")
```

Since `data` is supplied as a list, the group sizes are obtained from
the individual list entries. Alternatively, the group sizes can be
supplied explicitly using the `nv` argument.

``` r
test_covariance(X = data, nv = c(12, 27), hypothesis = "equal")
```

We can also test whether the two groups have equal correlation matrices.

``` r
test_correlation(X = data, hypothesis = "equal-correlated")
```

With the combined test, we can test for equality of both covariance and
correlation matrices.

``` r
test_combined(X = data, nv = c(12, 27))
```

The structure tests are one-sample tests. For example, we can test
whether the covariance matrix of the first group is diagonal.

``` r
test_covariance_structure(X = data[[1]], structure = "diag")
```

Banded structures require the additional `bandwidth` argument.

``` r
test_covariance_structure(
  X = data[[1]],
  structure = "banded",
  bandwidth = 2
)
```

Analogously, we can test structures of correlation matrices.

``` r
test_correlation_structure(
  X = data[[1]],
  structure = "banded-toeplitz",
  bandwidth = 2
)
```

## Methodological references

The following articles provide the statistical methodology implemented
in `CovCorTest`:

- **Covariance matrix hypotheses (`test_covariance`):**  
  Sattler, P., Bathke, A. C. & Pauly, M. (2022). Testing hypotheses
  about covariance matrices in general MANOVA designs. *Journal of
  Statistical Planning and Inference* 219, 134–146.
  <https://doi.org/10.1016/j.jspi.2021.12.001>

- **Correlation matrix hypotheses and the combined test
  (`test_correlation`, `test_combined`):**  
  Sattler, P. & Pauly, M. (2024). Testing hypotheses about correlation
  matrices in general MANOVA designs. *TEST* 33, 496–516.
  <https://doi.org/10.1007/s11749-023-00906-6>

- **Covariance and correlation structure tests
  (`test_covariance_structure`, `test_correlation_structure`):**  
  Sattler, P. & Dobler, D. (2026). Testing for patterns and structures
  in covariance and correlation matrices. *Journal of Multivariate
  Analysis* 211, 105517. <https://doi.org/10.1016/j.jmva.2025.105517>
