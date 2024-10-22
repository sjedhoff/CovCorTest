
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CovCorTest

## Statistical Tests for Covariance and Correlation Matrices and their Structures

<!-- badges: start -->
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

The offical version of CovCorTest can be installed using the R Console:

``` r
install.packages("CovCorTest")
#> Installiere Paket nach 'C:/Users/sjedh/AppData/Local/Temp/RtmpCYyLbV/temp_libpath282074287f3f'
#> (da 'lib' nicht spezifiziert)
#> Warning: package 'CovCorTest' is not available for this version of R
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

You can install the development version of CovCorTest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sjedhoff/CovCorTest")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo sjedhoff/CovCorTest@HEAD
#> sandwich   (3.1-0  -> 3.1-1 ) [CRAN]
#> mvtnorm    (1.2-4  -> 1.3-1 ) [CRAN]
#> Rcpp       (1.0.10 -> 1.0.13) [CRAN]
#> rbibutils  (2.2.16 -> 2.3   ) [CRAN]
#> data.table (1.14.8 -> 1.16.2) [CRAN]
#> multcomp   (1.4-25 -> 1.4-26) [CRAN]
#> plotrix    (3.8-2  -> 3.8-4 ) [CRAN]
#> plyr       (1.8.8  -> 1.8.9 ) [CRAN]
#> Installing 8 packages: sandwich, mvtnorm, Rcpp, rbibutils, data.table, multcomp, plotrix, plyr
#> Installiere Pakete nach 'C:/Users/sjedh/AppData/Local/Temp/RtmpCYyLbV/temp_libpath282074287f3f'
#> (da 'lib' nicht spezifiziert)
#> Paket 'sandwich' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'mvtnorm' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'Rcpp' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'rbibutils' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'data.table' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'multcomp' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'plotrix' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'plyr' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> 
#> Die heruntergeladenen Binärpakete sind in 
#>  C:\Users\sjedh\AppData\Local\Temp\RtmpiYrDPt\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\sjedh\AppData\Local\Temp\RtmpiYrDPt\remotes2eb047671f86\sjedhoff-CovCorTest-9f3daa6235923210e940154cb30d510cfc20b70a/DESCRIPTION' ...     checking for file 'C:\Users\sjedh\AppData\Local\Temp\RtmpiYrDPt\remotes2eb047671f86\sjedhoff-CovCorTest-9f3daa6235923210e940154cb30d510cfc20b70a/DESCRIPTION' ...   ✔  checking for file 'C:\Users\sjedh\AppData\Local\Temp\RtmpiYrDPt\remotes2eb047671f86\sjedhoff-CovCorTest-9f3daa6235923210e940154cb30d510cfc20b70a/DESCRIPTION' (415ms)
#>       ─  preparing 'CovCorTest': (728ms)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  installing the package to process help pages
#>      Lade nötigen Namensraum: CovCorTest
#>       ─  saving partial Rd database (1.5s)
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>       ─  building 'CovCorTest_0.0.1.tar.gz'
#>      
#> 
#> Installiere Paket nach 'C:/Users/sjedh/AppData/Local/Temp/RtmpCYyLbV/temp_libpath282074287f3f'
#> (da 'lib' nicht spezifiziert)
```

## Structure of the package

The package is structures in tests regarding the covariance matrix and
the correlation matrix and their structures. For this, each of the
matrices, covariance and correlation, has eac three different test
functions. The best approach is to start with the simple functions:
`TestCovariance_simple` and `TestCorrelation_simple`. These function
take the dataset, the group sizes (when testing for multiple groups) and
the hypothesis, which should be tested. Since the hypothesis can be
chosen using a character string like “equal”, no further knowledge about
the matrices used to test the hypotheses is needed. The advanced methods
are `TestCovariance_base` and `TestCorrelation_base`, where the
hypothesis is selected using a hypothesis matrix and a corresponding
vector. This can be used to test all forms of hypotheses, but in-depth
knowledge is necessary.

The structures of the covariance and correlation matrices can be tested
using `TestCovariance_structure`and `TestCorrelation_structure`
respectively. Instead of a hypothesis, a structure can be selected using
a string, which will then be tested.

## Example

``` r
library(CovCorTest)
## basic example code
```

## Literature
