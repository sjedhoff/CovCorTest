% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_functions.R
\name{CovTest}
\alias{CovTest}
\title{CovTest Object}
\usage{
CovTest()
}
\value{
An object of class \code{CovTest}.
}
\description{
This help page describes the structure of the \code{CovTest} class,
which is used to represent the results of a covariance and correlation test.
}
\details{
A \code{CovTest} object is a list with the following components:
\describe{
\item{\code{method}}{Character. Either 'Covariance' or 'Correlation'.}
\item{\code{pvalue}}{Numeric. The p-value of the test.}
\item{\code{Teststatistic}}{Numeric. The test statistic.}
\item{\code{CovarianceMatrix}}{Matrix. The covariance estimator for the
teststatistic.}
\item{\code{C}}{Numeric. A constant or vector of constants used in the
test.}
\item{\code{Xi}}{Numeric. A parameter related to the test.}
\item{\code{resampling_method}}{Character. The resampling method used in
the test.}
\item{\code{repetitions}}{Integer. The number of repetitions used in
resampling.}
\item{\code{hypothesis}}{Character. The hypothesis being tested.}
\item{\code{nv}}{Numeric. The sample size or the number of variables.}
}
}
\keyword{internal}
