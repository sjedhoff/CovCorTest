#' CovTest Object
#'
#' This help page describes the structure of the \code{CovTest} class,
#' which is used to represent the results of a covariance and correlation test.
#'
#' A \code{CovTest} object is a list with the following components:
#' \describe{
#'   \item{\code{method}}{Character. Either 'Covariance' or 'Correlation'.}
#'   \item{\code{pvalue}}{Numeric. The p-value of the test.}
#'   \item{\code{Teststatistic}}{Numeric. The test statistic.}
#'   \item{\code{CovarianceMatrix}}{Matrix. The estimated covariance or correlation matrix.}
#'   \item{\code{C}}{Numeric. A constant or vector of constants used in the test.}
#'   \item{\code{Xi}}{Numeric. A parameter related to the test.}
#'   \item{\code{resampling_method}}{Character. The resampling method used in the test.}
#'   \item{\code{repetitions}}{Integer. The number of repetitions used in resampling.}
#'   \item{\code{hypothesis}}{Character. The hypothesis being tested.}
#'   \item{\code{nv}}{Numeric. The sample size or the number of variables.}
#' }
#'
#' @return An object of class \code{CovTest}.
#' @export
#' @keywords internal
CovTest <- function() {
  structure(list("method" = character(1),
                 "pvalue" = numeric(1),
                 "Teststatistic" = numeric(1),
                 "CovarianceMatrix" = matrix(),
                 "C" = numeric(1),
                 "Xi" = numeric(1),
                 "resampling_method" = character(1),
                 "repetitions" = integer(1),
                 "hypothesis" = character(1),
                 "nv" = numeric(1)),
            class = "CovTest")
}




#' Print function for CovTest object
#'
#' @param x an \code{\link{CovTest}}  object
#' @param ... additional parameters
#'
#' @return no return, just print
#' @exportS3Method
print.CovTest <- function(x, ...){
  method_print <- ifelse(x$resampling_method == "MC", "Monte-Carlo-technique",
                         ifelse(x$resampling_method == "BT", "Bootstrap", "Taylor-based Monte-Carlo-approach"))
  group_text <- ifelse(length(x$nv) == 1, "one group", paste0("",length(x$nv), "  groups"))

  cat("\n
       \t ",x$method," Test \n \t    ",group_text,"\n\n Hypothesis: \t\t",
  x$hypothesis,
  "\n Teststatistic value: \t",
  round(x$Teststatistic, digits = 4),
  "\n p-value: \t \t",
  round(x$pvalue, digits = 4),
  "\n \n p-value computed using ", method_print, " and n=", x$repetitions, "\n", sep = "")
}


