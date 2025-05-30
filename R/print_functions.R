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
#'   \item{\code{CovarianceMatrix}}{Matrix. The covariance estimator for the
#'   teststatistic.}
#'   \item{\code{C}}{Numeric. A constant or vector of constants used in the
#'   test.}
#'   \item{\code{Xi}}{Numeric. A parameter related to the test.}
#'   \item{\code{resampling_method}}{Character. The resampling method used in
#'   the test.}
#'   \item{\code{repetitions}}{Integer. The number of repetitions used in
#'   resampling.}
#'   \item{\code{hypothesis}}{Character. The hypothesis being tested.}
#'   \item{\code{nv}}{Numeric. The sample size or the number of variables.}
#' }
#'
#' @return An object of class \code{CovTest}.
#' @export
#' @keywords internal
CovTest <- function() {
  structure(list("method" = character(1),
                 "pvalue" = 0.1,
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



#' CombTest Object
#'
#' This help page describes the structure of the \code{CombTest} class,
#' which is used to represent the results the combined covariance and
#' correlation test.
#'
#' A \code{CombTest} object is a list with the following components:
#' \describe{
#'   \item{\code{method}}{Character. Either 'Covariance' or 'Correlation'.}
#'   \item{\code{pvalue-Variances}}{Numeric. The p-value of the test regarding
#'   the covariances.}
#'   \item{\code{pvalue-Correlations}}{Numeric. The p-value of the test
#'   regarding the correlations.}
#'   \item{\code{pvalue-Total}}{Numeric. The p-value of the whole test of the
#'   global hypothesis.}
#'   \item{\code{Teststatistic}}{Numeric. The test statistic.}
#'   \item{\code{resampling_method}}{Character. The resampling method used in
#'   the test.}
#'   \item{\code{repetitions}}{Integer. The number of repetitions used in
#'   resampling.}
#'   \item{\code{nv}}{Numeric. The sample size or the number of variables.}
#' }
#'
#' @return An object of class \code{CombTest}.
#' @export
#' @keywords internal
CombTest <- function() {
  structure(list("pvalue_Variances" = 0.1,
                 "pvalue_Correlations" = 0.1,
                 "pvalue_Total" = 0.1,
                 "Teststatistic" = numeric(1),
                 #"CovarianceMatrix" = matrix(),
                 "repetitions" = integer(1),
                 "nv" = numeric(1)),
            class = "CombTest")
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
                         ifelse(x$resampling_method == "BT", "Bootstrap",
                                "Taylor-based Monte-Carlo-approach"))


  pval <- ifelse(x$pvalue < 10^(-4), paste0("p < 1e-", log10(x$repetitions)),
                 paste0("p = ", round(x$pvalue, digits = 4)))

  group_text <- ifelse(length(x$nv) == 1, "one group", paste0("",length(x$nv),
                                                              " groups"))

  cat("\n
       \t ",x$method," Test \n \t    ",group_text,"\n\n Hypothesis: \t\t",
      x$hypothesis,
      "\n Teststatistic value: \t",
      round(x$Teststatistic, digits = 4),
      "\n p-value: \t \t",
      pval,
      "\n \n p-value computed using ", method_print, " with B=", x$repetitions,
      " repetitions \n",
      sep = "")
}

#' Print function for CombTest object
#'
#' @param x an \code{\link{CombTest}}  object
#' @param ... additional parameters
#'
#' @return no return, just print
#' @exportS3Method
print.CombTest <- function(x, ...){

  pvalCov <- ifelse(x$pvalue_Variances < 10^(-4), paste0("p < 1e-",
                log10(x$repetitions)),
                paste0("p = ", round(x$pvalue_Variances, digits = 4)))
  pvalCorr <- ifelse(x$pvalue_Correlations < 10^(-4), paste0("p < 1e-",
                log10(x$repetitions)),
                paste0("p = ", round(x$pvalue_Correlations, digits = 4)))
  pvalTotal <- ifelse(x$pvalue_Total < 10^(-4), paste0("p < 1e-",
                log10(x$repetitions)),
                paste0("p = ", round(x$pvalue_Total, digits = 4)))



  cat("\n
       \t Combined Test \n",
      "\n p-value-Variances: \t \t",
      pvalCov,
      "\n p-value-Correlations: \t \t",
      pvalCorr,
      "\n p-value-Total: \t \t",
      pvalTotal,
      "\n \n p-values computed using Taylor-based Monte-Carlo-approach with B=",
      x$repetitions, " repetitions \n",
      sep = "")
}

