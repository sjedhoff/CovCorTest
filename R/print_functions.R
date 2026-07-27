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
#'    \item{\code{C}}{Numeric matrix. The hypothesis matrix used for the
#'   computation of the test statistic. If \code{AM = 1}, this may be the
#'   alternative hypothesis matrix used internally.}
#'   \item{\code{Xi}}{Numeric vector. The hypothesis vector used together with
#'   \code{C}. If \code{AM = 1}, this may be the transformed hypothesis vector
#'   used internally.}
#'   \item{\code{AM}}{Integer. Indicates whether the alternative hypothesis
#'   matrix was used.}
#'   \item{\code{C_original}}{Numeric matrix. The original hypothesis matrix
#'   before the optional transformation by \code{AM}, if applicable.}
#'   \item{\code{Xi_original}}{Numeric vector. The original hypothesis vector
#'   before the optional transformation by \code{AM}, if applicable.}
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
                 "AM" = integer(1),
                 "C_original" = matrix(),
                 "Xi_original" = numeric(1),
                 "resampling_method" = character(1),
                 "repetitions" = integer(1),
                 "hypothesis" = character(1),
                 "nv" = numeric(1)),
            class = "CovTest")
}



#' Combined covariance and correlation test result
#'
#' Constructs an object of class \code{CombTest}. Objects of this class
#' represent results from the combined test for equality of covariance and
#' correlation matrices.
#'
#' A \code{CombTest} object contains the following components:
#' \describe{
#'   \item{\code{pvalue_Variances}}{
#'     Numeric. The p-value for equality of the variances.
#'   }
#'   \item{\code{pvalue_Correlations}}{
#'     Numeric. The p-value for equality of the correlations.
#'   }
#'   \item{\code{pvalue_Total}}{
#'     Numeric. The p-value for the global combined hypothesis.
#'   }
#'   \item{\code{Teststatistic}}{
#'     Numeric. The value of the combined test statistic.
#'   }
#'   \item{\code{repetitions}}{
#'     Integer. The number of resampling repetitions.
#'   }
#'   \item{\code{nv}}{
#'     Integer vector. The sample sizes of the two groups.
#'   }
#' }
#'
#' @param pvalue_Variances Numeric. The p-value for equality of the variances.
#' @param pvalue_Correlations Numeric. The p-value for equality of the
#'   correlations.
#' @param pvalue_Total Numeric. The p-value for the global combined hypothesis.
#' @param Teststatistic Numeric. The value of the combined test statistic.
#' @param repetitions Integer. The number of resampling repetitions.
#' @param nv Integer vector containing the sample sizes of the two groups.
#'
#' @return An object of class \code{CombTest}.
#'
#' @examples
#' result <- CombTest(
#'   pvalue_Variances = 0.12,
#'   pvalue_Correlations = 0.08,
#'   pvalue_Total = 0.08,
#'   Teststatistic = 2.4,
#'   repetitions = 1000L,
#'   nv = c(20L, 25L)
#' )
#'
#' @export
CombTest <- function(
    pvalue_Variances = NA_real_,
    pvalue_Correlations = NA_real_,
    pvalue_Total = NA_real_,
    Teststatistic = NA_real_,
    repetitions = NA_integer_,
    nv = integer()
) {
  structure(
    list(
      pvalue_Variances = pvalue_Variances,
      pvalue_Correlations = pvalue_Correlations,
      pvalue_Total = pvalue_Total,
      Teststatistic = Teststatistic,
      repetitions = repetitions,
      nv = nv
    ),
    class = "CombTest"
  )
}

#' Format a resampling-based p-value
#'
#' Formats a p-value obtained from a finite number of resampling repetitions.
#' If none of the resampled statistics is at least as extreme as the observed
#' statistic, the estimated p-value is zero and is displayed using the
#' resolution bound \code{1 / repetitions}.
#'
#' @param pvalue A single numeric p-value between zero and one.
#' @param repetitions A single positive integer specifying the number of
#'   resampling repetitions.
#' @param digits A single positive integer specifying the number of significant
#'   digits used for a nonzero p-value.
#'
#' @return A character string containing the formatted p-value.
#'
#' @keywords internal
format_resampling_pvalue <- function(pvalue, repetitions, digits = 4L) {
  if (
    !is.numeric(pvalue) ||
    length(pvalue) != 1L ||
    (!is.na(pvalue) && (
      !is.finite(pvalue) ||
      pvalue < 0 ||
      pvalue > 1
    ))
  ) {
    stop(
      "pvalue must be NA or a single finite number between 0 and 1."
    )
  }

  if (
    !is.numeric(digits) ||
    length(digits) != 1L ||
    is.na(digits) ||
    !is.finite(digits) ||
    digits <= 0 ||
    digits != floor(digits)
  ) {
    stop("digits must be a single positive integer.")
  }

  digits <- as.integer(digits)

  # Permit empty result objects such as CombTest().
  if (is.na(pvalue)) {
    return("p = NA")
  }

  # The number of repetitions is only needed for a zero p-value.
  if (pvalue == 0) {
    repetitions <- CheckRepetitions(
      repetitions,
      minimum_recommended = 1L
    )

    bound <- format(
      1 / repetitions,
      scientific = TRUE,
      digits = digits,
      trim = TRUE
    )

    return(paste0("p < ", bound))
  }

  paste0(
    "p = ",
    formatC(
      pvalue,
      format = "fg",
      digits = digits,
      flag = "#"
    )
  )
}


#' Print function for CovTest object
#'
#' @param x an \code{\link{CovTest}}  object
#' @param ... additional parameters
#'
#' @return no return, just print
#' @export
print.CovTest <- function(x, ...){
  method_print <- ifelse(x$resampling_method == "MC", "Monte-Carlo-technique",
                         ifelse(x$resampling_method == "BT", "Bootstrap",
                                "Taylor-based Monte-Carlo-approach"))


  pval <- format_resampling_pvalue(
    pvalue = x$pvalue,
    repetitions = x$repetitions
  )

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
#' @export
print.CombTest <- function(x, ...){

  pvalCov <- format_resampling_pvalue(
    pvalue = x$pvalue_Variances,
    repetitions = x$repetitions
  )

  pvalCorr <- format_resampling_pvalue(
    pvalue = x$pvalue_Correlations,
    repetitions = x$repetitions
  )

  pvalTotal <- format_resampling_pvalue(
    pvalue = x$pvalue_Total,
    repetitions = x$repetitions
  )



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

