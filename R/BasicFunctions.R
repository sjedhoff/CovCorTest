

#' Anova-Type-Statistic with weighted sum
#'
#' @description
#' Calculation of a Anova-Type-Statistic using
#'
#' @param A a matrix
#' @param repetitions a scalar, number of runs
#' @return a vector of the length of repetitions
#'
#' @export
#' @keywords internal
ATSwS <- function(A, repetitions){
  Chi <- matrix(stats::rchisq(dim(A)[1] * repetitions, df = 1), ncol = repetitions)
  return(colSums(crossprod(eigen(A, only.values = 1)$value, Chi)) / sum(diag(A)))
}


#' Anova-Type-statistic
#'
#' @param N number of observations
#' @param vVarData a matrix of vectorized covariance/correlation data
#' @param C hypothesis matrix for calculating the ATS
#' @param HatCov covariance matrix
#' @param Xi a vector defining together with C the investigated hypothesis
#'
#' @return a vector
#'
#' @export
#' @keywords internal
ATS <- function(N, vVarData, C, HatCov, Xi = 0){
  CDiff <-  C %*% vVarData - Xi
  statisticATS <-  N * crossprod(CDiff) / (sum(diag(QF(C, HatCov))))
  return(as.numeric(statisticATS))
}

#' Function to generate bootstrap observations
#'
#' @param WSigma weight matrix
#' @param nv number of observations in the groups
#'
#' @return a matrix
#'
#' @export
#' @keywords internal
generateData <- function(WSigma, nv){
  data <- WSigma %*% matrix(stats::rnorm(dim(WSigma)[1] * nv), ncol = nv)
  return(data)
}
