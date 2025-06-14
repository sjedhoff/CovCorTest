#' @title Root transformation of the vectorized correlation matrix
#'
#' @description A function calculating the roots of a vectorized covariance
#' matrix. The roots increasing, so square root for the first secondary
#' diagonals, third root for the second secondary diagonal and so on. For roots
#' with even order the absolute value of the argument is used, since the
#' arguments can be negative.
#'
#' @param v vectorized correlation matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' correlation matrix
#' @param d dimension of the correlation matrix
#' @return a transformed vector
#'
#' @keywords internal
#'
subdiagonal_mean_ratio_cor <- function(v, a, d){
  ratio <- rep(0, d - 2)
  ae <-  c(a, a[d] + 1)
  for(l in 3:(d)){
    ratio[l - 2] <-  mean(v[(ae[l]-d):(ae[l + 1] - 1-d)]) /
      mean(v[(ae[l - 1]-d):(ae[l] - 1-d)])
  }
  return(c(v,ratio))
}







#' Function for the Taylor-based Monte-Carlo-approximation for one group
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for one group. After receiving some
#' auxiliary matrices and data, the Monte-Carlo observations are generated and
#' different parts of the final sum are defined. Based on this a number of the
#' Taylor-based ATS are calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the
#' approximation
#' @param C the used hypothesis matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param CorData the calculated correlation matrix
#' @param MvrH an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M4 a auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param n1 the total sample size, a scalar
#' @param Atilde an auxiliary matrix for the transformation from row-wise
#' vectorization
#' @param Jacobi the Jacobian matrix of the transformation function applied
#' for the diagonal vectorized correlation  to diagonalwise vectorization. used
#' for the transformation function 'subdiagonal_mean_ratio_cor', else NULL
#' @return a matrix containing the values of the Taylor ATS for a number of
#' repetitions
#'
#' @keywords internal
#'
Tayapp1G <- function(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M4, L,
                     P, Q, n1, Atilde = NULL, Jacobi = NULL){
  vechCorData <- matrixcalc::vech(CorData)
  DvechCorDataM <- as.vector(vechCorData) * M4
  XPB <- generateData(MSrootStUpsi, repetitions)

  PUX <- P %*% XPB
  MUX0 <- M4 %*% Q %*% XPB

  HX <- Qvech(PUX, repetitions)

  Part1 <- MvrH %*% XPB
  Part2 <- 1 / (4) * as.vector(vechCorData) * HX
  Part3 <- 1 / (2) * XPB * MUX0
  Part4 <- 3 / (8) * DvechCorDataM %*% HX

  # without the transformation
  if(is.null(Atilde)){
    CXTay <- C %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1)) / sqrt(n1)
    Result <- n1 * apply(CXTay, 2, crossprod) / Trace
  }
  # with the transformation
  else{
    XTaydv <- Atilde %*% (Part1 + L%*%(Part2 - Part3 + Part4)/sqrt(n1))/sqrt(n1)
    CXTaydv <- C %*% Jacobi %*% XTaydv
    Result <- n1 * apply(CXTaydv, 2, crossprod) / Trace
  }
  return(Result)
}

#' Function for the Taylor-based Monte-Carlo-approximation for multiple groups
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for multiple groups. After receiving
#' some auxiliary matrices and data, the Monte-Carlo observations are generated
#' and different parts of the final sum are defined. Based on this a number of
#' the Taylor-based ATS are calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the
#' approximation
#' @param C the used hypothesis matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param CorData the calculated correlation matrix
#' @param MvrH an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param MvrH an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M4 an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param nv vector of sample sizes
#' @return a matrix containing the values of the Taylor ATS for a number
#' of repetitions
#'
#' @keywords internal
#'
#'
TayappMG <- function(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M4, L,
                     P, Q, nv){
  vechCorData <- lapply(CorData, matrixcalc::vech)
  DvechCorDataM <- lapply(lapply(vechCorData, as.vector), "*", M4)

  XPB <- mapply(generateData, MSrootStUpsi, repetitions, SIMPLIFY = FALSE)

  PUX <- lapply(XPB, function(X) P %*% X )
  MUX0 <- lapply(XPB, function(X) (M4 %*% Q) %*% X )
  HX <- lapply(PUX, Qvech, repetitions)

  Part1 <- mapply(function(A, B) A %*% B, MvrH, XPB, SIMPLIFY = FALSE)
  Part2 <- mapply("*", lapply(vechCorData, as.vector), HX, SIMPLIFY = FALSE)
  Part3 <- mapply("*", XPB, MUX0, SIMPLIFY = FALSE)
  Part4 <- mapply("%*%", DvechCorDataM, HX, SIMPLIFY = FALSE)
  XTay <- mapply(function(Part1, L, Part2, Part3, Part4, nv){
    return((Part1 + L %*% (1 / 4 * Part2 - 1 / 2 * Part3 + 3 / 8 * Part4) /
              sqrt(nv)) / sqrt(nv))
  }, Part1, list(L), Part2, Part3, Part4, nv, SIMPLIFY = FALSE)

  XTayMatrix <- do.call(rbind, XTay)

  CXTay <- C %*% XTayMatrix
  Result <- sum(nv) *apply(CXTay, 2, crossprod) / Trace
  return(Result)
}


