#' @title  Transformation of the vectorized covariance matrix by
#' quotients of means
#'
#' @description A function which calculates the mean of the secondary diagonals
#' and divide them through the next one. Since the elements can be negative, for
#' the denominator absolute values are used.
#' @param v vectorized covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#'
#' @keywords internal
#'
subdiagonal_mean_ratio_fct <- function(v, a, d){
  ratio <- rep(0, d - 1)
  ae <-  c(a, a[d] + 1)
  for(l in 2:(d)){
    ratio[l - 1] <-  mean(v[ae[l]:(ae[l + 1] - 1)]) /
      mean(v[ae[l - 1]:(ae[l] - 1)])
  }
  return(c(v, ratio))
}


#' @title Jacobian matrix for transformation functions
#'
#' @description A function which calculates the Jacobian matrix for a given
#' transformation function \code{\link{subdiagonal_mean_ratio_cor}} or
#' \code{\link{subdiagonal_mean_ratio_fct}}
#' @param X vectorized covariance matrix for which the Jacobian matrix is
#' applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorized matrix
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or
#'  \code{\link{subdiagonal_mean_ratio_cor}}
#' @return the Jacobian matrix applied for the given vector
#'
#' @keywords internal
#'
Jacobian <- function(X, a, d, p, fun){
    if(fun == "subdiagonal_mean_ratio_fct"){
      J <-  matrix(0, d - 1, p)
      for(l in 1:(d - 1)){
        S1 <- sum((X[a[l]:(a[l] + d - l)]))
        S2 <- sum(X[a[l + 1]:(a[l + 1] + d - l - 1)])

        J[l, a[l + 1] + 0:(d - l - 1)] <- (d - l + 1) / (d - l) / S1
        J[l, a[l] + 0:(d - l)] <- (d - l + 1) / (d - l)  * (-S2 / (S1^2))
      }
      return(rbind(diag(1, p, p), J))
    }
    else{
      if(fun == "subdiagonal_mean_ratio_cor"){
        J <-  matrix(0, d - 2, p-d)
        for (l in 2:(d - 1)){
          S1 <- sum((X[(a[l]-d):(a[l] - l)]))
          S2 <- sum(X[(a[l + 1]-d):(a[l + 1]- l - 1)])

          J[l-1, (a[l + 1] + 0:(d - l - 1))-d] <- (d - l + 1) / (d - l) / S1
          J[l-1, (a[l] + 0:(d - l))-d] <- (d - l + 1) / (d - l)  *(-S2 / (S1^2))
        }
        return(rbind(diag(1, p-d, p-d), J))

      }

      else{
        stop("fun must be 'subdiagonal_mean_ratio_fct',
             'subdiagonal_mean_ratio_cor'")
      }
    }

}


#' @title ATS for transformed vectors
#'
#' @description A function which calculates the Anova-type-statistic based on
#' a transformation function
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param v vectorized empirical covariance matrix of the original data
#' @param a vector containing the indices which belongs to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorized matrix
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or
#'  \code{\link{subdiagonal_mean_ratio_cor}}
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#'
ATS_fun <- function(N, X, C, v, a, d, p, fun){
  Xmean <-  rowMeans(X)
  CDiff <- C %*% (do.call(fun, list(Xmean, a, d)) - do.call(fun, list(v, a, d)))
  Jacobi <-  Jacobian(Xmean, a, d, p, fun)
  HatCov <-  stats::var(t(X))
  Trace <-  sum(diag(QF(C %*% Jacobi, HatCov)))
  return(c(N * crossprod(CDiff) / Trace))
}


#' @title Bootstrap using transformation for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix HatCov, which matrix root MSrootHatCov is given and
#' expectation vector vX. For the generated bootstrap sample the value of the
#' ATS based on a transformation is calculated
#' @param N.sim control variable for using sapply
#' @param n1 a scalar, declaring the sample size for the bootstrap sample
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorized matrix
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix root of the covariance matrix HatCov, to generate
#' the bootstrap sample
#' @param vX the expectation vector for the bootstrap sample
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or
#' \code{\link{subdiagonal_mean_ratio_cor}}
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#'
Bootstrap_trans <- function(N.sim, n1, a, d, p, C, MSrootHatCov, vX, fun){
  XPB <- generateData(MSrootHatCov, n1) + vX
  return(ATS_fun(n1, XPB, C, vX, a, d, p, fun))
}


#' @title Bootstrap for one and multiple groups
#'
#' @description This function generates normal distributed random vectors.
#' For one group, nv random vectors with covariance matrix HatCov are generated
#' and the corresponding value of the ATS is generated. For multiple groups the
#' corresponding sample sizes from nv are used.
#' The weighted sum of covariance matrices is calculated and used to calculate
#' the value of the ATS.
#' @param N.sim control variable for using sapply
#' @param nv scalar (one group) or vector (multiple groups) of sample sizes for
#' the bootstrap samples
#' @param C hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix (one group) or list of matrices (multiple groups)
#' of roots of the covariance matrices, to generate
#' the bootstrap sample
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#'
Bootstrap <- function(N.sim, nv, C, MSrootHatCov){
  # one group
  if(length(nv) == 1){
    XPB <- generateData(MSrootHatCov, nv)
    PBHatCov <- stats::var(t(XPB))
    return(ATS(nv, rowMeans(XPB), C, PBHatCov))
  }
  # multiple groups
  else{
    N <- sum(nv)
    kappainvv <- N / nv

    DataPB <- mapply(generateData, MSrootHatCov, nv, SIMPLIFY = FALSE)
    PBHatCov <- WDirect.sumL(lapply(DataPB,
                                    function(X) stats::var(t(X))), kappainvv)
    return(ATS(N, unlist(lapply(DataPB, rowMeans)), C, PBHatCov))
  }
}
