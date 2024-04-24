htilde <- function(x, a, d){
  xt <- x
  for(i in 3:d){
    if ((i %% 2 == 1)){
      xt[(0:(d - i)) + a[i - 1] - (i - 2)] <- abs(x[(0:(d - i)) + a[i - 1] - (i -2)]) ^ (1 / (i - 1))
    }
    if ((i %% 2 == 0)){
      xt[(0:(d - i)) + a[i - 1] - (i - 2)] <- (x[(0:(d - i)) + a[i - 1] - (i -2)] <= 0)
      * (-abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^ (1 / (i - 1)))
      + (x[(0:(d -i)) + a[i - 1] - (i - 2)] > 0) * (abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^(1 / (i - 1)))
    }
  }
  return(xt)
}

#' @title Jacobian matrix for the function htilde
#'
#' @description A function which calculates the Jacobian matrix for the root
#' transformation h
#' @param x vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return a transformed vector
#' @export
Jacobianhtilde <- function(x, a, d, p){
  E <- rep(1, p - d)
  for(i in 3:d){
    if((i %% 2 == 1)){
      E[(0:(d - i)) + a[i - 1] - (i - 2)] <-  x[(0:(d - i)) + a[i - 1] - (i - 2)] /
        ((i - 1) * abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^ (2 - 1 / (i - 1)))
    }
    if((i %% 2 == 0)){
      E[(0:(d - i)) + a[i - 1] - (i - 2)] <- 1 / ((i - 1) * abs(x[(0:(d - i)) +
                                                                   a[i - 1] - (i - 2)]) ^ (1 - (1 / (i - 1))))
    }
  }
  return(diag(E, p - d, p - d))
}


Bootstraphtilde <- function(N.sim, n1, a, d, p, C, Upsidv, vCorData){
  XPB <- gData(Upsidv, n1) + vCorData
  return(ATShtilde(n1, XPB, C, vCorData, a, d, p))
}

#' @title ATS for vectors transformed with the function htilde
#'
#' @description A function which calculates the Anova type statistic based on
#' the transformation htilde
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param r vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return a transformed vector
#' @export
ATShtilde <- function(N, X, C, r, a, d, p){
  Xmean <- rowMeans(X)
  CDiff <- C %*% (htilde(Xmean, a, d) - htilde(r, a, d))
  Jacobi <- Jacobianhtilde(Xmean, a, d, p)
  HatCov <- tvar(X)
  Trace <- sum(diag(QF(C %*% Jacobi, HatCov)))
  return(c(N * crossprod(CDiff) / Trace))
}


#' Function for the Taylor-based Monte-Carlo-approximation for one group
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for one group.  After receiving some
#' auxiliary matrices and data, the Monte-Carlo observations are generated and
#' different parts of the final sum are defined. Based on this a number of the
#' Taylor-based ATS are calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the approximation
#' @param C the used hypothesis matrix
#' @param p the dimension of the vectorized matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param Cordata the calculated correlation matrix
#' @param MvrH an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M a auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param N the total sample size, a scalar
#' @return a matrix containing the values of the Taylor ATS for a number of repetitions
Tayapp1G <- function(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M, L, P, Q, n1){
  vechCorData <- vech(CorData)
  DvechCorDataM <- as.vector(vechCorData) * M
  XPB <- gData(MSrootStUpsi, repetitions)

  PUX <- P %*% XPB
  MUX0 <- M %*% Q %*% XPB

  HX <- Qvech(PUX, repetitions)

  Part1 <- MvrH %*% XPB
  Part2 <- 1 / (4) * as.vector(vechCorData) * HX
  Part3 <- 1 / (2) * XPB * MUX0
  Part4 <- 3 / (8) * DvechCorDataM %*% HX
  CXTay <- C %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1)) / sqrt(n1)

  Result <- n1 * apply(CXTay, 2, crossprod) / Trace
  return(Result)
}

#' Function for the Taylor-based Monte-Carlo-approximation for multiple groups
#'
#' @description An auxiliary function to calculate the values for the Taylor-based
#' Monte-Carlo-approximation for one group.  After receiving some auxiliary matrices
#' and data, the Monte-Carlo observations are generated and different parts of the
#' final sum are defined. Based on this a number of the Taylor-based ATS are
#' calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the approximation
#' @param C the used hypothesis matrix
#' @param p the dimension of the vectorized matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param Cordata the calculated correlation matrix
#' @param MvrH an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param MvrH an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param nv vector of sample sizes
#' @return a matrix containing the values of the Taylor ATS for a number of repetitions
TayappMG <- function(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M, L, P, Q, nv){
  vechCorData <- lapply(CorData, vech)
  DvechCorDataM <- lapply(lapply(vechCorData, as.vector), "*", M)

  XPB <- mapply(gData, MSrootStUpsi, repetitions, SIMPLIFY = FALSE)

  PUX <- lapply(XPB, MprodBackward, P)
  MUX0 <- lapply(XPB, MprodBackward, (M %*% Q))
  HX <- lapply(PUX, Qvech, repetititons, SIMPLIFY = FALSE)

  Part1 <- mapply(MprodBackward, XPB, MvrH, SIMPLIFY = FALSE)
  Part2 <- mapply("*", lapply(vechCorData, as.vector), HX, SIMPLIFY = FALSE)
  Part3 <- mapply("*", XPB, MUX0, SIMPLIFY = FALSE)
  Part4 <- mapply("%*%", DvechCorDataM, HX, SIMPLIFY = FALSE)
  XTay <- mapply(function(Part1, L, Part2, Part3, Part4, nv) {
    return((Part1 + L %*% (1 / 4 * Part2 - 1 / 2 * Part3 + 3 / 8 * Part4) /
              sqrt(nv)) / sqrt(nv))
  },
  Part1, list(L), Part2, Part3, Part4, nv, SIMPLIFY = FALSE)
  CXTaymean <- C %*% unlist(lapply(XTay, rowMeans))
  Result <- sum(nv) * tcrossprod(CXTaymean) / Trace
  return(Result)
}






#' Function for the test regarding correlation matrices for one group.
#'
#' @description This function conduces the test for hypotheses regarding the
#' correlation matrix for one group. Depending on the chosen method a bootstrap
#' or Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param C hypothesis matrix for calculating the ATS
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or Monte-Carlo-technique("MC")
#' or a Taylor-based Monte-Carlos-approach("Tay") is used.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1.000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @export
TestCorrelation1G <- function(X, C, Xi, method, repetitions = 1000, seed = 0){
  if(seed != 0){
    (set.seed(seed))
  }

  n1 <- dim(X)[2]
  d <- dim(X)[1]
  p <- *(d+1)/2
  a <- cumsum(c(1,(d):2))
  H <- matrix(rep(a,d),d,d,byrow=TRUE)
  L <- diag(1,p,p)[-a,]
  M <- matrix(0,p,p)
  for(i in 1:p){
    M[i,c(vechp(H)[i],vechp(t(H))[i])] <- 1
  }
  M1 <- L%*%M
  VarData <- tvar(X)
  CorData <- cov2cor(VarData)
  vCorData <- vechp(CorData)
  Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod),nrow=p,ncol=n1)
  HatCov <- tvar(Xq)
  MvrH1 <- (L-1/2*vCorData*M1)
  MvrH2 <- sqrt(diag(as.vector(1/vtcrossprod(matrix(vech(VarData)[a]))),p,p))
  Upsi <- QF(MvrH1%*%MvrH2,HatCov)
  if(method=="MC"){ ResamplingResult <- ATSwS(QF(C,Upsi),repetitions)}
  if(method=="BT"){ ResamplingResult <- sapply(1:repetitions,Bootstrap1G,n1,C,MSroot(Upsi))}
  if(method=="Tay"){
    P <- diag(1,p,p)[a,]
    Q <- diag(as.vector(vech(diag(1,d,d))),p,p)
    StUpsi <- QF(MvrH2,HatCov)
    Trace <- sum(diag(QF(C,Upsi)))
    ResamplingResult <- Tayapp1G(repetitions,C,MSroot(StUpsi),CorData,MvrH1,Trace,M,L,P,Q,n1)
  }

  Teststatistic <- ATS(n1,vCorData,C,Upsi,Xi)
  pvalue <- mean(ResamplingResult<Teststatistic)
  if(seed!=0){
    (set.seed(NULL))
  }
  return(list("pvalue"=pvalue,"Teststatistic"=Teststatistic,"CovarianceMatrix"=Upsi))}


#' Simplified call of function for the test regarding correlation matrices for one group.
#'
#' @description This function is for more applied users so no hypothesis matrix
#' or corresponding vector is necessary. This is replaced by predefined
#' hypotheses, from which is chosen. From this C and Xi are built and the
#' function Testcorrelation1G is used.
#' @param X a matrix containing the observation vectors as columns
#' @param hypothesis a character to choose one of the predefined hypotheses which are
#' "equal-variances","uncorrelated","given-trace" and "given-matrix".
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#' is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1.000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @export
TestCorrelation1Gsimple <- function(X, hypothesis, method = "BT", repetitions = 1000, seed = 0){
  d <-  dim(X)[1]
  p <-  d * (d + 1) / 2
  if(d == 1){
    stop("Correlation is only defined for dimension higher than one")
  }
  if(d > 1){
    if((1 - (hypothesis %in% c("equal-correlated", "uncorrelated"))) == 1) {
      stop("no predefined hypothesis")
    }
    else{
      if(hypothesis == "equal-correlated"){
        C <- Pd(p - d)
        Xi <- rep(0, p - d)
        return(TestCorrelation1G(X, C, Xi, method, repetitions, seed))
      }
      if(hypothesis == "uncorrelated"){
        C <- diag(1, p - d, p - d)
        Xi <- rep(0, p - d)
        return(TestCorrelation1G(X, C, Xi, method, repetitions, seed))
      }
    }
  }
}


#' Function for the test regarding correlation matrices for multiple groups.
#'
#' @description This function conducts the test for hypotheses regarding the
#' correlation matrix for more than one group. Depending on the chosen method a
#' bootstrap, a Monte-Carlo-technique or a Taylor-base Monte-Carlo-approach is
#' used to calculate p-value of the ATS based on a specified number of runs.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix
#' @param nv vector of sample sizes for the bootstrap samples
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#'  is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1.000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @export
TestCorrelationMGinner <- function(Data, nv, C, Xi, method, repetitions = 1000,seed = 0){
  if(seed != 0){
    (set.seed(seed))
  }
  d <- dim(X)[1]
  p <- d * (d + 1) / 2
  a <- cumsum(c(1, (d):2))
  N <- sum(nv)
  kappainvv <- N / nv
  H <- matrix(rep(a, d), d, d, byrow = TRUE)
  L <- diag(1, p, p)[-a, ]
  M <- matrix(0, p, p)
  for(l in 1:p){
    M[l, c(vechp(H)[l], vechp(t(H))[l])] <- 1
  }
  M1 <- L %*% M

  Datac <- lapply(Data, centering)
  VarData <- lapply(Data, tvar)
  CorData <- lapply(VarData, cov2cor)
  vCorData <- lapply(CorData, vechp)
  DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
  HatCov <- lapply(DataQ, tvar)
  MvrH <- StUpsi <- list()

  for (i in 1:length(nv)){
    MvrH[[i]] <- (L - 1 / 2 * vCorData[[i]] * M1)
    StUpsi[[i]] <- QF(sqrt(diag(as.vector(1 / vtcrossprod(matrix(vech(VarData[[i]])[a]))), p, p)), HatCov[[i]])
  }
  Upsi <- mapply(QF, MvrH, StUpsi, SIMPLIFY = FALSE)

  MSrootUpsi <- lapply(Upsi, MSroot)
  UpsiCombined <- WDirect.sumL(Upsi, kappainvv)

  if(method == "MC"){
    ResamplingResult <- ATSwS(QF(C, UpsiCombined), repetitions)
  }
  if(method == "BT"){
    ResamplingResult <- sapply(1:repetitions, BootstrapMG, nv, N, kappainvv, C, MSrootUpsi)
  }
  if(method == "Tay"){
    P <- diag(1, p, p)[a, ]
    Q <- diag(as.vector(vech(diag(1, d, d))), p, p)
    MSrootStUpsi <- lapply(StUpsi, MSroot)
    Trace <- sum(diag(QF(C, UpsiCombined)))
    ResamplingResult <- TayappMG(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M, L, P, Q, nv)
  }

  Teststatistic <- ATS(N, unlist(vCorData), C, UpsiCombined, Xi)
  pvalue <- mean(ResamplingResult < Teststatistic)
  if (seed != 0){
    (set.seed(NULL))
  }
  return(list("pvalue" = pvalue,
              "Teststatistic" = Teststatistic,
              "CovarianceMatrix" = UpsiCombined))
}

#' Classical call of function for the test regarding correlation matrices for multiple groups.
#'
#' @description This function is for more mathematical users specifying the
#' hypothesis through hypothesis matrix C and appropriate vector Xi for which the
#' function TestCorrelationeMGinner is used. Previously in case data are given
#' as a list it is checked whether all dimensions coincide, and in case data is
#' given as a matrix it is transferred to a list.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix
#' @param nv vector of sample sizes for the bootstrap samples
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#' is used.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1.000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#'  X=matrix(rnorm(5*20),5,20)
#' Y=matrix(rnorm(5*30),5,30)
#' XY=cbind(X,Y)
#' nv=c(20,30)
#' C=Pd(2)%x%diag(1,10,10)
#' Xi=matrix(0,20,1)
#' TestCorrelationMG(XY,nv,C,Xi,method="MC")
#' @export
TestCorrelationMG <- function(X, nv, C, Xi, method, repetitions = 1000,seed = 0){
  Data <- Listcheck(X, nv)
  dimensions <- sapply(Data, dim)[1,]
  if(max(dimensions) != mean(dimensions)){
    stop("dimensions do not accord")
  }
  else{
    d <- dimensions[1]
    if(d == 1){
      stop("Correlation is only defined for dimension higher than one")
    }
    if(d > 1) {
      return(TestCorrelationMGinner(Data, nv, C, Xi, method, repetitions = 1000,seed = 0))
    }
  }
}



#' Simplified call of function for the test regarding correlation matrices for multiple groups.
#'
#' @description This function is for more applied users so no hypothesis matrix
#' or corresponding vector is necessary, and the most relevant hypothesis of
#' equal correlation matrices is checked. Appropriate C and Xi are built and the
#' function TestCorrelationMGinner is used.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix
#' @param nv vector of sample sizes for the bootstrap samples
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#'  is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1.000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' Y=matrix(rnorm(5*30),5,30)
#' XY=cbind(X,Y)
#' nv=c(20,30)
#' TestCorrelationeMGsimple(XY,nv,hypothesis="equality",method="MC")
#'
#' @export
TestCorrelationMGsimple <-function(X, nv, method = "BT", repetitions = 1000, seed = 0){
  Data <- Listcheck(X, nv)
  dimensions <- sapply(Data, dim)[1,]
  if(max(dimensions) != mean(dimensions)){
    stop("dimensions do not accord")
  }
  else{
    d <- dimensions[1]
    if(d == 1){
      stop("Correlation is only defined for dimension higher than one")
    }
    if(d > 1){
      p <- d * (d + 1) / 2
      groups <- length(nv)
      C <- Pd(groups) %x% diag(1, p - d, p - d)
      Xi <- rep(0, groups * (p - d))
      return(TestCorrelationMGinner(Data, nv, C, Xi, method, repetitions, seed))
    }
  }
}

#' Taylor-based Monte-Carlo-approximation with transformation htilde
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation with the correlation transformation
#' htilde .  After receiving some auxiliary matrices and data, the Monte-Carlo
#' observations are generated and different parts of the final sum are defined.
#' Based on this a number of the  Taylor-based ATS are calculated, where the
#' number can be chosen.
#' @param repetitions a number specifying the number of runs for the approximation
#' @param C the used hypothesis matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param Cordata the calculated correlation matrix
#' @param Jacobi the Jacobian matrix of the function hthilde applied for the
#' diagonal vectorised correlation
#' @param MvrH1 an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Atilde an auxiliary matrix for the transformation from row-wise vectorisation
#' to diagonalwise vectorisation
#' @param n1 the sample size, a scalar
#' @return a matrix containing the values of the Taylor ATS for a number of repetitions
Tayapphtilde <- function(repetitions, C, MSrootStUpsi, CorData, Jacobi, MvrH1,
                            Trace, M, L, P, Q, Atilde, n1){
  vechCorData <- vech(CorData)
  DvechCorDataM <- as.vector(vechCorData) * M
  XPB <- gData(MSrootStUpsi, repetitions)

  PUX <- P %*% XPB
  MUX0 <- M %*% Q %*% XPB

  HX <- Qvech(PUX, repetitions)

  Part1 <- MvrH1 %*% XPB
  Part2 <- 1 / (4) * as.vector(vechCorData) * HX
  Part3 <- 1 / (2) * XPB * MUX0
  Part4 <- 3 / (8) * DvechCorDataM %*% HX
  XTaydv <- Atilde %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1))
  CXTaydv <- C %*% Jacobi %*% XTaydv
  Result <- apply(CXTaydv, 2, crossprod) / Trace
  return(Result)
}



#' @title Test for the correlation matrix of data regarding their structure
#'
#' @description With this function the covariance matrix of data can be checked
#' for one of the usual structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param structure a character specifying the structure regarding them the
#' covariance matrix should be checked. Options are "Hautoregressive",
#' "diagonal", "Hcompoundsymmetry" and "Htoeplitz", where H stands for Heterogenous.
#' @param method a character, to chose whether bootstrap("BT"),
#'  Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#'  is used, while bootstrap is the predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' TestCorStructure(X,structure="Htoeplitz",method="MC")
#' @export
TestCorStructure <- function(X, structure, method = "BT",repetitions = 1000, seed = 0){
  if(seed != 0){
    (set.seed(seed))
  }
  n1 <- dim(X)[2]
  d <- dim(X)[1]
  if(d == 1){
    stop("Correlation is only defined for dimension higher than one")
  }
  if((1 - (structure %in% c("Hautoregressive","diagonal","Hcompoundsymmetry","Htoeplitz"))) == 1){
    stop("no predefined hypothesis")
  }
  else{
    if(d > 1){
      p <- d * (d + 1) / 2
      pu <- d * (d - 1) / 2
      a <- cumsum(c(1, (d):2))
      H <- matrix(rep(a, d), d, d, byrow = TRUE)
      L <- diag(1, p, p)[-a, ]
      M <- matrix(0, p, p)
      Q <- diag(as.vector(vech(diag(1, d, d))), p, p)
      for(i in 1:p){
        M[i, c(vech(t(H))[i], vech(H)[i])] <- 1
      }
      M1 <- L %*% (M + Q)
      VarData <- tvar(X)
      CorData <- cov2cor(VarData)
      vCorData <- dvechp(CorData, a, d, pu)
      Xq <- matrix(apply(X - rowMeans(X), 2, vtcrossprod), nrow = p, ncol = n1)
      HatCov <- tvar(Xq)
      MvrH1 <- (L - 1 / 2 * vechp(CorData) * M1)
      MvrH2 <- sqrt(diag(as.vector(1 / vtcrossprod(matrix(vech(VarData)[a]))), p, p))
      Atilde <- matrix(0, pu, pu)
      for (l in 1:(d - 1)){
        for (k in 1:(d - l)){
          Atilde[a[l] + k - l, a[k] - (k - l)] <- 1
        }
      }
      Upsidv <- QF(Atilde %*% MvrH1 %*% MvrH2, HatCov)

      Xi <- rep(0, pu)
      if(structure == "Hautoregressive"){
        Jacobi <- Jacobianhtilde(vCorData, a, d, p)
        Upsidvhtilde <- QF(Jacobi, Upsidv)
        C <- Pd(pu)
        Teststatistic <- ATS(n1, htilde(vCorData, a, d), C, Upsidvhtilde, Xi)
        if(method == "MC"){
          ResamplingResult <- ATSwS(QF(C, Upsidvhtilde), repetitions)
        }
        if(method == "BT"){
          ResamplingResult <- sapply(1:repetitions, Bootstraphtilde, n1, a, d,
                                     p, C, MSroot(Upsidv), vCorData)
        }
        if(method == "Tay"){
          P <- diag(1, p, p)[a, ]
          StUpsi <- QF(MvrH2, HatCov)
          Trace <- sum(diag(QF(C, Upsidvhtilde)))

          ResamplingResult <- Tayapphtilde(repetitions, C, MSroot(StUpsi), CorData, Jacobi,
                                            MvrH1, Trace, M, L, P, Q, Atilde, n1)
        }
      }


      else{
        if(structure == "diagonal"){
          C <- diag(1, pu, pu)
        }
        if(structure == "Hcompoundsymmetry"){
          C <- Pd(pu)
        }
        if(structure == "Htoeplitz"){
          C <- Pd(d - 1)
          for (l in 3:d){
            C <- direct.sum(C, Pd(d - l + 1))
          }
        }
        Teststatistic <- ATS(n1, vCorData, C, Upsidv, Xi)
        if(method == "MC"){
          ResamplingResult <- ATSwS(QF(C, Upsidv), repetitions)
        }
        if(method == "BT"){
          ResamplingResult <- sapply(1:repetitions, Bootstrap, n1, C, MSroot(Upsidv))
        }
        if(method == "Tay"){
          P <- diag(1, p, p)[a, ]
          StUpsi <- QF(MvrH2, HatCov)
          Trace <- sum(diag(QF(C, Upsidv)))
          ResamplingResult <- Tayapp(repetitions, C, MSroot(StUpsi), CorData,
                                    MvrH1, Trace, M, L, P, Q, Atilde,n1)
        }
      }
      pvalue <- mean(ResamplingResult < Teststatistic)
      if(seed != 0){
        set.seed(NULL)
      }
      return(list("pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsidv))
    }
  }
}
