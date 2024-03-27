#' @title Root transformation of the vectorised covariance matrix
#'
#' @description A function calculating the roots of a vectorised covariance matrix.
#' The roots increasing, so square root for the first secondary diagonals, third
#' root for the second secondary diagonal and so on. For roots with even order the
#' absolute value of the argument is used, since the arguments can be negative.
#'
#' @param x vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#' @export
h <- function(x, a, d){
  for(i in 3:d){
    if(i %% 2 == 1){
      x[(0:(d - i)) + a[i]] <- abs( x[(0:(d - i)) + a[i]] )^(1/(i - 1))
    }
    if((i %% 2 == 0)){
      x[(0:(d - i)) + a[i]] <- (x[(0:(d - i)) + a[i]] <= 0) * (-abs(x[(0:(d - i)) + a[i]])^(1/(i - 1))) +
                                  (x[(0:(d - i)) + a[i]] > 0) * (abs(x[(0:(d - i)) + a[i]])^(1/(i - 1)))
    }
  }
  return(x)
}


#' @title  Transformation of the vectorised covariance matrix by quotients of means
#'
#' @description A function which calculates the mean of the secondary diagonals
#' and divide them through the next one. Since the elements can be negative, for
#' the denominator absolute values are used.
#' @param v vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#' @export
g <- function(v, a, d){
  ratio <- rep(0, d - 1)
  ae <-  c(a, a[d] + 1)
  for(l in 2:(d)){
    ratio[l - 1] <-  mean(v[ae[l]:(ae[l + 1] - 1)]) / mean(abs(v[ae[l - 1]:(ae[l] - 1)]))
  }
  return(c(v, ratio))
}

#' @title Jacobian matrix for transformation functions
#'
#' @description A function which calculates the Jacobian matrix for a given
#' transformation function h or g
#' @param X vectorised covariance matrix for which the Jacobian matrix is applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param fun transformation function, that should be used. 'g' or 'h'
#' @return the Jacobian matrix applied for the given vector
#' @export
Jacobian <- function(X, a, d, p, fun){
  if(fun == "h"){
    E <- rep(1, p)
    for(i in 3:d){
      if(i %% 2 == 1){
        E[(0:(d - i)) + a[i]] <-  X[(0:(d - i)) + a[i]] / ((i - 1) * abs(X[(0:(d - i)) + a[i]])^(2 - 1/(i - 1)))
      }
      if(i %% 2 == 0){
        E[(0:(d - i)) + a[i]] <-  1/((i - 1) * abs(X[(0:(d - i)) + a[i]])^(1 - (1/(i - 1))))
      }
    }
    return(diag(E, p, p))
  }
  else{
    if(fun == "g"){
      J <-  matrix(0, d - 1, p)
      for (l in 1:(d - 1)){
        S1 <- sum(abs(X[a[l]:(a[l] + d - l)]))
        S2 <- sum(X[a[l + 1]:(a[l + 1] + d - l - 1)])

        J[l, a[l + 1] + 0:(d - l - 1)] <- (d - l + 1) / (d - l) / S1
        J[l, a[l] + 0:(d - l)] <- (d - l + 1) / (d - l) * sign(X[a[l] + 0:(d - l)]) * (-S2 / (S1) ^ 2)
      }
      return(rbind(diag(1, p, p), J))
    }
    else{
      stop("fun must be 'g' or 'h'")
    }
  }
}

#' @title ATS for transformed vectors
#'
#' @description A function which calculates the Anova type statistic based on
#' a transformation function
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param v vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belongs to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param fun transformation function, that should be used. 'g' or 'h'
#' @return a scalar, the value of the ATS
#' @export
ATS_fun <- function(N, X, C, v, a, d, p, fun){
  Xmean <-  rowMeans(X)
  CDiff <- C %*% (do.call(fun, list(Xmean, a, d)) - do.call(fun, list(v, a, d)))
  Jacobi <-  Jacobian(Xmean, a, d, p, fun)
  HatCov <-  tvar(X)
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
#' @param p dimension of the vectorised matrix
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix root of the covariance matrix HatCov, to generate
#' the bootstrap sample
#' @param vX the expectation vector for the bootstrap sample
#' @param fun transformation function, that should be used. 'g' or 'h'
#' @return a scalar, the value of the ATS
#' @export
Bootstrap_trans <- function(N.sim, n1, a, d, p, C, MSrootHatCov, vX, fun){
  XPB <- gData(MSrootHatCov, n1) + vX
  return(ATS_fun(n1, XPB, C, vX, a, d, p, fun))
}


#' @title Bootstrap for one and multiple groups
#'
#' @description This function generates normal distributed random vectors. For one group,
#' nv random vectors with covariance matrix HatCov are generated and the corresponding value
#' of the ATS is generated. For multiple groups the corresponding sample sizes from nv are used.
#' The weighted sum of covariance matrices is calculated and used to calculate the value
#' of the ATS.
#' @param N.sim control variable for using sapply
#' @param nv scalar (one group) or vector (multiple groups) of sample sizes for the bootstrap samples
#' @param C hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix (one group) or list of matrices (multiple groups) of roots of the covariance matrices, to generate
#' the bootstrap sample
#' @return a scalar, the value of the ATS
#' @export
Bootstrap <- function(N.sim, nv, C, MSrootHatCov){
  # one group
  if(length(nv) == 1){
    XPB <- gData(MSrootHatCov, nv)
    PBHatCov <- tvar(XPB)
    return(ATS(nv, rowMeans(XPB), C, PBHatCov))
  }
  # multiple groups
  else{
    N <- sum(nv)
    kappainvv <- N / nv

    DataPB <- mapply(gData, MSrootHatCov, nv, SIMPLIFY = FALSE)
    PBHatCov <- WDirect.sumL(lapply(DataPB, tvar), kappainvv)
    return(ATS(N, unlist(lapply(DataPB, rowMeans)), C, PBHatCov))

  }
}


#' @title Inner function for test regarding covariance matrices with one group.
#'
#' @description This function conducts the test for hypotheses regarding the
#' covariance matrix for one group. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X a matrix containing the observation vectors as columns
#' @param C hypothesis matrix for calculating the ATS
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or Monte-Carlo-technique("MC")
#' is used.
#' @param repetitions a scalar, indicates the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X <- matrix(rnorm(5*20),5,20)
#' d <- dim(X)[1]
#' p <- d*(d+1)/2
#' a <- 1
#' C <- diag(1,p,p)[(1:d),]
#' Xi <- rep(0,d)
#' TestCovariance1G(X,C,Xi,method="MC")
#' @export
TestCovariance1G <- function(X, C, Xi, method, repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  n1 <- dim(X)[2]
  vX <- matrixcalc::vech(tvar(X))
  Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod),ncol=n1)
  HatCov <- tvar(Xq)
  if(method=="MC"){  ResamplingResult <- ATSwS(QF(C,HatCov),repetitions) }
  if(method=="BT"){  ResamplingResult <- sapply(1:repetitions,Bootstrap,n1,C,MSroot(HatCov)) }

  Teststatistic <- ATS(n1, vX, C, HatCov, Xi)
  pvalue <- mean(ResamplingResult < Teststatistic)

  return(list("pvalue"=pvalue, "Teststatistic"=Teststatistic, "CovarianceMatrix"=HatCov))
}


#' @title Inner function for test regarding covariance matrices for multiple groups.
#'
#' @description This function conducts the test for hypotheses regarding the
#' covariance matrix for more than one group. Depending on the chosen method a
#' bootstrap  or Monte-Carlo-technique is used to calculate p-value of the ATS
#' based on a specified number of runs.
#' @param X a list containing the observation vectors. Each matrix in this list
#' is another group, where the observation vectors are the columns
#' @param nv vector of sample sizes for the bootstrap samples
#' @param C hypothesis matrix for calculating the ATS
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @export
TestCovarianceMGinner <- function(X, nv, C, Xi, method, repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  N <- sum(nv)
  kappainvv <- N / nv

  Datac <- lapply(X, centering) #Centered random variables
  VarData <- lapply(X, tvar)
  vVarData <-  unlist(lapply(VarData, matrixcalc::vech))
  DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
  HatCov <- lapply(DataQ, tvar)

  MSrootHatCov <- lapply(HatCov, MSroot)
  HatCovCombined <- WDirect.sumL(HatCov, kappainvv)

  if(method == "MC"){
    ResamplingResult <- ATSwS(QF(C, HatCovCombined), repetitions)
  }
  if(method == "BT"){
    ResamplingResult <- sapply(1:repetitions, Bootstrap, nv, C, MSrootHatCov)
  }

  Teststatistic <- ATS(N, vVarData, C, HatCovCombined, Xi)
  pvalue <- mean(ResamplingResult < Teststatistic)
  return(list("pvalue" = pvalue,"Teststatistic" = Teststatistic,"CovarianceMatrix" = HatCov))

}


#' @title Simplified call for test regarding covariance matrices with one group.
#'
#' @description This function is for more applied users so no hypothesis matrix or
#' corresponding vector is necessary. This is replaced by predefined hypotheses,
#' from which is chosen. From this C and Xi are built and the function
#'  \code{\link{Testcovariance1G}} is used.
#' @param X a matrix containing the observation vectors as columns
#' @param hypothesis a character to choose one of the predefined hypotheses which are
#' "equal-variances", "uncorrelated", "given-trace" and "given-matrix".
#' @param Xi a scalar or square matrix to specify the hypothesis in case of
#' "given-trace" or "given-matrix". The value is predefined as 0, which means
#' trace of 1 resp the identity as given matrix.
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' #' @examples
#' X <- matrix(rnorm(5*20),5,20)
#' d <- dim(X)[1]
#' p <- d*(d+1)/2
#' a <- 1
#' C <- diag(1,p,p)[a ,]
#' Xi <- rep(0,d)
#' TestCovariance1Gsimple(X, hypothesis="equal-variances")
#' @export
TestCovariance1Gsimple <- function(X, hypothesis, Xi = 0, method = "BT",
                                    repetitions = 1000, seed = NULL){
  d <-  dim(X)[1]
  p <-  d * (d + 1) / 2 #dimension vectorized covariance matrix
  if (d > 1){
    a <- cumsum(c(1, (d):2))
  }
  else{
    a <- 1
  }
  if( !(hypothesis %in% c("equal-variances","uncorrelated","given-trace","given-matrix"))){
    stop("no predefined hypothesis")
  }
  else{
    if (hypothesis == "equal-variances"){
      C <-  diag(1, p, p)[a,]
      Xi <-  rep(0, d)
      return(TestCovariance1G(X, C, Xi, method, repetitions, seed))
    }
    if(hypothesis == "given-trace"){
      Tracevector <-  matrix(0, 1, p)
      Tracevector[1, a] <-  1
      C <- Tracevector
      if(Xi == 0){
        Xi <-  1
      }
      return(TestCovariance1G(X, C, Xi, method, repetitions, seed))
    }
    if(hypothesis == "uncorrelated"){
      C <-  diag(1, p, p)[-a,]
      Xi <-  rep(0, p - d)
      return(TestCovariance1G(X, C, Xi, method, repetitions, seed))
    }
    if (hypothesis == "given-matrix"){
      C <-  diag(1, p, p)
      CHECK <-  0
      if (max(abs(Xi)) == 0 &
          length(Xi) == 1) {
        CHECK <-  1
        Xi <- matrixcalc::vech(diag(1, d, d))
      }
      if(max(abs(Xi)) != 0 & length(Xi) == 1 & d == 1){
        CHECK <- 1
      }
      if(is.matrix(Xi)){
        if(matrixcalc::is.square.matrix(Xi) & dim(Xi)[1] == d){
          CHECK <-  1
          Xi <- matrixcalc::vech(Xi)
        }
      }

      if (CHECK == 0) {
        stop("incorrect given matrix")
      }
      else{
        return(TestCovariance1G(X, C, Xi, method, repetitions, seed))
      }
    }
  }
}


#' @title Simplified call for test regarding covariance matrices for multiple groups.
#'
#' @description This function is for more applied users so no hypothesis matrix
#' or corresponding vector is necessary. This is replaced by predefined hypotheses,
#' from which is chosen. From this C and Xi are built and the function
#' TestcovarianceMGinner is used. Previously in case data are given as a list it
#' is checked whether all dimensions coincide,#' and in case data is given as
#' matrix it is transferred to a list.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix
#' @param nv vector of sample sizes for the bootstrap samples
#' @param hypothesis a character to choose one of the predefined hypotheses which are
#' "equal-variances", "uncorrelated", "given-trace" and "given-matrix".
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' Y=matrix(rnorm(5*30),5,30)
#' XY=cbind(X,Y)
#' nv=c(20,30)
#' TestCovarianceMGsimple(XY,nv,hypothesis="equal",method="MC")
#' @export
TestCovarianceMGsimple <- function(X, nv, hypothesis, method = "BT",
                                    repetitions = 1000, seed = NULL){
  Data <-  Listcheck(X, nv)
  dimensions <-  sapply(Data, dim)[1,]
  if(max(dimensions) != mean(dimensions)){
    stop("dimensions do not accord")
  }
  else{
    d <-  dimensions[1]
    p <-  d * (d + 1) / 2#dimension vectorized covariance matrix
    if(d > 1){
      a <-  cumsum(c(1, (d):2))
    }
    else{
      a <-1
    }
    groups <- length(nv)
    if(!(hypothesis %in% c("equal", "equal-trace", "equal-diagonals"))){
      stop("no predefined hypothesis")
    }
    else{
      if (hypothesis == "equal"){
        C <- Pd(groups) %x% diag(1, p, p)
        Xi <- rep(0, p * groups)
        return(TestCovarianceMGinner(Data, nv, C, Xi, method, repetitions))
      }
      if(hypothesis == "equal-trace"){
        Tracevector <- matrix(0, 1, p)
        Tracevector[1, a] <- 1
        C <- Pd(groups) %x% Tracevector
        Xi <- rep(0, groups)
        return(TestCovarianceMGinner(Data, nv, C, Xi, method, repetitions))
      }
      if(hypothesis == "equal-diagonals"){
        C <- Pd(groups) %x% diag(1, p, p)[a,]
        Xi <- rep(0, times = groups * d)
        return(TestCovarianceMGinner(Data, nv, C, Xi, method, repetitions, seed))
      }
    }
  }
}



#' @title Classical call for test regarding covariance matrices for multiple groups.
#'
#' @description This function is for more mathematical users specifying the hypothesis
#' through hypothesis matrix C and appropriate vector Xi for which the function
#' TestCovarianceMGinner is used. Previously in case data are given as a list it
#' is checked whether all dimensions coincide, and in case data is given as matrix
#' it is transferred to a list.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix
#' @param nv vector of sample sizes for the bootstrap samples
#' @param C hypothesis matrix for calculating the ATS
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#'
#' X <- matrix(rnorm(5*20),5,20)
#' Y <- matrix(rnorm(5*30),5,30)
#' XY <- cbind(X,Y)
#' nv <- c(20,30)
#' C <- Pd(2)%x%diag(1,15,15)
#' Xi <- matrix(0,30,1)
#' TestCovarianceMG(XY,nv,C,Xi,method="MC")
#' @export
TestCovarianceMG <- function(X, nv, C, Xi, method = "BT", repetitions = 1000, seed = NULL){
  Data <-  Listcheck(X, nv)
  dimensions <-  sapply(Data, dim)[1,]
  if(max(dimensions) != mean(dimensions)){
    stop("dimensions do not accord")
  }
  else{
    return(TestCovarianceMGinner(Data, nv, C, Xi, method, repetitions, seed))
  }
}




#' @title Test for structure of data's covariance matrix
#'
#' @description With this function the covariance matrix of data can be checked
#' for one of the usual structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param structure a character specifying the structure regarding them the
#' covariance matrix should be checked. Options are "autoregressive", "FO-autoregressive"
#' "diagonal", "sphericity", "compoundsymmetry" and "toeplitz".
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' TestCovStructure(X,structure="toeplitz",method="MC")
#' @export
TestCovStructure <- function(X, structure, method, repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  n1 <- dim(X)[2]
  d <- dim(X)[1]
  if(d==1){ stop("Structures can be only investigated for more than one dimension") }
  if(!(structure %in% c("autoregressive", "FO-autoregressive", "diagonal", "sphericity", "compoundsymmetry", "toeplitz") )){
    stop("no predefined hypothesis")
  }
  else{
    if(d>1){
      p <- d*(d+1)/2 #dimension vectorized covariance matrix
      a <- cumsum(c(1,(d):2))

      vX <- dvech(tvar(X),a,d,p, inc_diag = TRUE)
      Xq <- apply(X-rowMeans(X),2,vdtcrossprod,a,d,p)
      HatCov <- tvar(Xq)
      if(structure %in% c("autoregressive","FO-autoregressive")){
        if(structure == "autoregressive"){
          C <- matrixcalc::direct.sum(diag(1,d,d),Pd(p-d))
          Xi <- c(rep(1,times=d),rep(0,times=p-d))
          Jacobi <- Jacobian(vX,a,d,p, 'h')
          HatCovh <- QF(Jacobi,HatCov)

          if(method == "MC"){ ResamplingResult=ATSwS(QF(C,HatCovh),repetitions) }
          if(method == "BT"){
            ResamplingResult <- sapply(1:repetitions,Bootstrap_trans,n1,a,d,p,C,MSroot(HatCov),vX,"h")
          }
          Teststatistic <- ATS(n1,h(vX,a,d),C,HatCovh,Xi)
          pvalue <- mean(ResamplingResult < Teststatistic)
          return(list("pvalue"=pvalue, "Teststatistic"=Teststatistic, "CovarianceMatrix"=HatCov))
        }
        if(structure=="FO-autoregressive"){
          C <- Pd(d)
          for(l in 2:d){
            C <- matrixcalc::direct.sum(C,Pd(d-l+1))
            }
          C <- matrixcalc::direct.sum(C,Pd(d-1))
          Xi <- rep(0,times=p+d-1)
          Jacobi <- Jacobian(vX,a,d,p, 'g')
          HatCovg <- QF(Jacobi,HatCov)

          if(method=="MC"){ ResamplingResult <- ATSwS(QF(C,HatCovg),repetitions) }
          if(method=="BT"){
            ResamplingResult=sapply(1:repetitions,Bootstrap_trans,n1,a,d,p,C,MSroot(HatCov),vX,"g")
          }
          Teststatistic <- ATS(n1,g(vX,a,d),C,HatCovg,Xi)
            pvalue <- mean(ResamplingResult<Teststatistic)
          return(list("pvalue"=pvalue, "Teststatistic"=Teststatistic, "CovarianceMatrix"=HatCov))
        }
      }
      else{
          Xi <- rep(0,p)
        if(structure=="diagonal"){
          C <- matrixcalc::direct.sum(matrix(0,d,d),diag(1,p-d,p-d))
        }
          if(structure=="sphericity"){
          C <- matrixcalc::direct.sum(Pd(d),diag(1,p-d,p-d))
        }
        if(structure=="compoundsymmetry"){
          C <- matrixcalc::direct.sum(Pd(d),Pd(p-d))
          }
        if(structure=="toeplitz"){
          C <- Pd(d)
          for(l in 2:d){
              C <- matrixcalc::direct.sum(C,Pd(d-l+1))
          }
        }
        if(method=="MC"){
          ResamplingResult <- ATSwS(QF(C,HatCov),repetitions)
          }
        if(method=="BT"){
          ResamplingResult <- sapply(1:repetitions,Bootstrap,n1,C,MSroot(HatCov))
          }

        Teststatistic <- ATS(n1,vX,C,HatCov,Xi)
        pvalue <- mean(ResamplingResult<Teststatistic)

        return(list("pvalue"=pvalue, "Teststatistic"=Teststatistic, "CovarianceMatrix"=HatCov))
      }
    }
  }
}

