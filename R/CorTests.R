#' @title Base function for testing correlation matrices
#'
#' @description This function conducts the test for hypotheses regarding the
#' correlation matrix. Depending on the chosen method a
#' bootstrap or Monte-Carlo-technique is used to calculate p-value of the ATS
#' based on a specified number of runs.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix and nv is used to indicate
#' the group sizes. For one group, nv is not necessary
#' @param nv vector of sample sizes for the bootstrap samples per group
#' @param C hypothesis matrix for calculating the ATS
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") or Taylor-based Monte-Carlos-approach("Tay") is used.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @param hypothesis character or NULL, will be displayed in the print call
#' @return an object of the class 'CovTest'
#'
#' @export
TestCorrelation_base <- function(X, nv = NULL, C, Xi, method, repetitions = 1000,
                                 seed = NULL, hypothesis = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }

  # one group
  if((length(nv) == 1) | is.null(nv)){
    nv <- dim(X)[2]
  }

  d <- dim(X)[1]
  p <- d*(d+1)/2
  a <- cumsum(c(1,(d):2))
  N <- sum(nv)
  H <- matrix(rep(a,d), d, d, byrow=TRUE)
  L <- diag(1,p,p)[-a,]
  M <- matrix(0,p,p)

  for(i in 1:p){
    M[i,c(vechp(H)[i],vechp(t(H))[i])] <- 1
  }
  M1 <- L%*%M

  # one group
  if(length(nv) == 1){
    VarData <- stats::var(t(X))
    CorData <- stats::cov2cor(VarData)
    vCorData <- vechp(CorData)
    Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod),nrow=p,ncol=nv)
    HatCov <- stats::var(t(Xq))
    MvrH1 <- (L-1/2*vCorData*M1)
    MvrH2 <- sqrt(diag(as.vector(1/vtcrossprod(matrix(matrixcalc::vech(VarData)[a]))),p,p))
    Upsi <- QF(MvrH1%*%MvrH2,HatCov)
    MSrootUpsi <- MSroot(Upsi)
    StUpsi <- QF(MvrH2,HatCov)
    MSrootStUpsi <- MSroot(StUpsi)
  }
  # multiple groups
  else{
    Datac <- lapply(X, centering)
    VarData <- lapply(X, function(X) stats::var(t(X)))
    CorData <- lapply(VarData, stats::cov2cor)
    vCorData <- lapply(CorData, vechp)
    DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
    HatCov <- lapply(DataQ, function(X) stats::var(t(X)))
    MvrH <- StUpsi <- list()

    for (i in 1:length(nv)){
      MvrH[[i]] <- (L - 1 / 2 * vCorData[[i]] * M1)
      StUpsi[[i]] <- QF(sqrt(diag(as.vector(1 / vtcrossprod(matrix(matrixcalc::vech(VarData[[i]])[a]))), p, p)), HatCov[[i]])
    }
    Upsi_list <- mapply(QF, MvrH, StUpsi, SIMPLIFY = FALSE)
    MSrootStUpsi <- lapply(StUpsi, MSroot)
    MSrootUpsi <- lapply(Upsi_list, MSroot)
    kappainvv <- N / nv
    Upsi <- WDirect.sumL(Upsi_list, kappainvv)
  }

  if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsi), repetitions) }
  if(method == "BT"){ ResamplingResult <- sapply(1:repetitions, Bootstrap, nv, C, MSrootUpsi) }
  if(method == "Tay"){
    P <- diag(1,p,p)[a,]
    Q <- diag(as.vector(matrixcalc::vech(diag(1,d,d))),p,p)
    Trace <- sum(diag(QF(C,Upsi)))
    if(length(nv) == 1){ # one group
      ResamplingResult <- Tayapp1G(repetitions, C, MSrootStUpsi, CorData, MvrH1, Trace, M, L, P, Q, nv, NULL)
    }
    else{ #multiple groups
      ResamplingResult <- TayappMG(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M, L, P, Q, nv)
    }
  }


  Teststatistic <- ATS(N, unlist(vCorData), C, Upsi, Xi)
  pvalue <- mean(ResamplingResult < Teststatistic)

  CovTest <- list("method" = "Correlation",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsi,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = hypothesis,
                  "nv" = nv)

  class(CovTest) <- "CovTest"

  return(CovTest)
}


#' @title Simplified call for test regarding correlation matrices
#'
#' @description This function is for more applied users so no hypothesis matrix or
#' corresponding vector is necessary. This is replaced by predefined hypotheses,
#' from which is chosen. From this C and Xi are built and the function
#'  \code{\link{TestCorrelation_base}} is used.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix and nv is used to indicate
#' the group sizes. For one group, nv is not necessary
#' @param nv vector of sample sizes for the bootstrap samples per group
#' @param hypothesis a character to choose one of the predefined hypotheses which are
#' "equal-correlated" and "uncorrelated" (only possible for one group)
#' @param method a character, to chose whether bootstrap("BT") or Taylor-based Monte-Carlos-approach("Tay")
#' Monte-Carlo-technique("MC") or  is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return an object of the class 'CovTest'
#'
#'
#' @export
TestCorrelation_simple <- function(X, nv = NULL, hypothesis, method = "BT", repetitions = 100, seed = NULL){
  # one group
  if(is.null(nv) | length(nv) == 1){
    d <- dim(X)[1]
  }
  # multiple groups
  else{
    X <- Listcheck(X, nv)
    dimensions <- sapply(X, dim)[1,]
    if(max(dimensions) != mean(dimensions)){
      stop("dimensions do not accord")
    }
    d <- dimensions[1]
    groups <- length(nv)
  }

  if(d == 1){
    stop("Correlation is only defined for dimensions higher than one")
  }
  p <- d * (d+1) / 2
  if(!(hypothesis %in% c("equal-correlated", "uncorrelated"))){
    stop("no predefined hypothesis")
  }

  if(hypothesis == "equal-correlated"){
    if(is.null(nv) | length(nv) == 1){
      C <- Pd(p-d)
      Xi <- rep(0, p-d)
    }
    else{
      C <- Pd(groups) %x% diag(1, p - d, p - d)
      Xi <- rep(0, groups * (p - d))
    }
    return(TestCorrelation_base(X, nv, C, Xi, method, repetitions = 1000,
                                 seed , hypothesis))
  }
  if(hypothesis == "uncorrelated"){
    if(!is.null(nv) | length(nv) != 1){
      stop("the hypothesis 'uncorrelated' can only be tested for a single group")
    }
    C <- diag(1, p - d, p - d)
    Xi <- rep(0, p - d)
    return(TestCorrelation_base(X, nv, C, Xi, method, repetitions, seed, hypothesis))
  }

}



#' @title Test for structure of data's correlation matrix
#'
#' @description With this function the correlation matrix of data can be checked
#' for one of the usual structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param structure a character specifying the structure regarding them the
#' covariance matrix should be checked. Options are "Hautoregressive", "diagonal",
#' "Hcompoundsymmetry" and "Htoeplitz".
#' @param method a character, to chose whether bootstrap("BT") or Taylor-based Monte-Carlos-approach("Tay")
#' Monte-Carlo-technique("MC") or is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set. A chosen seed is deleted at the end.
#' @return an object of the class 'CovTest'
TestCorrelation_structure <- function(X, structure, method = "BT", repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }

  n1 <- dim(X)[2]
  d <- dim(X)[1]

  if(d == 1){
    stop("Correlation is only defined for dimension higher than one")
  }

  if(!(structure %in% c("Hautoregressive","diagonal","Hcompoundsymmetry","Htoeplitz"))){
    stop("no predefined hypothesis")
  }

  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))
  H <- matrix(rep(a, d), d, d)
  L <- diag(1, p, p)[-a, ]
  M <- matrix(0, p, p)
  Q <- diag(as.vector(matrixcalc::vech(diag(1, d, d))), p, p)
  for(i in 1:p){
    M[i, c(matrixcalc::vech(t(H))[i], matrixcalc::vech(H)[i])] <- 1
  }
  M1 <- L %*% (M + Q)
  VarData <- stats::var(t(X))
  CorData <- stats::cov2cor(VarData)
  vCorData <- dvech(CorData, a, d, pu, inc_diag = FALSE)
  Xq <- matrix(apply(X - rowMeans(X), 2, vtcrossprod), nrow = p, ncol = n1)
  HatCov <- stats::var(t(Xq))
  MvrH1 <- (L - 1 / 2 * vechp(CorData) * M1)

  MvrH2 <- sqrt(diag(as.vector(1 / vtcrossprod(matrix(matrixcalc::vech(VarData)[a]))), p, p))
  Atilde <- matrix(0, pu, pu)
  for(l in 1:(d - 1)){
    for(k in 1:(d - l)){
      Atilde[a[l] + k - l, a[k] - (k - l)] <- 1
    }
  }

  Upsidv <- QF(Atilde %*% MvrH1 %*% MvrH2, HatCov)
  Xi <- rep(0, pu)

  if(structure == "Hautoregressive"){
    Jacobi <-   Jacobian(vCorData, a, d, p, fun = "ascending_root_fct_cor")
    Upsidvhtilde <- QF(Jacobi, Upsidv)
    C <- Pd(pu)
    Teststatistic <- ATS(n1, ascending_root_fct_cor(vCorData, a, d), C, Upsidvhtilde, Xi)
    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidvhtilde), repetitions) }
    if(method == "BT"){ ResamplingResult <- sapply(1:repetitions, Bootstrap_trans, n1, a, d,
                                                   p, C, MSroot(Upsidv), vCorData, fun = "ascending_root_fct_cor") }
    if(method == "Tay"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidvhtilde)))

      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                       MvrH1, Trace, M, L, P, Q, n1, Atilde, Jacobi)
    }
  }
  else{
    if(structure == "diagnoal"){
      C <- diag(1, pu, pu)
    }
    if(structure == "Hcompoundsymmetry"){
      C <- Pd(pu)
    }
    if(structure == "Htoeplitz"){
      C <- Pd(d - 1)
      for (l in 3:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
    }
    Teststatistic <- ATS(n1, vCorData, C, Upsidv, Xi)

    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidv), repetitions) }
    if(method == "BT"){ ResamplingResults <- sapply(1:repetitions, Bootstrap, n1, C, MSroot(Upsidv)) }
    if(method == "Tay"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidv)))
      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                 MvrH1, Trace, M, L, P, Q, Atilde,n1)
    }
  }

  pvalue <- mean(ResamplingResult < Teststatistic)


  CovTest <- list("method" = "Correlation",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsidv,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = structure,
                  "nv" = n1)

  class(CovTest) <- "CovTest"

  return(CovTest)


}
