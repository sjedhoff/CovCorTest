#' @title Base function for testing covariance matrices
#'
#' @description This function conducts the test for hypotheses regarding the
#' covariance matrix. Depending on the chosen method a
#' bootstrap  or Monte-Carlo-technique is used to calculate p-value of the Anova-Type-Statistic
#' based on a specified number of runs.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix and nv is used to indicate
#' the group sizes. For one group, nv is not necessary
#' @param nv vector of group sizes
#' @param C hypothesis matrix for calculating the Anova-Type-Statistic
#' @param Xi a vector defining together with C the investigated hypothesis
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used.
#' @param repetitions a scalar, indicates the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set.
#' @param hypothesis character or NULL, will be displayed in the print call
#' @return an object of the class \code{\link{CovTest}}
#'
#' @references \insertRef{sattler_cov_2020}{CovCorTest}
#' @export
#'
#' @import MANOVA.RM
#' @importFrom Rdpack reprompt
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' vars <- colnames(EEGwide)[1:6]
#'
#' X <- t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "AD",vars])
#'
#' # Testing the trace
#' C <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1), nrow = 1, ncol = 21)
#' Xi <- 2
#'
#' TestCovariance_base(X = X, nv = NULL, C = C, Xi = Xi, method = "BT", repetitions = 1000,
#'       seed = 31415, hypothesis = "Trace")
#'
TestCovariance_base <- function(X, nv = NULL, C, Xi, method, repetitions = 1000,
                                seed = NULL, hypothesis = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }

  method <- toupper(method)
  if(!(method == "MC" | method == "BT")){
    stop("method must be bootstrap ('BT') or Monte-Carlo-technique('MC')")
  }

  listcheck <- Listcheck(X,nv)
  X <- listcheck[[1]]
  nv <- listcheck[[2]]

  # just one group is nv = NULL
  if(is.null(nv)){
    dimensions <- dim(X)[1]
    groups <- 1


    nv <- dim(X)[2]
    vX <- matrixcalc::vech(stats::var(t(X)))
    Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod),ncol=nv)
    HatCov <- stats::var(t(Xq))
    MSrootHatCov <- MSroot(HatCov)


  }
  # multiple groups
  else{
    dimensions <- unlist(lapply(X, nrow))
    groups <- length(nv)


    N  <- sum(nv)
    kappainvv <- N / nv

    Datac <- lapply(X, function(x) x-rowMeans(x))
    VarData <- lapply(X, function(X) stats::var(t(X)))
    vX <-  unlist(lapply(VarData, matrixcalc::vech))
    DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
    HatCov_list <- lapply(DataQ, function(X) stats::var(t(X)))

    MSrootHatCov <- lapply(HatCov_list, MSroot)
    HatCov <- WDirect.sumL(HatCov_list, kappainvv)
  }

  # Check dimensions of C and Xi
  d <- dimensions[1]
  p <- d * (d+1) / 2
  ifelse(d > 1, a <- cumsum(c(1, d:2)), a <- 1)
  if(!is.matrix(C)){
    stop("C must be a matrix")
  }
  if( (nrow(C) != length(Xi)) | (ncol(C) != groups*p) ){
    stop("dimensions of C and Xi do not align")
  }

  if(method == "MC"){
    ResamplingResult <- ATSwS(QF(C, HatCov), repetitions)
  }
  else if(method == "BT"){
    ResamplingResult <- sapply(1:repetitions, Bootstrap, nv, C, MSrootHatCov)
  }

  Teststatistic <- ATS(sum(nv), vX, C, HatCov, Xi)
  pvalue <- mean(ResamplingResult > Teststatistic)

  if(is.null(hypothesis)){
    hypothesis <- "C * v = Xi"
  }

  CovTest <- list("method" = "Covariance",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = HatCov,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = hypothesis,
                  "nv" = nv)

  class(CovTest) <- "CovTest"

  return(CovTest)

}



#' @title Simplified call for test regarding covariance matrices
#'
#' @description This function conducts tests for hypotheses regarding the covariance
#' matrix for a more applied user. A hypothesis can be selected out of a group of
#' predefined hypotheses. From this C and Xi are built and the function
#'  \code{\link{TestCovariance_base}} is used.
#' @param X a list or matrix containing the observation vectors. In case of a list,
#' each matrix in this list is another group, where the observation vectors are the
#' columns. For a matrix, all groups are together in one matrix and nv is used to indicate
#' the group sizes. For one group, nv is not necessary
#' @param nv vector of group sizes
#' @param hypothesis a character to choose one of the predefined hypotheses which are
#' "equal", for multiple groups "equal-trace" and "equal-diagonals" and for a single group
#' "given-trace", "given-matrix" and "uncorrelated"
#' @param A a scalar or square matrix to specify the hypothesis in case of
#' "given-trace" or "given-matrix". The value is predefined as NULL, which means
#' trace of 1 respectively the identity as given matrix.
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set.
#' @return an object of the class \code{\link{CovTest}}
#'
#' @references \insertRef{sattler_cov_2020}{CovCorTest}
#'
#' @import MANOVA.RM
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' vars <- colnames(EEGwide)[1:6]
#'
#' # Part the data into six groups of sex and diagnosis
#' X_list <- list(t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "AD",vars]),
#'                t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "MCI",vars]),
#'                t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "SCC",vars]),
#'                t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",vars]),
#'                t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "MCI",vars]),
#'                t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "SCC",vars]))
#'
#' nv <- unlist(lapply(X_list, ncol))
#'
#' TestCovariance_simple(X = X_list, nv = nv, hypothesis = "equal-trace", method = "MC",
#'                      repetitions = 1000, seed = NULL)
#'
#' TestCovariance_simple(X_list[[1]], hypothesis = "given-trace", A = 3)
#'
#' @export
TestCovariance_simple <- function(X, nv = NULL, hypothesis, A = NULL, method = "MC",
                                  repetitions = 1000, seed = NULL){
  hypothesis <- tolower(hypothesis)
  method <- toupper(method)
  if(!(method == "MC" | method == "BT")){
    stop("method must be bootstrap ('BT') or Monte-Carlo-technique('MC')")
  }

  listcheck <- Listcheck(X,nv)
  X <- listcheck[[1]]
  nv <- listcheck[[2]]

  # multiple groups
  if(!is.null(nv)){
    dimensions <- unlist(lapply(X, nrow))
    groups <- length(nv)
  }
  # one group
  else{
    dimensions <- dim(X)[1]
    groups <- 1
  }

  d <- dimensions[1]
  p <- d * (d+1) / 2
  ifelse(d > 1, a <- cumsum(c(1, d:2)), a <- 1)

  if(!(hypothesis %in%
       c("equal", "equal-trace", "equal-diagonals", "given-trace", "given-matrix", "uncorrelated"))){
    stop("no predefined hypothesis")
  }
  if(!is.null(A) & !(hypothesis %in% c("given-trace", "given-matrix"))){
    warning(paste0("the input argument A is not used, since the selected hypothesis is '", hypothesis, "'"))
  }

  if(hypothesis == "equal"){
    if(groups == 1){
      C <- Pd(d) %*% diag(1, p, p)[a,]
      Xi <- rep(0, d)
    }
    else{
      C <- Pd(groups) %x% diag(1, p, p)
      Xi <- rep(0, p*groups)
    }
    return(TestCovariance_base(X, nv = nv, C = C, Xi = Xi, method = method,
                               repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }
  if(hypothesis == "given-trace"){
    if(groups > 1){
      stop("the hypothesis 'given-trace' can only be tested for one group")
    }
    tracevec <- matrix(0, 1, p)
    tracevec[1,a] <- 1
    C <- tracevec
    if(is.null(A)){
      A <- 1
      warning("since no input A for a trace to be tested is given, a trace of 1 is tested")
    }
    if(!is.numeric(A) | length(A) != 1){
      stop("for testing the trace the input A must be a scalar")
    }
    return(TestCovariance_base(X, nv = nv, C = C, Xi = A, method = method,
                               repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }
  if(hypothesis == "given-matrix"){
    if(groups > 1){
      stop("the hypothesis 'given-matrix' can only be tested for one group")
    }
    C <- diag(1, p, p)
    if(is.null(A)){
      Xi <- matrixcalc::vech(diag(1, d, d))
      warning("since no input A for a matrix to be tested is given, the identity matrix is tested")
    }
    else{
      if(!is.matrix(A)){
        stop("the given matrix A must be a matrix with dimensions d x d")
      }
      if(matrixcalc::is.square.matrix(A) & dim(A)[1] == d){
        Xi <- matrixcalc::vech(A)
      }
      else{
        stop("the given matrix A must be a square matrix with dimensions d x d")
      }
    }
    return(TestCovariance_base(X, nv = nv, C = C, Xi = Xi, method = method,
                               repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }
  if(hypothesis == "uncorrelated"){
    if(groups > 1){
      stop("the hypothesis 'uncorrelated' can only be tested for one group")
    }
    C <- diag(1, p, p)[-a,]
    Xi <- rep(0, p - d)
    return(TestCovariance_base(X = X, nv = nv, C = C, Xi = Xi, method = method,
                            repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }
  if(hypothesis == "equal-trace"){
    if(groups == 1){
      stop("the hypothesis 'equal-trace' can only be tested for multiple groups")
    }
    tracevec <- matrix(0, 1, p)
    tracevec[1, a] <- 1
    C <- Pd(groups) %x% tracevec
    Xi <- rep(0, groups)
    return(TestCovariance_base(X = X, nv = nv, C = C, Xi = Xi, method = method,
                               repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }
  if(hypothesis == "equal-diagonals"){
    if(groups == 1){
      stop("the hypothesis 'equal-diagonals' can only be tested for multiple groups")
    }
    C <- Pd(groups) %x% diag(1, p, p)[a,]
    Xi <- rep(0, times = groups * d)
    return(TestCovariance_base(X = X, nv = nv, C = C, Xi = Xi, method = method,
                               repetitions = repetitions, seed = seed, hypothesis = hypothesis))
  }

}

#' @title Test for structure of data's covariance matrix
#'
#' @description This function conducts the test for the covariance matrix of data regarding
#' structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the Anova-Type-Statistics based on a
#' specified number of runs.
#' @param X a matrix containing the observation vectors as columns (one group only)
#' @param structure a character specifying the structure regarding the
#' covariance matrix should be checked. Options are "autoregressive" ("ar"), "FO-autoregressive" ("FO-ar"),
#' "diagonal" ("diag"), "sphericity" ("spher"), "compoundsymmetry" ("cs") and "toeplitz" ("toep").
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is NULL, which means no seed is set.
#' @return an object of the class \code{\link{CovTest}}
#'
#'
#' @references \insertRef{sattler_structures_2024}{CovCorTest}
#'
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' # Select only the males with the diagnosis AD
#' X <- as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
#'                           c("brainrate_temporal", "brainrate_frontal","brainrate_central",
#'                             "complexity_temporal","complexity_frontal", "complexity_central")])
#'
#' TestCovariance_structure(X = X, structure = "diagonal", method = "MC")
#'
#' @export
TestCovariance_structure <- function(X, structure, method = "BT", repetitions = 1000, seed = NULL){

  structure <- tolower(structure)
  method <- toupper(method)
  if(!(method == "MC" | method == "BT")){
    stop("method must be bootstrap ('BT') or Monte-Carlo-technique('MC')")
  }

  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  if(is.list(X)){
    if(length(X) > 1){
      warning("The input X must be a matrix but is a list. Only the first element of the list is used.")
      X <- X[[1]]
    }
    if(length(X) == 1){
      X <- X[[1]]
    }

  }

  n1 <- dim(X)[2]
  d <- dim(X)[1]
  if(d==1){ stop("Structures can be only investigated for more than one dimension") }
  if(!(structure %in% c("autoregressive", "ar", "fo-autoregressive", "fo-ar", "diagonal",
                        "diag", "sphericity", "spher", "compoundsymmetry", "cs", "toeplitz", "toep") )){
    stop("no predefined hypothesis")
  }

  if(d > 1){
    p <- d * (d + 1) / 2
    a <- cumsum(c(1, (d):2))

    vX <- dvech(stats::var(t(X)), a, d, p, inc_diag = TRUE)
    Xq <- apply(X - rowMeans(X), 2, vdtcrossprod, a, d, p)
    HatCov <- stats::var(t(Xq))

    if(structure == "autoregressive" | structure == "ar"){
      C <- matrixcalc::direct.sum(diag(1, d, d), Pd(p - d))
      Xi <- c(rep(1, times = d), rep(0, times = p - d))
      Jacobi <- Jacobian(vX, a, d, p, 'ascending_root_fct')
      HatCovh <- QF(Jacobi, HatCov)

      if(method == "MC") {
        ResamplingResult <- ATSwS(QF(C, HatCovh), repetitions)
      }
      if(method == "BT") {
        ResamplingResult <- sapply(1:repetitions, Bootstrap_trans, n1, a ,d, p,
                                   C, MSroot(HatCov), vX, 'ascending_root_fct')
      }
      Teststatistic <- ATS(n1, ascending_root_fct(vX, a, d), C, HatCovh, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)
    }

    if(structure == "fo-autoregressive" | structure == "fo-ar"){
      C <- Pd(d)
      for(l in 2:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
      C <- matrixcalc::direct.sum(C, Pd(d - 1))
      Xi <- rep(0, times = p + d - 1)
      Jacobi <- Jacobian(vX, a, d, p, 'subdiagonal_mean_ratio_fct')
      HatCovg <- QF(Jacobi, HatCov)

      if(method == "MC"){
        ResamplingResult <- ATSwS(QF(C, HatCovg), repetitions)
      }
      if(method == "BT"){
        ResamplingResult <- sapply(1:repetitions, Bootstrap_trans, n1, a, d, p,
                                   C, MSroot(HatCov), vX, 'subdiagonal_mean_ratio_fct')
      }
      Teststatistic <- ATS(n1, subdiagonal_mean_ratio_fct(vX, a, d), C, HatCovg, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)
    }

    if(structure %in% c("diagonal", "diag", "sphericity", "spher", "compoundsymmetry", "cs", "toeplitz", "toep")){

      Xi <- rep(0, p)

      if(structure == "diagonal" | structure == "diag"){
        C <- matrixcalc::direct.sum(matrix(0, d, d), diag(1, p - d, p - d))
      }
      if(structure == "sphericity" | structure == "spher"){
        C <- matrixcalc::direct.sum(Pd(d), diag(1, p - d, p - d))
      }
      if(structure == "compoundsymmetry" | structure == "cs"){
        C <- matrixcalc::direct.sum(Pd(d), Pd(p - d))
      }
      if(structure == "toeplitz" | structure == "toep"){
        C <- Pd(d)
        for(l in 2:d){
          C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
        }
      }

      if(method == "MC"){
        ResamplingResult <- ATSwS(QF(C, HatCov), repetitions)
      }
      if(method == "BT"){
        ResamplingResult <- sapply(1:repetitions, Bootstrap, n1, C, MSroot(HatCov))
      }

      Teststatistic <- ATS(n1, vX, C, HatCov, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)

    }
  }
  CovTest <- list("method" = "Covariance",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = HatCov,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = structure,
                  "nv" = 1)

  class(CovTest) <- "CovTest"

  return(CovTest)
}


