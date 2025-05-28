#' @title Test for Covariance Matrices
#'
#' @description This function conducts statistical tests for hypotheses regarding covariance matrices.
#' Users can either select from predefined hypotheses (e.g., equal covariance, equal trace, etc.) or
#' provide their own contrast matrix `C` and vector `Xi` for custom hypotheses. It supports both
#' bootstrap and Monte Carlo resampling methods to obtain the p-value of the ANOVA-type statistic (ATS).
#'
#' @param X A list or a matrix containing the observation vectors. If a list, each entry is a group,
#'   with observations as columns. If a matrix, all groups are combined, and `nv` must be used to indicate group sizes.
#' @param nv (Optional) A vector indicating group sizes, needed when `X` is a combined matrix or for multiple groups.
#' @param hypothesis A character specifying one of the predefined hypotheses:
#'   \itemize{
#'     \item `"equal"` — equal covariance matrices
#'     \item `"equal-trace"` — equal traces across groups
#'     \item `"equal-diagonals"` — equal variances across groups
#'     \item `"given-trace"` — test against a given trace (single group only)
#'     \item `"given-matrix"` — test against a given covariance matrix (single group only)
#'     \item `"uncorrelated"` — test if variables are uncorrelated (single group only)
#'   }
#'   If `C` and `Xi` are provided, this can be set to `NULL`.
#' @param A Optional scalar or matrix to define the hypothesis value when `hypothesis` is `"given-trace"` (scalar)
#'   or `"given-matrix"` (matrix). Ignored for other hypotheses.
#' @param C (Optional) A user-defined contrast matrix for testing custom hypotheses. Must match dimensions with `Xi`.
#' @param Xi (Optional) A numeric vector used in combination with `C` to specify a custom hypothesis.
#' @param method A character indicating the resampling method: `"BT"` (Bootstrap) or `"MC"` (Monte Carlo).
#' @param repetitions Number of repetitions to use for the resampling method (default: 1000, should be >= 500).
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class \code{\link{CovTest}}.
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
#' C <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1),
#'             nrow = 1, ncol = 21)
#' Xi <- 2
#'
#' test_covariance(X = X, nv = NULL, C = C, Xi = Xi, method = "BT",
#'             repetitions = 1000, seed = 31415)
test_covariance <- function(X, nv = NULL, C = NULL, Xi = NULL,
                           hypothesis = NULL, A = NULL,
                           method = "MC", repetitions = 1000, seed = NULL) {

  method <- toupper(method)
  if(!(method %in% c("MC", "BT"))){
    stop("method must be 'MC' or 'BT'")
  }

  listcheck <- Listcheck(X, nv)
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

  # check for hypothesis
  if(!is.null(hypothesis)){
    hypothesis <- tolower(hypothesis)
    if(!(hypothesis %in%
         c("equal", "equal-trace", "equal-diagonals", "given-trace",
           "given-matrix", "uncorrelated"))){
      stop("no predefined hypothesis")
    }

  if(!is.null(A) & !(hypothesis %in% c("given-trace", "given-matrix"))){
    warning(paste0("the input argument A is not used, since the selected
                   hypothesis is '", hypothesis, "'"))
  }

    switch(hypothesis,
           "equal" = {
             if(groups == 1){
               C <- Pd(d) %*% diag(1, p, p)[a,]
               Xi <- rep(0, d)
             }
             else{
               C <- Pd(groups) %x% diag(1, p, p)
               Xi <- rep(0, p*groups)
             }
           },
           "equal-trace" = {
             if(groups == 1){
               stop("the hypothesis 'equal-trace' can only be tested for
           multiple groups")
             }
             tracevec <- matrix(0, 1, p)
             tracevec[1, a] <- 1
             C <- Pd(groups) %x% tracevec
             Xi <- rep(0, groups)
           },
           "equal-diagonals" = {
             if(groups == 1){
               stop("the hypothesis 'equal-diagonals' can only be tested for
           multiple groups")
             }
             C <- Pd(groups) %x% diag(1, p, p)[a,]
             Xi <- rep(0, times = groups * d)
           },
           "given-trace" = {
             if(groups > 1){
               stop("the hypothesis 'given-trace' can only be tested for one group")
             }
             tracevec <- matrix(0, 1, p)
             tracevec[1,a] <- 1
             C <- tracevec
             if(is.null(A)){
               A <- 1
               warning("since no input A for a trace to be tested is given, a trace of
              1 is tested")
             }
             if(!is.numeric(A) | length(A) != 1){
               stop("for testing the trace the input A must be a scalar")
             }
             Xi <- A
           },
           "given-matrix" = {
             if(groups > 1){
               stop("the hypothesis 'given-matrix' can only be tested for one group")
             }
             C <- diag(1, p, p)
             if(is.null(A)){
               Xi <- matrixcalc::vech(diag(1, d, d))
               warning("since no input A for a matrix to be tested is given, the
              identity matrix is tested")
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
           },
           "uncorrelated" = {
             if(groups > 1){
               stop("the hypothesis 'uncorrelated' can only be tested for one group")
             }
             C <- diag(1, p, p)[-a, , drop = FALSE]
             Xi <- rep(0, p - d)
           }
    )
  }

  # seed
  if(!is.null(seed)){
    if(exists(".Random.seed")){
      old_seed <- .Random.seed
      on.exit({ .Random.seed <<- old_seed })
    }
    if(!exists(".Random.seed")){
      on.exit({ set.seed(NULL) })
    }
    set.seed(seed)
  }

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
    ResamplingResult <- vapply(1:repetitions, Bootstrap, numeric(1), nv, C,
                               MSrootHatCov)
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


#' @title Test for structure of data's covariance matrix
#'
#' @description This function conducts the test for the covariance matrix of
#' data regarding structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the
#' Anova-type-statistic(ATS) based on a specified number of runs.
#' @param X a matrix containing the observation vectors as columns
#' (one group only)
#' @param structure a character specifying the structure regarding the
#' covariance matrix should be checked. Options are "autoregressive" ("ar"),
#' "FO-autoregressive" ("FO-ar"), "diagonal" ("diag"), "sphericity" ("spher"),
#' "compoundsymmetry" ("cs") and "toeplitz" ("toep").
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the
#' predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen
#' method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed a seed, if it should be set for reproducibility. Predefined value
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
#'                           c("brainrate_temporal", "brainrate_frontal",
#'                           "brainrate_central","complexity_temporal",
#'                           "complexity_frontal", "complexity_central")])
#'
#' test_covariance_structure(X = X, structure = "diagonal", method = "MC")
#'
#' @export
test_covariance_structure <- function(X, structure, method = "BT",
                                     repetitions = 1000, seed = NULL){

  structure <- tolower(structure)
  method <- toupper(method)
  if(!(method == "MC" | method == "BT")){
    stop("method must be bootstrap ('BT') or Monte-Carlo-technique('MC')")
  }

  if(!is.null(seed)){
    if(exists(".Random.seed")) {
      old_seed <- .Random.seed
      on.exit({ .Random.seed <<- old_seed })
    }
    else{
      on.exit({ set.seed(NULL) })
    }
    set.seed(seed)
  }

  if(is.list(X)){
    if(length(X) > 1){
      warning("The input X must be a matrix but is a list. Only the first
              element of the list is used.")
      X <- X[[1]]
    }
    if(length(X) == 1){
      X <- X[[1]]
    }

  }

  n1 <- dim(X)[2]
  d <- dim(X)[1]
  if(d==1){ stop("Structures can be only investigated for more than one
                 dimension") }
  if(!(structure %in% c("autoregressive", "ar", "fo-autoregressive", "fo-ar",
                        "diagonal", "diag", "sphericity", "spher",
                        "compoundsymmetry", "cs", "toeplitz", "toep") )){
    stop("no predefined hypothesis")
  }

  if(d > 1){
    p <- d * (d + 1) / 2
    a <- cumsum(c(1, (d):2))

    vX <- dvech(stats::var(t(X)), a, d, p, inc_diag = TRUE)
    Xq <- apply(X - rowMeans(X), 2, vdtcrossprod, a, d, p)
    HatCov <- stats::var(t(Xq))

    if(structure == "autoregressive" | structure == "ar"){
      C <- diag(1,d,d)
      for(l in 2:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
      C <- matrixcalc::direct.sum(C, Pd(d - 1))
      Xi <- c(rep(1,d),rep(0, times = p - 1))
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
        ResamplingResult <- vapply(1:repetitions, Bootstrap_trans,
                                   FUN.VALUE = numeric(1),
                                   n1, a, d, p,
                                   C, MSroot(HatCov), vX,
                                   'subdiagonal_mean_ratio_fct')
      }
      Teststatistic <- ATS(n1, subdiagonal_mean_ratio_fct(vX, a, d), C,
                           HatCovg, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)
    }

    if(structure %in% c("diagonal", "diag", "sphericity", "spher",
                        "compoundsymmetry", "cs", "toeplitz", "toep")){

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
        ResamplingResult <- vapply(1:repetitions, Bootstrap,
                                   FUN.VALUE = numeric(1),
                                   n1,
                                   C, MSroot(HatCov))
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


