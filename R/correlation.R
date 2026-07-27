#' Test for Correlation Matrices
#'
#' @description This function conducts statistical tests for hypotheses
#' regarding correlation matrices.
#' Users can either select from predefined hypotheses or
#' provide their own contrast matrix `C` and vector `Xi` for custom hypotheses.
#' It supports both bootstrap and Monte Carlo resampling methods to
#' obtain the p-value of the ANOVA-type statistic (ATS).
#'
#' @param X A list or a matrix containing the observation vectors. If a list,
#'  each entry is a group, with observations as rows. If a matrix,
#'  all groups are combined, and `nv` must be used to indicate group sizes.
#' @param nv (Optional) A vector indicating group sizes, needed when
#'  `X` is a combined matrix or for multiple groups.
#' @param hypothesis A character specifying one of the predefined hypotheses:
#' \itemize{
#'     \item `"equal-correlated"` — equal correlation matrices
#'     \item `"uncorrelated"` — test if variables are uncorrelated
#'     (single group only)
#'   }
#'  If `C` and `Xi` are provided, this can be set to `NULL`.
#' @param C A hypothesis matrix specifying the null hypothesis.
#'  Optional if \code{hypothesis} is provided.
#' @param Xi A numeric vector specifying the expected values of the contrasts
#' under the null hypothesis. Optional if \code{hypothesis} is provided.
#' @param AM Binary variable indicating whether the computations should use an
#' alternative hypothesis matrix as proposed by
#' \insertCite{Sattler2025;textual}{CovCorTest}. This matrix has fewer rows than
#'  the original hypothesis matrix while leaving the resulting test statistics
#'  unchanged. By default, this alternative matrix is used.
#' @param method Character string indicating the resampling method to use.
#'  One of \code{"BT"} (bootstrap), \code{"MC"} (Monte Carlo), or
#'  \code{"TAY"} (Taylor approximation).
#' @param repetitions Integer. Number of resampling repetitions
#' (default is 1000).
#'
#' @return An object of class \code{"CovTest"}.
#'
#' @references
#' Sattler, P. and Pauly, M. (2024). Testing hypotheses about correlation matrices in general MANOVA designs. \emph{TEST}, 33(2), 496--516. \doi{10.1007/s11749-023-00906-6}
#'
#' @import MANOVA.RM
#' @examples
#' # Example with one group:
#' set.seed(31415)
#' X <- matrix(rnorm(5 * 100), ncol = 5)
#' test_correlation(X, hypothesis = "uncorrelated",
#'                   method = "BT", repetitions = 100)
#'
#' @export
test_correlation <- function(X, nv = NULL,
                             C = NULL, Xi = NULL,
                             AM = 1,
                             hypothesis = NULL,
                             method = "BT",
                             repetitions = 1000) {
  method <- toupper(method)
  if(!(method %in% c("MC", "BT", "TAY"))){
    stop("method must be bootstrap ('BT'), Monte-Carlo-technique('MC') or
         Taylor-based Monte-Carlo-approach('TAY')")
  }
  if (!(AM %in% c(0, 1))) {
    stop("AM must be either 0 or 1.")
  }

  listcheck <- Listcheck(X, nv)
  X <- listcheck[[1]]
  nv <- listcheck[[2]]

  if(!is.null(hypothesis)){
    # infer d and p
    if (is.null(nv)) {
      d <- dim(X)[1]
    } else {
      d <- nrow(X[[1]])
      groups <- length(nv)
    }
    p <- d * (d + 1) / 2
    hypothesis <- tolower(hypothesis)

    if (!(hypothesis %in% c("equal-correlated", "uncorrelated"))) {
      stop("hypothesis must be one of 'equal-correlated' or 'uncorrelated'")
    }

    if (hypothesis == "equal-correlated") {
      if (is.null(nv)) {
        C <- Pd(p - d)
        Xi <- rep(0, p - d)
      } else {
        C <- Pd(groups) %x% diag(1, p - d)
        Xi <- rep(0, groups * (p - d))
      }
    }

    if (hypothesis == "uncorrelated") {
      if (!is.null(nv)) {
        stop("the hypothesis 'uncorrelated' can only be tested
             for a single group")
      }
      C <- diag(1, p - d)
      Xi <- rep(0, p - d)
    }
  } else {
    # if no hypothesis, use C and Xi
    if (is.null(C) || is.null(Xi)) {
      stop("Either provide 'hypothesis' or both 'C' and 'Xi'")
    }
  }




  # one group
  if(is.null(nv)){
    nv <- dim(X)[2]
    d <- dim(X)[1]
    groups <- 1
  }
  # multiple groups
  else{
    d <- unlist(lapply(X, nrow))[1]
    groups <- length(nv)
  }

  # Check dimensions of C and Xi
  pu <- d*(d-1)/2
  if(!is.matrix(C)){
    stop("C must be a matrix")
  }
  if( (nrow(C) != length(Xi)) || (ncol(C) != groups*pu) ){
    stop("dimensions of C and Xi do not align")
  }

  p <- d*(d+1)/2
  a <- cumsum(c(1,(d):2))
  N <- sum(nv)
  H <- matrix(rep(a,d), d, d, byrow=TRUE)
  L <- diag(1,p,p)[-a,]
  M4 <- matrix(0,p,p)

  for(i in 1:p){
    M4[i, c(matrixcalc::vech(t(H))[i], matrixcalc::vech(H)[i])] <- 1
  }
  M1 <- L%*%M4

  # one group
  if(length(nv) == 1){
    VarData <- stats::var(t(X))
    CorData <- stats::cov2cor(VarData)
    vCorData <- vechp(CorData)
    Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod), nrow=p, ncol=nv)
    HatCov <- stats::var(t(Xq))
    MvrH1 <- (L-1/2*vCorData*M1)
    MvrH2 <- sqrt(diag(as.vector(1/vtcrossprod(matrix(
      matrixcalc::vech(VarData)[a]))),p,p))
    Upsi <- QF(MvrH1%*%MvrH2,HatCov)
    MSrootUpsi <- MSroot(Upsi)
    StUpsi <- QF(MvrH2,HatCov)
    MSrootStUpsi <- MSroot(StUpsi)
  }
  # multiple groups
  else{
    Datac <- lapply(X, function(x) x-rowMeans(x))
    VarData <- lapply(X, function(X) stats::var(t(X)))
    CorData <- lapply(VarData, stats::cov2cor)
    vCorData <- lapply(CorData, vechp)
    DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
    HatCov <- lapply(DataQ, function(X) stats::var(t(X)))
    MvrH <- StUpsi <- list()

    for (i in seq_along(nv)){
      MvrH[[i]] <- (L - 1 / 2 * vCorData[[i]] * M1)
      StUpsi[[i]] <- QF(sqrt(diag(as.vector(1 / vtcrossprod(matrix(
        matrixcalc::vech(VarData[[i]])[a]))), p, p)), HatCov[[i]])
    }
    Upsi_list <- mapply(QF, MvrH, StUpsi, SIMPLIFY = FALSE)
    MSrootStUpsi <- lapply(StUpsi, MSroot)
    MSrootUpsi <- lapply(Upsi_list, MSroot)
    kappainvv <- N / nv
    Upsi <- WDirect.sumL(Upsi_list, kappainvv)
  }
  C_original <- C
  Xi_original <- Xi

  if (AM == 1) {
    AlternativeHypothesis <- CompanionHypothesis(C, Xi)
    C <- AlternativeHypothesis$L
    Xi <- AlternativeHypothesis$ytilde
  }


  if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsi), repetitions) }
  if(method == "BT"){ ResamplingResult <- vapply(X = 1:repetitions,
                                                 FUN = Bootstrap,
                                                 FUN.VALUE = numeric(1),
                                                 nv, C, MSrootUpsi)
  }
  if(method == "TAY"){
    P <- diag(1,p,p)[a,]
    Q <- diag(as.vector(matrixcalc::vech(diag(1,d,d))),p,p)
    Trace <- sum(diag(QF(C,Upsi)))
    if(length(nv) == 1){ # one group
      ResamplingResult <- Tayapp1G(repetitions, C, MSrootStUpsi, CorData,
                                   MvrH1, Trace, M4, L, P, Q, nv, NULL, NULL)
    }
    else{ #multiple groups
      ResamplingResult <- TayappMG(repetitions, C, MSrootStUpsi, CorData,
                                   MvrH, Trace, M4, L, P, Q, nv)
    }
  }


  Teststatistic <- ATS(N, unlist(vCorData), C, Upsi, Xi)
  pvalue <- mean(ResamplingResult > Teststatistic)

  CovTest <- list("method" = "Correlation",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsi,
                  "C" = C,
                  "Xi" = Xi,
                  "AM" = AM,
                  "C_original" = C_original,
                  "Xi_original" = Xi_original,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = hypothesis,
                  "nv" = nv)

  class(CovTest) <- "CovTest"

  return(CovTest)
}


#' @title Test for structure of data's correlation matrix
#'
#' @description With this function the correlation matrix of data can be checked
#' for one of the predefined structures. Depending on the chosen method a
#' bootstrap, the Taylor-based Monte-Carlo approach or  Monte-Carlo-technique
#' is used to calculate the p-value of the Anova-type-statistic(ATS) based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as rows (one group)
#' @param structure a character specifying the structure regarding them the
#' correlation matrix should be checked. Options are "Hautoregressive" ("Har"),
#' "diagonal" ("diag"), "Hcompoundsymmetry" ("Hcs") and "Htoeplitz" ("Htoep"),
#' "banded" ("band"), and "banded-toeplitz" ("band-toep"). For banded structures,
#' all correlations outside the specified bandwidth, i.e. from the
#' bandwidth + 1-th off-diagonal onward, are assumed to be zero.
#' @param AM Binary variable indicating whether the computations should use an
#' alternative hypothesis matrix as proposed by
#' \insertCite{Sattler2025;textual}{CovCorTest}. This matrix has fewer rows than
#' the original hypothesis matrix while leaving the resulting test statistics
#' unchanged. By default, this alternative matrix is used.
#' @param method a character, to chose whether bootstrap("BT") or Taylor-based
#' Monte-Carlo-approach("TAY") or Monte-Carlo-technique("MC") is used, while
#' bootstrap is the predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen
#' method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param bandwidth Optional positive integer specifying the number of
#' off-diagonals allowed to be non-zero for banded correlation structures.
#' Only used for \code{"banded"} and \code{"banded-toeplitz"} structures,
#' with default \code{NA}.
#' @return an object of the class \code{\link{CovTest}}
#'
#'
#' @references \insertRef{sattler_structures_2024}{CovCorTest}
#'
#' @import MANOVA.RM
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' # Select only the males with the diagnosis AD
#' X <- as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
#'              c("brainrate_temporal", "brainrate_frontal","brainrate_central",
#'              "complexity_temporal","complexity_frontal",
#'              "complexity_central")])
#'
#' test_correlation_structure(X = X, structure = "diagonal", method = "MC")
#'
#' @export
test_correlation_structure <- function(X, structure, AM = 1, method = "BT",
                                      repetitions = 1000, bandwidth = NA){

  structure <- tolower(structure)
  method <- toupper(method)
  if(!(method == "MC" || method == "BT" || method == "TAY")){
    stop("method must be bootstrap ('BT'), Monte-Carlo-technique('MC') or
         Taylor-based Monte-Carlos-approach('Tay')")
  }
  if (!(AM %in% c(0, 1))) {
    stop("AM must be either 0 or 1.")
  }

  if(is.list(X) & !is.data.frame(X)){
    if(length(X) > 1){
      warning("The input X must be a matrix or data frame but is a list. Only the first element of the list is used.")
    }
    listcheck <- Listcheck(X[[1]])
  }
  else{
    listcheck <- Listcheck(X)
  }

  X <- listcheck[[1]]
  n1 <- dim(X)[2]
  d <- dim(X)[1]

  if(!(structure %in% c("hautoregressive", "har",
                        "diagonal", "diag",
                        "hcompoundsymmetry", "hcs",
                        "htoeplitz", "htoep",
                        "banded", "band",
                        "banded-toeplitz", "band-toep"))){
    stop("no predefined hypothesis")
  }
  if (structure %in% c("banded", "band", "banded-toeplitz", "band-toep")) {
    bandwidth <- CheckBandwidth(bandwidth, d)
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
  M1 <- L %*% M
  VarData <- stats::var(t(X))
  CorData <- stats::cov2cor(VarData)
  vCorData <- dvech(CorData, a, d, p, inc_diag = FALSE)
  Xq <- matrix(apply(X - rowMeans(X), 2, vtcrossprod), nrow = p, ncol = n1)
  HatCov <- stats::var(t(Xq))
  MvrH1 <- (L - 1 / 2 * vechp(CorData) * M1)

  MvrH2 <- sqrt(diag(as.vector(1 / vtcrossprod(matrix(
                                        matrixcalc::vech(VarData)[a]))), p, p))
  Atilde <- matrix(0, pu, pu)
  for(l in 1:(d - 1)){
    for(k in 1:(d - l)){
      Atilde[a[l] + k - l, a[k] - (k - l)] <- 1
    }
  }

  Upsidv <- QF(Atilde %*% MvrH1 %*% MvrH2, HatCov)
  Xi <- rep(0, pu)

  if(structure == "hautoregressive" || structure == "har"){
    Jacobi <- Jacobian(vCorData, a, d, p, fun = "subdiagonal_mean_ratio_cor")
    Upsidvhtilde <- QF(Jacobi, Upsidv)
    C <- Pd(d - 1)
    for(l in 3:d){
      C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
    }
    C <- matrixcalc::direct.sum(C, Pd(d - 2))
    Xi <- rep(0, p - 2)

    C_original <- C
    Xi_original <- Xi

    if (AM == 1) {
      AlternativeHypothesis <- CompanionHypothesis(C, Xi)
      C <- AlternativeHypothesis$L
      Xi <- AlternativeHypothesis$ytilde
    }

    Teststatistic <- ATS(n1, subdiagonal_mean_ratio_cor(vCorData, a, d), C,
                         Upsidvhtilde, Xi)
    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidvhtilde),
                                                  repetitions) }
    if(method == "BT"){ ResamplingResult <- vapply(X = 1:repetitions,
                                            FUN = Bootstrap_trans,
                                            FUN.VALUE = numeric(1),
                                            n1, a, d,
                                            p, C, MSroot(Upsidv),
                                            vCorData,
                                            fun = "subdiagonal_mean_ratio_cor")
    }
    if(method == "TAY"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidvhtilde)))

      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                  MvrH1, Trace, M, L, P, Q, n1, Atilde, Jacobi)
    }
  }
  else{
    if(structure == "diagonal" || structure == "diag"){
      C <- diag(1, pu, pu)
    }
    if(structure == "hcompoundsymmetry" || structure == "hcs"){
      C <- Pd(pu)
    }
    if(structure == "htoeplitz" | structure == "htoep"){
      C <- Pd(d - 1)
      for(l in 3:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
    }

    if(structure == "banded" | structure == "band"){

      allowed <- bandwidth * d - bandwidth * (bandwidth + 1) / 2
      outside_positions <- (allowed + 1):pu

      C <- diag(1, pu, pu)[outside_positions, , drop = FALSE]
      Xi <- rep(0, nrow(C))
    }

    if(structure == "banded-toeplitz" | structure == "band-toep"){

      allowed <- bandwidth * d - bandwidth * (bandwidth + 1) / 2
      n_outside <- pu - allowed

      C <- Pd(d - 1)

      if(bandwidth >= 2){
        for(l in 2:bandwidth){
          C <- matrixcalc::direct.sum(C, Pd(d - l))
        }
      }

      C <- matrixcalc::direct.sum(C, diag(1, n_outside, n_outside))

      Xi <- rep(0, nrow(C))
    }
    C_original <- C
    Xi_original <- Xi

    if (AM == 1) {
      AlternativeHypothesis <- CompanionHypothesis(C, Xi)
      C <- AlternativeHypothesis$L
      Xi <- AlternativeHypothesis$ytilde
    }
    Teststatistic <- ATS(n1, vCorData, C, Upsidv, Xi)

    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidv), repetitions) }
    if(method == "BT"){ ResamplingResult <- vapply(X = 1:repetitions,
                                                   FUN = Bootstrap,
                                                   FUN.VALUE = numeric(1),
                                                   n1,
                                                   C, MSroot(Upsidv)) }
    if(method == "TAY"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidv)))
      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                   MvrH1, Trace, M, L, P, Q,n1)
    }
  }

  pvalue <- mean(ResamplingResult > Teststatistic)


  CovTest <- list("method" = "Correlation",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsidv,
                  "C" = C,
                  "Xi" = Xi,
                  "AM" = AM,
                  "C_original" = C_original,
                  "Xi_original" = Xi_original,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = structure,
                  "nv" = n1)

  class(CovTest) <- "CovTest"

  return(CovTest)

}
