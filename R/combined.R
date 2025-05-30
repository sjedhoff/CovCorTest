#' @title The Taylor-based Monte-Carlo-approximation for a combined test
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for the combined test. After receiving
#'  some auxiliary matrices and data, the Monte-Carlo observations are
#'  generated, and different parts of the final sum are defined. Based on this,
#'  a number of the Taylor-based vectors are calculated, where the number can be
#'  chosen.
#' @param repetitions a number specifying the number of runs for the
#' approximation
#' @param MSrootHatCov the matrix root of the covariance matrix for the Taylor
#' observations
#' @param CorData the calculated correlation matrix
#' @param MvrH1 an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param MvrH2 an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param M4 an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorized
#'  covariances to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorized
#' covariances to vectorized correlations
#' @param nv vector of sample sizes
#' @return a matrix containing the values of the test vector for a number of
#' repetitions
#'
#' @export
TaylorCombined <- function(repetitions, MSrootHatCov, CorData, MvrH1, MvrH2,
                           M4, L, P, Q, nv) {
  XPB <- mapply(generateData, MSrootHatCov, repetitions, SIMPLIFY = FALSE)
  XTaylorUpper <- lapply(XPB, function(X)
    P %*% X)
  XPB <- mapply("%*%", MvrH2, XPB, SIMPLIFY = FALSE)

  vechCorData <- lapply(CorData, matrixcalc::vech)
  DvechCorDataM <- lapply(lapply(vechCorData, as.vector), "*", M4)

  PUX <- lapply(XPB, function(X) P %*% X)
  MUX0 <- lapply(XPB, function(X) (M4 %*% Q) %*% X)
  HX <- mapply(Qvech, PUX, repetitions, SIMPLIFY = FALSE)

  Part1 <- mapply(function(A, B)
    A %*% B, MvrH1, XPB, SIMPLIFY = FALSE)
  Part2 <- mapply("*", lapply(vechCorData, as.vector), HX, SIMPLIFY = FALSE)
  Part3 <- mapply("*", XPB, MUX0, SIMPLIFY = FALSE)
  Part4 <- mapply("%*%", DvechCorDataM, HX, SIMPLIFY = FALSE)

  XTaylorLower <- mapply(function(Part1, L, Part2, Part3, Part4, nv) {
    return((Part1 + L %*% (1 / 4 * Part2 - 1 / 2 * Part3 + 3 / 8 * Part4) /
              sqrt(nv)))
  }, Part1, list(L), Part2, Part3, Part4, nv, SIMPLIFY = FALSE)

  TTaylor <- sqrt(sum(nv)) * (
    1 / sqrt(nv[[1]]) * rbind(XTaylorUpper[[1]], XTaylorLower[[1]]) - 1 /
      sqrt(nv[[2]]) *rbind(XTaylorUpper[[2]], XTaylorLower[[2]])
  )
  return(TTaylor)
}


#' @title Combined test for equality of covariance matrices and correlation
#' matrices
#'
#' @description For two groups a combined test for equality of covariance
#' matrices and equality of correlation matrices between these groups is
#' conducted. Both hypotheses can be rejected or only the larger one, the
#' equality of the covariance matrices
#' @param X  a list or matrix containing the observation vectors.
#' In case of a list,each matrix in this list is another group, where the
#' observation vectors are the columns. For a matrix, all groups are together in
#' one matrix and nv is used to indicate the group sizes.
#' @param nv vector of group sizes
#' @param repetitions a scalar,  indicate the number of runs for the chosen
#' method. The predefined value is 1.000, and the number should not be
#' below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined
#' values is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return an object of the class \code{\link{CovTest}}.
#'
#' @import MANOVA.RM
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' vars <- colnames(EEGwide)[1:6]
#'
#' # Part the data into six groups of sex and diagnosis
#' X_list <- list(t(EEGwide[EEGwide$sex=="M" & EEGwide$diagnosis=="AD",vars]),
#'                t(EEGwide[EEGwide$sex=="M" & EEGwide$diagnosis=="MCI",vars]))
#'
#' nv <- unlist(lapply(X_list, ncol))
#'
#' test_combined(X_list,nv,seed=31415)
#'
#' @export
test_combined <- function(X, nv = NULL, repetitions = 1000, seed = NULL) {
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

  listcheck <- Listcheck(X, nv)

  X <- listcheck[[1]]
  nv <- listcheck[[2]]
  groups <- length(nv)

  dimensions <- sapply(X, dim)[1, ]
  if (max(dimensions) != mean(dimensions)) {
    stop("dimensions do not accord")
  }
  else{
    d <- dimensions[1]
    if(d == 1) {
      stop("Correlation is only defined for dimension higher than one")
    }
    if(groups != 2) {
      stop("This test is only defined for exactly two groups")
    }
    if(d > 1) {
      p <- d * (d + 1) / 2
      a <- cumsum(c(1, (d):2))
      N <- sum(nv)
      H <- matrix(rep(a, d), d, d, byrow = TRUE)
      L <- diag(1, p, p)[-a, ]
      M <- matrix(0, p, p)
      for (l in 1:p) {
        M[l, c(vechp(H)[l], vechp(t(H))[l])] <- 1
      }
      M1 <- L %*% M

      Datac <- lapply(X, function(x)
        x - rowMeans(x))
      VarData <- lapply(X, function(X)
        stats::var(t(X)))
      CorData <- lapply(VarData, stats::cov2cor)
      vCorData <- lapply(CorData, vechp)
      Teststatistic <- sqrt(N) * (c(diag(VarData[[1]]), vCorData[[1]]) -
                                    c(diag(VarData[[2]]), vCorData[[2]]))
      P <- diag(1, p, p)[a, ]
      Q <- diag(as.vector(matrixcalc::vech(diag(1, d, d))), p, p)

      DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
      HatCov <- lapply(DataQ, function(X)
        stats::var(t(X)))
      MvrH1 <- list()
      MvrH2 <- list()
      for (i in seq_along(nv)) {
        MvrH1[[i]] <- (L - 1 / 2 * vCorData[[i]] * M1)
        MvrH2[[i]] <- sqrt(diag(as.vector(1 / vtcrossprod(
          matrix(matrixcalc::vech(VarData[[i]])[a])
        )), p, p))
      }

      MSrootHatCov <- lapply(HatCov, MSroot)
      ResamplingResult <- TaylorCombined(repetitions, MSrootHatCov, CorData,
                                         MvrH1, MvrH2, M, L, P, Q, nv)
      beta <- 1 - c(max(rowMeans((ResamplingResult < Teststatistic)[1:d, ])),
                    max(rowMeans((ResamplingResult < Teststatistic)[-(1:d), ])))
      Tquantile <- apply(ResamplingResult, 1, stats::quantile, 1 - beta,
                         Type = 2)
      alpha <- c(mean(apply((ResamplingResult > Tquantile[1, ]), 2, max)),
                 mean(apply((ResamplingResult > Tquantile[2, ]), 2, max)))

      if (!is.null(seed))
        (set.seed(NULL))

      CombTest <- list("pvalue_Variances" = alpha[1],
                       "pvalue_Correlations" = alpha[2],
                       "pvalue_Total" = min(alpha),
                       "Teststatistic" = max(Teststatistic),
                       "repetitions" = repetitions,
                       "nv" = nv)

      class(CombTest) <- "CombTest"
      return(CombTest)
    }
  }
}
