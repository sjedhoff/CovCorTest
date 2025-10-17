#' Extend a matrix to full rank using identity columns
#'
#' This function extends a matrix `V` (assumed to have fewer than `p` columns)
#' to full rank by iteratively appending standard basis vectors from the identity matrix,
#' ensuring that each added column increases the rank.
#'
#' @param V A numeric matrix with `p` rows and fewer than `p` columns.
#'
#' @return A numeric matrix of dimension `p x p` with full column rank.
#'
#' @details
#' The function assumes that `nrow(V)` equals the intended dimension `p` of the parameter space.
#' It adds columns from the identity matrix (i.e., standard basis vectors) only when they increase the rank,
#' and stops when full rank is achieved.
#'
#'
get_extended_matrix <- function(V) {
  p <- nrow(V)  # Number of rows of V (assumed to be the dimension of the parameter space)

  # If V is square, it's already p x p, and extending it makes no sense in the model
  if(p == ncol(V)){
    stop("V is already a square matrix.")  # Error: V should not already be square
  }

  unit_vec <- diag(p)  # Identity matrix of size p (used to extract standard basis vectors)
  i <- 1  # Index for standard basis vectors

  rk <- qr(V)$rank  # Compute the current rank of V

  # Keep adding standard basis vectors as columns until full rank is reached
  while(rk < p){
    V_new <- cbind(V, unit_vec[, i])  # Add the i-th standard basis vector to V

    # Only accept the new column if it increases the rank
    if (qr(V_new)$rank > rk) {
      V <- V_new           # Update V
      rk <- qr(V_new)$rank  # Update rank
    }

    i <- i + 1  # Move to the next basis vector
  }

  return(V)  # Return the extended matrix
}



#' Construct hypothesis matrix and vector from linear covariance model structure
#'
#' Computes a hypothesis matrix `C` and hypothesis vector `zeta` based on a given
#' parameter vector `v0` and a matrix `V` representing the model structure
#' (e.g., vectorised components of a linear covariance structure model).
#'
#' @param v0 A numeric vector of length `p` (number of parameters). Represents the
#'   parameter vector at which the hypothesis is to be evaluated.
#' @param V A numeric matrix of size `p x q`, representing the structured design or constraint matrix
#'   for the model. Here, the vectorised matrices from the linear covariance structure model build the columns of
#'   V.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`Hypothesenmatrix`}{A numeric matrix `C` such that the hypothesis can be written as `C %*% theta = zeta`}
#'     \item{`Hypothesenvector`}{The numeric vector `zeta`, computed as `C %*% v0`}
#'   }
#'
#' @details
#' The function extends `V` to full rank using \link{get_extended_matrix}, constructs a contrast
#' matrix `E` for the complement of the model-implied space, and computes the corresponding hypothesis matrix `C`.
#'
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' X <- t(as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
#'                          c("brainrate_temporal", "brainrate_frontal","brainrate_central",
#'                             "complexity_temporal","complexity_frontal", "complexity_central")]))
#' v0 <- rep(0,21)
#' v_auxiliary <- c(1, rep(0,5), 1, rep(0,4), 1, rep(0,3), 1, rep(0,2), 1, 0, 1)
#' V <- cbind(v_auxiliary, 1-v_auxiliary)
#' h <- get_hypothesis(v0,V)
#' set.seed(123)
#' test_covariance(X = X,C = h$hypothesis_matrix, Xi = h$hypothesis_vector,
#'                 method = "MC", repetitions = 1000)
#'
#' @export
# Function: get_hypothesis
# Purpose: Construct a hypothesis matrix C and a hypothesis vector zeta from a given vector v0 and matrix V,
# containing the vectorised matrices from the linear covariance structure model
#'
#' @references
#' Sattler, P. and Dobler, D. (2025). Testing for patterns and structures in covariance and correlation matrices. \emph{arXiv preprint} \url{https://arxiv.org/abs/2310.11799}
#'
get_hypothesis <- function(v0, V) {
  if(!is.matrix(V)){
    stop("V must be a matrix.")
  }
  if(!is.vector(v0)){
    stop("v0 must be a vector.")
  }
  if(!is.numeric(V) | !is.numeric(v0)){
    stop("V and v0 have to be numeric.")
  }
  p <- nrow(V)  # Number of parameters (rows in V)
  q <- ncol(V)  # Number of constraints or model dimensions (columns in V)
  if(q >= p){
    stop("The matrix V must have more columns than rows.")
  }
  # Error handling: v0 must be of the same dimension as the number of parameters
  if(length(v0) != p){
    stop("v0 must have the same length as rows in V.")
  }

  V_ext <- get_extended_matrix(V)  # Extend V to full rank (if not already)

  # Select the rows (q+1 to p) of the identity matrix to build the contrast matrix E
  # This corresponds to the "unexplained" part of the parameter space
  E <- diag(1, p, p)[(q + 1):p, ]

  # Compute the hypothesis matrix C by transforming E with the inverse of the extended V
  C <-  E %*% solve(V_ext)

  zeta <- C %*% v0  # Compute the transformed (reduced) hypothesis vector

  # Return the hypothesis matrix and vector
  return(list(
    "hypothesis_matrix" = C,
    # "C" matrix in the hypothesis H0: C * theta = zeta
    "hypothesis_vector" = zeta
  ))
}
