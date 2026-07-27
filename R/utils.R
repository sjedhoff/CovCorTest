################################################################################
##                              Basic Functions                                #
################################################################################



#' Anova-Type-Statistic with weighted sum
#'
#' @description
#' Calculation of a Anova-Type-Statistic using
#'
#' @param A a matrix
#' @param repetitions a scalar, number of runs
#' @return a vector of the length of repetitions
#'
#' @keywords internal
ATSwS <- function(A, repetitions){
  Chi <- matrix(stats::rchisq(dim(A)[1] * repetitions, df = 1),
                ncol = repetitions)
  return(colSums(crossprod(eigen(A, only.values = 1)$value, Chi))/sum(diag(A)))
}


#' Anova-Type-statistic
#'
#' @param N number of observations
#' @param vVarData a matrix of vectorized covariance/correlation data
#' @param C hypothesis matrix for calculating the ATS
#' @param HatCov covariance matrix
#' @param Xi a vector defining together with C the investigated hypothesis
#'
#' @return a vector
#'
#'
#' @keywords internal
ATS <- function(N, vVarData, C, HatCov, Xi = 0){
  CDiff <-  C %*% vVarData - Xi
  statisticATS <-  N * crossprod(CDiff) / (sum(diag(QF(C, HatCov))))
  return(as.numeric(statisticATS))
}

#' Function to generate bootstrap observations
#'
#' @param WSigma weight matrix
#' @param nv number of observations in the groups
#'
#' @return a matrix
#'
#'
#' @keywords internal
generateData <- function(WSigma, nv){
  data <- WSigma %*% matrix(stats::rnorm(dim(WSigma)[1] * nv), ncol = nv)
  return(data)
}

################################################################################
##                            Utility Functions                               ##
################################################################################






#' @title Centering matrix
#'
#' @description matrix Pd for testing equality of the d components of a vector
#' @param d a scalar, characterizing the matrix and set its dimension
#' @return a matrix
#'
#'
#' @noRd
Pd <- function(d){
  return( diag(1,d,d) - matrix(1/d,d,d) )
}
#' @title Function to transform the data into a list, if there are not already
#'
#' @param X object that should be checked
#' @param nv number of subjects per group
#'
#' @return A list with two components
#' \item{X}{ Dataset in the right format: for a single group, a single matrix.
#' For multiple groups, a list with an element for each group containing a
#' matrix.}
#' \item{nv}{ Number of subjects per group: NV for a single group and a vector
#' for multiple groups.}
#'
#'
#' @keywords internal
Listcheck <- function(X, nv = NULL) {

  ## Helpers ---------------------------------------------------------------

  is_plain_list <- function(x) {
    is.list(x) && !is.data.frame(x)
  }

  is_numeric_matrix_or_df <- function(x) {
    if (is.matrix(x)) {
      return(is.numeric(x))
    }

    if (is.data.frame(x)) {
      return(all(vapply(x, is.numeric, logical(1))))
    }

    FALSE
  }

  transpose_element <- function(x) {
    if (is.data.frame(x)) {
      t(as.matrix(x))
    } else {
      t(x)
    }
  }

  ## Input checks ----------------------------------------------------------

  if (is_plain_list(X)) {
    if (length(X) == 0) {
      stop("X must not be an empty list.")
    }

    if (!all(vapply(X, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
      stop("All list elements have to be matrices or data frames.")
    }

    if (!all(vapply(X, is_numeric_matrix_or_df, logical(1)))) {
      stop("All list elements have to be numeric matrices or numeric data frames.")
    }
  } else {
    if (!is.matrix(X) && !is.data.frame(X)) {
      stop("X has to be a numeric matrix, numeric data frame, or a list of those.")
    }

    if (!is_numeric_matrix_or_df(X)) {
      stop("X has to be a numeric matrix or numeric data frame.")
    }
  }

  if (!is.null(nv)) {
    if (!is.numeric(nv) || any(is.na(nv)) || any(nv <= 0) || any(nv != floor(nv))) {
      stop("nv has to contain positive integer group sizes.")
    }

    nv <- as.integer(nv)
  }

  ## Construct data --------------------------------------------------------

  # Case 1: no nv or a single nv value
  if (is.null(nv) || length(nv) == 1) {

    # Single matrix/data frame input
    if (!is_plain_list(X)) {
      if (!is.null(nv) && nv != nrow(X)) {
        warning(paste0(
          "the number of rows of X (", nrow(X),
          ") and the group size (", nv, ") do not align"
        ))
      }

      data <- transpose_element(X)
      nv_ <- NULL
    }

    # List input
    else {

      # List with one element: treat as one group
      if (length(X) == 1) {
        data_raw <- X[[1]]

        if (!is.null(nv) && nv != nrow(data_raw)) {
          warning(paste0(
            "the number of rows of X (", nrow(data_raw),
            ") and the group size (", nv, ") do not align"
          ))
        }

        data <- transpose_element(data_raw)
        nv_ <- NULL
      }

      # List with multiple elements: infer nv from list elements
      else {
        data <- lapply(X, transpose_element)
        nv_ <- vapply(data, ncol, integer(1))

        warning(paste0(
          "no nv or unfitting nv is given, will proceed with nv = c(",
          paste(nv_, collapse = " "), ")"
        ))
      }
    }
  }

  # Case 2: nv with multiple group sizes
  else {

    # List input
    if (is_plain_list(X)) {

      # List with one element: treat as one group and warn if nv has multiple entries
      if (length(X) == 1) {
        data <- transpose_element(X[[1]])

        warning(paste0(
          "X contains only one group but nv has length ", length(nv),
          "; will proceed with one group"
        ))

        nv_ <- NULL
      }

      # List with multiple elements
      else {
        data <- lapply(X, transpose_element)
        nv_ <- vapply(data, ncol, integer(1))

        if (!identical(as.numeric(nv), as.numeric(nv_))) {
          warning(paste0(
            "nv does not have the corresponding dimensions to X, will proceed with nv = c(",
            paste(nv_, collapse = " "), ")"
          ))
        }
      }
    }

    # Matrix/data frame input: split rows according to nv
    else {
      X_mat <- if (is.data.frame(X)) as.matrix(X) else X

      if (nrow(X_mat) != sum(nv)) {
        stop(paste0(
          "the number of rows (", nrow(X_mat),
          ") and the sum of group sizes (", sum(nv), ") do not align"
        ))
      }

      starts <- cumsum(c(1, nv[-length(nv)]))
      ends <- cumsum(nv)

      data <- vector("list", length(nv))

      for (i in seq_along(nv)) {
        data[[i]] <- t(X_mat[starts[i]:ends[i], , drop = FALSE])
      }

      nv_ <- nv
    }
  }

  ## Check dimensions before missing-value handling ------------------------

  if (!is.null(nv_)) {
    dimensions <- vapply(data, nrow, integer(1))

    if (length(unique(dimensions)) != 1) {
      stop("All groups must have the same number of variables.")
    }
  }

  ## Check for missing values ----------------------------------------------

  # One group
  if (is.null(nv_) && any(is.na(data))) {
    data_na <- data

    # Remove rows where all values are NA
    keep_rows <- !apply(data, 1, function(x) all(is.na(x)))
    data <- data[keep_rows, , drop = FALSE]

    if (nrow(data) < nrow(data_na)) {
      warning(paste0(
        nrow(data_na) - nrow(data),
        " row(s) with only NA values were removed"
      ))
    }

    data_after_rows <- data

    # Remove columns with at least one NA
    keep_cols <- !apply(data, 2, function(x) any(is.na(x)))
    data <- data[, keep_cols, drop = FALSE]

    if (ncol(data) < ncol(data_after_rows)) {
      warning(paste0(
        ncol(data_after_rows) - ncol(data),
        " subject(s) is/are removed due to missing values"
      ))
    }
  }

  # Multiple groups
  if (!is.null(nv_) && any(vapply(data, function(x) any(is.na(x)), logical(1)))) {
    data_na <- data

    # Remove variable rows where at least one group has only NA values
    na_rows_by_group <- vapply(
      data,
      function(mat) apply(mat, 1, function(row) all(is.na(row))),
      FUN.VALUE = logical(nrow(data[[1]]))
    )

    na_rows <- apply(na_rows_by_group, 1, any)

    data <- lapply(data, function(mat) mat[!na_rows, , drop = FALSE])

    if (any(na_rows)) {
      warning(paste0(
        sum(na_rows),
        " row(s) with only NA values were removed"
      ))
    }

    # Remove subjects/columns with at least one NA within each group
    data <- lapply(data, function(mat) {
      mat[, !apply(mat, 2, function(col) any(is.na(col))), drop = FALSE]
    })

    removed_subjects <- sum(
      vapply(data_na, ncol, integer(1)) -
        vapply(data, ncol, integer(1))
    )

    if (removed_subjects > 0) {
      warning(paste0(
        removed_subjects,
        " subject(s) is/are removed due to missing values"
      ))
    }

    nv_ <- vapply(data, ncol, integer(1))
  }

  ## Final dimension checks ------------------------------------------------

  # Multiple groups
  if (!is.null(nv_)) {
    dimensions <- vapply(data, nrow, integer(1))

    if (length(unique(dimensions)) != 1) {
      stop("All groups must have the same number of variables.")
    }

    if (any(vapply(data, ncol, integer(1)) == 1)) {
      stop("testing covariance/correlation not possible: at least one group has only one subject")
    }

    if (nrow(data[[1]]) == 1) {
      stop("testing covariance/correlation not possible: only one variable to test")
    }
  }

  # One group
  else {
    if (nrow(data) == 1) {
      stop("testing covariance/correlation not possible with only one variable")
    }

    if (ncol(data) == 1) {
      stop("testing covariance/correlation not possible with only one subject")
    }
  }

  return(list(X = data, nv = nv_))
}

#' @title Quadratic form for vectors and matrices
#'
#' @param A,B matrices or vectors
#' @return a metrix or vector
#'
#' @noRd
QF <- function(A, B){
  return( A %*% B %*% t(A) )
}

#' @title Square root of a matrix
#'
#' @param X matrix
#' @return matrix
#'
#' @noRd
MSroot <- function(X){
  if(length(X) == 1){
    MSroot <- matrix(sqrt(X),1,1)
  }
  else{
    SVD <- svd(X)
    MSroot <- SVD$u %*% ( tcrossprod(sqrt(diag(SVD$d)), (SVD$v)) )
  }
  return(MSroot)
}


#' @title Diagonal vectorization
#' @description  Diagonal vectorization of the upper triangular matrix
#' @param X quadratic matrix which should be diagonalized
#' @param a vector containing the indices which belong to the diagonal of the
#' matrix
#' @param d dimension of the matrix which should be vectorized
#' @param p dimension of the vectorized  matrix
#' @param inc_diag TRUE or FALSE: should the diagonal be included?
#'
#' @return vector
#'
#' @keywords internal
#'
dvech <- function(X, a, d, p, inc_diag){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }
  else{
    E <- rep(X[1,d],p)
    for(i in 1:(d-1)){
      E[a[i]:(a[i+1]-1)] <- diag(X[1:(d-i+1),i:d])
    }
    # without the diagonal
    if(!inc_diag){
      E <- E[-(1:d)]
    }

    return(E)
  }
}

#' Vectorization of the upper triangular part of the matrix
#'
#' @param X A square numeric matrix.
#'
#' @return A vector containing the upper triangular part of \code{X}.
#'
#' @keywords internal
#'
vechp <- function(X){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }

  return(as.vector(t(X)[!upper.tri(X,TRUE)]))
}

#' @title Weighted direct sums for lists
#' the corresponding components of a matrix w, containing the weights. These,
#' now weighted matrices are put together to one larger block-diagonal matrix.
#'
#' @param X matrix
#' @param w weight matrix
#'
#' @return matrix
#'
#' @keywords internal
#'
WDirect.sumL <- function(X, w){
  groups <- length(X)
  if(groups == 1){
    Result <- X*w
  }
  else{
    Result <- matrixcalc::direct.sum(w[1]*X[[1]], w[2]*X[[2]])
    if(groups > 2){
      for(i in 3:groups){
        Result <-  matrixcalc::direct.sum(Result, w[i]*X[[i]])
      }
    }
  }
  return(Result)
}



#' @title Function to calculate vech(X t(X))
#'
#' @param X matrix
#' @return vector
#'
#'
#' @keywords internal
vtcrossprod <- function(X){
  return(matrixcalc::vech(tcrossprod(X,X)))
}

#' @title Function to calculate dvech(X t(X))
#'
#' @param X matrix
#' @param a indices that belong to the diagonal of the matrix
#' @param d dimension of the matrix
#' @param p dimension of the vectorized  matrix
#'
#' @return vector
#'
#'
#' @keywords internal
vdtcrossprod <- function(X,a,d,p){
  return(dvech(tcrossprod(X,X),a,d,p, inc_diag = TRUE))
}



#' @title Auxiliary function to calculate the covariance of the vectorized
#' correlation matrix
#'
#' @param X matrix
#' @param n number of columns
#'
#' @return matrix
#'
#'
#' @keywords internal
Qvech <- function(X, n){
  return(matrix(apply(X,2,vtcrossprod), ncol=n))
}

#' Compact matrix square root
#'
#' @param X Numeric positive semidefinite matrix.
#' @param tol Numeric tolerance used for the rank calculation.
#'
#' @return A compact matrix square root L such that t(L) %*% L is approximately X.
#'
#' @noRd
.MSrootcompact <- function(X, tol = 1e-10) {

  r <- qr(X, tol = tol)$rank

  if (r == 0) {
    return(matrix(0, nrow = 0, ncol = ncol(X)))
  }

  SVD <- svd(X)

  L <- sqrt(diag(SVD$d[seq_len(r)], nrow = r, ncol = r)) %*%
    t(SVD$u[, seq_len(r), drop = FALSE])

  return(L)
}

#' Check hypothesis matrix and vector
#'
#' @param H Numeric hypothesis matrix.
#' @param y Numeric hypothesis vector.
#' @param tol Numeric tolerance.
#'
#' @return A list containing checked H and y.
#'
#' @noRd
HypoCheck <- function(H, y, tol = 1e-10) {

  if (!is.matrix(H)) {
    stop("The hypothesis must be specified by a matrix.")
  }

  if (!is.numeric(H)) {
    stop("The hypothesis matrix must be numeric.")
  }

  if (is.matrix(y) && ncol(y) == 1) {
    y <- as.numeric(y)
  }

  if (!is.vector(y) || !is.numeric(y)) {
    stop("The corresponding hypothesis vector must be numeric.")
  }

  y <- as.numeric(y)

  if (nrow(H) != length(y)) {
    stop("Dimensions of hypothesis matrix and vector must match.")
  }

  if (any(!is.finite(H)) || any(!is.finite(y))) {
    stop("Hypothesis matrix and vector must contain finite values only.")
  }

  qH <- qr(H, tol = tol)
  rH <- qH$rank

  qty <- qr.qty(qH, y)

  if (length(y) > rH && any(abs(qty[(rH + 1):length(y)]) > tol)) {
    stop("The hypothesis H v = y is inconsistent.")
  }

  return(list(H = H, y = y))
}


#' Companion hypothesis matrix
#'
#' Constructs an internal companion representation of a hypothesis matrix. If
#' the matrix has full row rank, the original matrix is returned. Otherwise, a
#' compact matrix square root of t(H) %*% H is computed, yielding a matrix with
#' fewer rows but equivalent quadratic forms. If a non-zero vector y is supplied,
#' it is transformed accordingly.
#'
#' @param H Numeric hypothesis matrix.
#' @param y Optional numeric vector corresponding to the right-hand side of the
#'   hypothesis. If omitted, a zero vector is used.
#' @param tol Numeric tolerance used for rank calculations.
#'
#' @return A list containing the matrix L and the transformed vector ytilde.
#'
#' @noRd
CompanionHypothesis <- function(H, y = NULL, tol = 1e-10) {

  if (is.null(y)) {
    y <- rep(0, nrow(H))
  }

  checked <- HypoCheck(H, y, tol = tol)
  H <- checked$H
  y <- checked$y

  if (qr(H, tol = tol)$rank == nrow(H)) {

    L <- H
    ytilde <- y

  } else {

    L <- .MSrootcompact(t(H) %*% H, tol = tol)

    if (all(y == 0)) {
      ytilde <- rep(0, nrow(L))
    } else {
      ytilde <- qr.solve(t(L), t(H) %*% y)
    }
  }

  return(list(
    L = L,
    ytilde = ytilde
  ))
}
#' Check bandwidth argument
#'
#' Checks whether the bandwidth argument is valid for banded covariance or
#' correlation structures. The bandwidth specifies the number of off-diagonals
#' that are allowed to be non-zero. It must be a positive integer and smaller
#' than the maximal possible bandwidth, since bandwidth zero corresponds to a
#' diagonal structure and the maximal bandwidth would impose no restriction.
#'
#' @param bandwidth A positive integer specifying the number of allowed
#'   off-diagonals.
#' @param d Integer. Dimension of the covariance or correlation matrix.
#'
#' @return The checked bandwidth as an integer.
CheckBandwidth <- function(bandwidth, d) {

  if (length(bandwidth) != 1) {
    stop("bandwidth must be a single value.")
  }

  if (is.na(bandwidth)) {
    stop("bandwidth must be specified for banded structures.")
  }

  if (!is.numeric(bandwidth) || !is.finite(bandwidth) ||
      bandwidth != floor(bandwidth)) {
    stop("bandwidth must be a positive integer.")
  }

  bandwidth <- as.integer(bandwidth)

  if (d < 3) {
    stop("A banded structure with positive bandwidth requires dimension d >= 3.")
  }

  if (bandwidth <= 0 || bandwidth >= d - 1) {
    stop("bandwidth must be between 1 and d - 2.")
  }

  return(bandwidth)
}

#' Validate the number of resampling repetitions
#'
#' Checks that the requested number of repetitions is a single, finite,
#' positive integer. A warning is issued when the number is below the
#' recommended minimum because the resulting p-value may be imprecise.
#'
#' @param repetitions A single positive integer specifying the number of
#'   resampling repetitions.
#' @param minimum_recommended A single positive integer specifying the number
#'   of repetitions below which a warning is issued. The default is 500.
#'
#' @return The validated number of repetitions as an integer.
#'
CheckRepetitions <- function(repetitions, minimum_recommended = 500L) {
  if (
    !is.numeric(repetitions) ||
    length(repetitions) != 1L ||
    is.na(repetitions) ||
    !is.finite(repetitions) ||
    repetitions <= 0 ||
    repetitions != floor(repetitions)
  ) {
    stop("repetitions must be a single positive integer.")
  }

  repetitions <- as.integer(repetitions)

  if (repetitions < minimum_recommended) {
    warning(
      paste0(
        "Fewer than ", minimum_recommended,
        " repetitions may result in an imprecise p-value."
      ),
      call. = FALSE
    )
  }

  repetitions
}
