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
#' @return list
#'
#' @export
#' @keywords internal
Listcheck <- function(X, nv){
  # no nv
  if(is.null(nv) | (length(nv) == 1)){
    # one group
    if(is.matrix(X)){
      if(!is.null(nv) && (nv != ncol(X))){
        warning(paste0("the number of columns of X (", ncol(X), ") and the group size (", nv, ") do not allign"))
      }
      data <- X
      nv_ <- NULL
    }
    # list with one element
    else{
      if(is.list(X) & (length(X) == 1)){
        data <- X[[1]]
        if(!is.null(nv) && (nv != ncol(data))){
          warning(paste0("the number of columns of X (", ncol(data), ") and the group size (", nv, ") do not allign"))
        }
        nv_ <- NULL
      }
      # list with more elements
      else{if(is.list(X) & (length(X) > 1)){
        data <- X
        nv_ <- unlist(lapply(X, ncol))
        warning(paste0("no nv or unfitting nv is given, will procede with nv = c(",paste0(nv_, collapse = " "),")"))
      }}
    }

  }

  # with nv
  else{
    # list
    if(is.list(X)){
      if(length(X) == 1){
        data <- X[[1]]
        if((length(nv) != 1) || (nv != ncol(data))){
          warning(paste0("the number of columns of X (", ncol(data), ") and the group size (",paste(nv, collapse = ","), ") do not allign"))
        }
        nv_ <- NULL
      }
      else{
        data <- X
        nv_ <- unlist(lapply(X, ncol))
        if(!identical(as.numeric(nv), as.numeric(nv_))){
          warning(paste0("nv does not have the corresponding dimensions to X, will procede with nv = c(", paste0(nv_, collapse = " "),")"))
        }
      }

    }
    # matrix
    else{
      if(ncol(X) != sum(nv)){
        stop(paste0("the number of columns (", ncol(X),") and the sum of group sizes (", sum(nv),") do not allign"))
      }
      v <- cumsum(c(1,nv))
      data <- list()
      for(i in 1:length(nv)){
        data[[i]] <- matrix(X[,v[i]:(v[i+1]-1)], ncol=nv[i])
      }
      nv_ <- nv
    }
  }

  ## For multiple groups: check dimensions
  if(!is.null(nv_)){
    dimensions <- unlist(lapply(data, nrow))
    if(max(dimensions) != mean(dimensions)){
      stop("dimensions do not accord")
    }
  }

  return(list(data,nv_))
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



#' @title Diagonal vectorisation
#' @description  Diagonal vectorisation of the upper triangular matrix
#' @param X quadratic matrix which should be diagonalized
#' @param a vector containing the indices which belong to the diagonal of the
#' matrix
#' @param d dimension of the matrix which should be vectorised
#' @param p dimension of the vectorised matrix
#' @param inc_diag TRUE or FALSE: should the diagonal be included?
#'
#' @return vector
#'
#' @keywords internal
#' @export
dvech <- function(X, a, d, p, inc_diag){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }
  else{
    E <- rep(X[1,d],p)
    if(inc_diag){ # with the diagonal
      for(i in 1:(d-1)){
        E[a[i]:(a[i+1]-1)] <- diag(X[1:(d-i+1),i:d])
      }
    }
    else{ # without the diagonal
      for(i in 2:(d-1)){
        E[(a[i]-d):(a[i+1]-1-d)] <- diag(X[1:(d-i+1),i:d])
      }
    }
    return(E)
  }
}

#' Vectorization of the upper triangular part of the matrix
#'
#' @param X
#'
#' @return vector
#'
#' @keywords internal
#' @export
vechp <- function(X){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }

  return(as.vector(t(X)[!upper.tri(X,TRUE)]))
}


#' @title Weighted direct sums for lists
#' @description Hereby the matrices which are part of a list are multiplied with
#' the corresponding components of a matrix w, containing the weights. These,
#' now weighted matrices are put together to one larger block-diagonal matrix.
#'
#' @param X matrix
#' @param w weight matrix
#'
#' @return matrix
#'
#' @keywords internal
#' @export
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
#' @export
#' @keywords internal
vtcrossprod <- function(X){
  return(matrixcalc::vech(tcrossprod(X,X)))
}

#' @title Function to calculate dvech(X t(X))
#'
#' @param X matrix
#' @param a indices that belong to the diagonal of the matrix
#' @param d dimension of the matrix
#' @param p dimension of the vectorised matrix
#'
#' @return vector
#'
#' @export
#' @keywords internal
vdtcrossprod <- function(X,a,d,p){
  return(dvech(tcrossprod(X,X),a,d,p, inc_diag = TRUE))
}



#' @title Auxiliary function to calculate the covariance of the vectorized correlation matrix
#'
#' @param X matrix
#' @param n number of columns
#'
#' @return matrix
#'
#' @export
#' @keywords internal
Qvech <- function(X, n){
  return(matrix(apply(X,2,vtcrossprod), ncol=n))
}

