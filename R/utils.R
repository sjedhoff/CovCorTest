#' @title Centering matrix
#'
#' @description matrix Pd for testing equality of the d components of a vector
#' @param d a scalar, characterizing the matrix and set its dimension
#' @return a matrix
#'
#
#'
Pd <- function(d){
  return( diag(1,d,d) - matrix(1/d,d,d) )
}

#' @title Function to transform the data into a list, if there are not already
#'
#' @param X object that should be checked
#' @param nv
#' @return list
#'
#' @noRd
Listcheck <- function(X,nv){
  if(typeof(X) != "list"){
    v <- cumsum(c(1,nv))
    Data <- list()
    for(i in 1:length(nv)){
      Data[[i]] <- matrix(X[,v[i]:(v[i+1]-1)], ncol=nv[i])
    }
  }
  else{
    Data <- X
  }
  return(Data)
}

#' @title Matrix product, where the order is switched
#'
#' @param X,M matrices that will be multiplied
#' @return matrix
#' @noRd
#'
MprodBackward <- function(X,M){
  return(M%*%X)
}


#' @title Quadratic form for vectors and matrices
#'
#' @param A,B matrices or vectors
#' @return
#' @noRd
QF <- function(A, B){
  return( A %*% B %*% t(A) )
}

#' @title Square root of a matrix
#'
#' @param X matrix
#' @return matrix
#' @export
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
#' @export
dvech <- function(X, a, d, p, inc_diag){
  if(!is.square.matrix(X)){
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



#' @title Weighted direct sums for lists
#' @description Hereby the matrices which are part of a list are multiplied with
#' the corresponding components of a matrix w, containing the weights. These,
#' now weighted matrices are put together to one larger block-diagonal matrix.
#'
#' @param X matrix
#' @param w weight matrix
#'
#' @return matrix
#' @export
WDirect.sumL <- function(X,w){
  groups <- length(X)
  if(groups==1){
    Result <- X*w
  }
  else{
    Result <- matrixcalc::direct.sum(w[1]*X[[1]], w[2]*X[[2]])
    if(groups>2){
      for (i in 3:groups){
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
vdtcrossprod <- function(X,a,d,p){
  return(dvech(tcrossprod(X,X),a,d,p))
}

#' @title Function to center observations
#'
#' @param X matrix that will be centered
centering <- function(X){
  return(X-rowMeans(X))
}


#' @title Function to calculate variance of transposed observations
#'
#' @param X matrix
#' @return matrix
tvar <- function(X){
  return(stats::var(t(X)))
}


#' @title Auxiliary function to calculate the covariance of the vectorized correlation matrix
#'
#' @param X matrix
#' @param n number of columns
#'
#' @return matrix
Qvech <- function(X, n){
  return(matrix(apply(X,2,vtcrossprod), ncol=n))
}

