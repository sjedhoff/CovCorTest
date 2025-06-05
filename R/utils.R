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
Listcheck <- function(X, nv){
  #List containing elements which are no matrices
  if(is.list(X) && (min(sapply(X,is.matrix))==0)){stop("all list elements have
                                                        to be matrices")}
  # no nv
  if(is.null(nv) || (length(nv) == 1)){
    # one group
    if(is.matrix(X)){
      if(!is.null(nv) && (nv != ncol(X))){
        warning(paste0("the number of columns of X (", ncol(X), ") and the group
                       size (", nv, ") do not allign"))
      }
      data <- X
      nv_ <- NULL
    }
    # list with one element
    else{
      if(is.list(X) && (length(X) == 1)){
        data <- X[[1]]
        if(!is.null(nv) && (nv != ncol(data))){
          warning(paste0("the number of columns of X (", ncol(data), ") and the
                         group size (", nv, ") do not allign"))
        }
        nv_ <- NULL
      }
      # list with more elements


      else{
        if(is.list(X) && (length(X) > 1)){
          data <- X
          nv_ <- unlist(lapply(X, ncol))
          warning(paste0("no nv or unfitting nv is given, will procede with nv =
                       c(",paste0(nv_, collapse = " "),")"))
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
          warning(paste0("the number of columns of X (", ncol(data), ") and
                         the group size (",paste(nv, collapse = ","), ") do not
                         allign"))
        }
        nv_ <- NULL
      }
      else{
        data <- X
        nv_ <- unlist(lapply(X, ncol))
        if(!identical(as.numeric(nv), as.numeric(nv_))){
          warning(paste0("nv does not have the corresponding dimensions to X,
                         will procede with nv = c(", paste0(nv_, collapse = " ")
                         ,")"))
        }
      }

    }
    # matrix
    else{
      if(ncol(X) != sum(nv)){
        stop(paste0("the number of columns (", ncol(X),") and the sum of group
                    sizes (", sum(nv),") do not allign"))
      }
      v <- cumsum(c(1,nv))
      data <- list()
      for(i in seq_along(nv)){
        data[[i]] <- matrix(X[,v[i]:(v[i+1]-1)], ncol=nv[i])
      }
      nv_ <- nv
    }
  }




  ## Check for missing values
  # one group
  if(is.null(nv_) && any(is.na(data))){
    data_na <- data
    # remove rows with NA all the way
    data <- data[!apply(data, 1, function(x) all(is.na(x))), , drop = FALSE]
    if(nrow(data) < nrow(data_na)){
      warning(paste0(nrow(data_na) - nrow(data), " row(s) with only NA values
                       were removed"))
    }
    # remove columns with at least one NA
    data <- data[, !apply(data, 2, function(x) any(is.na(x))), drop = FALSE]
    if(ncol(data) < ncol(data_na)){
      warning(paste0(ncol(data_na) - ncol(data), " subject(s) is/are removed
                       due to missing values"))
    }
  }
  # more groups
  if(!is.null(nv_) && any(unlist(lapply(X, function(x) any(is.na(x)))))){
    data_na <- data
    # rows, where at least in one group only missing values are pres
    na_rows <- apply(vapply(data, function(mat) apply(mat, 1, function(row)
      all(is.na(row))), FUN.VALUE = logical(nrow(X[[1]]))), 1, any)
    # remove these rows
    data <- lapply(data, function(mat) mat[!na_rows, , drop=FALSE])
    if(any(na_rows)){
      warning(paste0(sum(na_rows)," row(s) with only NA values were removed"))
    }
    # remove columns with at least one NA
    data <- lapply(data, function(mat) mat[, !apply(mat, 2,
                                                    function(col)
                                                      any(is.na(col))),
                                           drop=FALSE])
    if(sum(unlist(lapply(data_na, ncol)) - unlist(lapply(data, ncol))) > 0){
      warning(paste0(sum(unlist(lapply(data_na, ncol)) -
                           unlist(lapply(data, ncol))),
                     " subject(s) is/are removed due to missing values"))
    }
    nv_ <- unlist(lapply(data, ncol))
  }


  ## Check dimensions: multiple groups
  if(!is.null(nv_)){
    dimensions <- unlist(lapply(data, nrow))
    if(max(dimensions) != mean(dimensions)){
      stop("dimensions do not accord")
    }
    if(any(unlist(lapply(data, ncol)) == 1)){
      stop("testing covariance/correlation not possible: at least one group has
           only one subject")
    }
    if(nrow(data[[1]]) == 1){
      stop("testing covariance/correlation not possible: only one variable to
           test")
    }
  }
  else{
    if(nrow(data) == 1){
      stop("testing covariance/correlation not possible with only one subject")
    }
    if(ncol(data) == 1){
      stop("testing covariance/correlation not possible with only one variable")
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
#' @param X
#'
#' @return vector
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

