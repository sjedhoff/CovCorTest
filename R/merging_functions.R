TestCovariance_base <- function(X, nv = NULL, C, Xi, method, repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }

  # just one group is nv = NULL or nv has length 1
  if(is.null(nv) | (length(nv) == 1)){
    nv <- dim(X)[2]
    vX <- matrixcalc::vech(tvar(X))
    Xq <- matrix(apply(X-rowMeans(X),2,vtcrossprod),ncol=nv)
    HatCov <- tvar(Xq)
    MSrootHatCov <- MSroot(HatCov)
  }
  # multiple groups
  else{
    N  <- sum(nv)
    kappainvv <- N / nv

    Datac <- lapply(X, centering) #Centered random variables
    VarData <- lapply(X, tvar)
    vX <-  unlist(lapply(VarData, matrixcalc::vech))
    DataQ <- mapply(Qvech, Datac, nv, SIMPLIFY = FALSE)
    HatCov_list <- lapply(DataQ, tvar)

    MSrootHatCov <- lapply(HatCov_list, MSroot)
    HatCov <- WDirect.sumL(HatCov_list, kappainvv)
  }

  if(method == "MC"){
    ResamplingResult <- ATSwS(QF(C, HatCov), repetitions)
  }
  else if(method == "BT"){
    ResamplingResult <- sapply(1:repetitions, Bootstrap, nv, MSrootHatCov)
  }
  else{
    stop("method must be 'MC' or 'BT'")
  }

  Teststatistic <- ATS(sum(nv), vX, C, HatCov, Xi)
  pvalue <- mean(ResamplingResult < Teststatistic)

  return(list("pvalue"=pvalue, "Teststatistic"=Teststatistic, "CovarianceMatrix"=HatCov))

}


TestCovariance_simple <- function(X, nv = NULL, hypothesis, method = "MC",
                                  repetitions = 1000, seed = NULL){

  X <- Listcheck(X,nv)
  # multiple groups
  if(!is.null(nv)){
    dimensions <- sapply(X, dim)[1,]
    if(max(dimensions) != mean(dimensions)){
      stop("dimensions do not accord")
    }
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

  if(!(hypothesis %in% c("equal", "equal-trace", "equal-diagonals"))){
    stop("no predefined hypothesis")
  }

  if(hypothesis == "equal"){
    if(groups == 1){
      C <- diag(1, p, p)[a,]
      Xi <- rep(0, d)
    }
    else{
      C <- Pd(groups) %x% diag(1, p, p)
      Xi <- rep(0, p*d)
    }
    return(TestCovariance_base(X, nv = nv, C = C, Xi = Xi, method = method,
                               repetitions = repetitions, seed = seed))


  }

}
