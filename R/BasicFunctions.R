

#' ATS with weighted sum
ATSwS <- function(A, repetitions){
  Chi <- matrix(stats::rchisq(dim(A)[1] * repetitions, df = 1), ncol = repetitions)
  return(colSums(crossprod(eigen(A, only.values = 1)$value, Chi)) / sum(diag(A)))
}


#' Anova-Type-statistic
ATS <- function(N, vVarData, C, HatCov, Xi = 0){
  CDiff <-  C %*% vVarData - Xi
  statisticATS <-  N * crossprod(CDiff) / (sum(diag(QF(C, HatCov))))
  return(as.numeric(statisticATS))
}

#' Function to generate bootstrap observations
gData <- function(WSigma, nv, d){
  Data <-  WSigma %*% matrix(stats::rnorm(dim(WSigma)[1] * nv), ncol = nv)
  return(Data)
}
