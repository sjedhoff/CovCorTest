## bekannt aus dem Teil zu linearen Modellen:
Baum.lm <- lm(Volume ~ Girth, data = datasets::trees)
## neues Objekt (und Objektklasse) konstruieren:
object <- list(data = datasets::trees[, c("Girth" , "Volume" )],
               lm.obj = Baum.lm)
class(object) <- "nurSoEinName"
print.nurSoEinName <- function(x, ...){
  cat(" Data:\n=====\n\n" )
  print(x$data, ...)
  cat(" \n\nEstimated model:\n================\n" )
  print(x$lm.obj, ...)
  invisible(x)
}
object

print.CovTest <- function(x, ...){
  cat("Hypothesis:")
  print(x$hypothesis)
  cat("p-Value:")
  print(x$pvalue)
}

x <- TestCovariance_simple(X_list[[1]], hypothesis = "uncorrelated")

x_test <- t.test(rnorm(10))
str(x_test)
