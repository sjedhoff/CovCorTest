

#' @exportS3Method
print.CovTest <- function(x, ...){
  method_print <- ifelse(x$method == "MC", "Monte-Carlo-technique", "Bootstrap")
  group_text <- ifelse(length(x$nv) == 1, "one group", paste0("",length(x$nv), "  groups"))

  cat("\n
       \t Covariance Test \n \t    ",group_text,"\n\n Hypothesis: \t\t",
  x$hypothesis,
  "\n Teststatistic value: \t",
  round(x$Teststatistic, digits = 4),
  "\n p-value: \t \t",
  round(x$pvalue, digits = 4),
  "\n \n p-value computed using ", method_print, " and n=", x$repetitions, "\n", sep = "")
}
