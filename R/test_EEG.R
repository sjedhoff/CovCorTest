### Test Data
data("EEG", package = "MANOVA.RM")

str(EEG)

table(EEG$sex, EEG$diagnosis)

data("EEGwide", package = "MANOVA.RM")
str(EEGwide)

table(EEGwide$sex, EEGwide$diagnosis)

vars <- colnames(EEGwide)[1:6]

var(EEGwide[,vars])
d <- 6
p <- d*(d+1)/2

X <- t(dplyr::arrange(EEGwide, sex, diagnosis)[,vars])

X_list <- list(t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "AD",vars]),
          t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "MCI",vars]),
          t(EEGwide[EEGwide$sex == "M" & EEGwide$diagnosis == "SCC",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "MCI",vars]),
          t(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "SCC",vars]))
nv <- c(12,27,20,24,30,47)

TestCovariance_simple(X = X, nv = nv, hypothesis = "uncorrelated", method = "MC",
                                   repetitions = 1000, seed = NULL)
# $pvalue
# [1] 0.954
#
# $Teststatistic
# [1] 2.938429

TestCovariance_simple(X = X_list[[1]], hypothesis = "equal")
# $pvalue
# [1] 0.956
#
# $Teststatistic
# [1] 3.518496


TestCovariance_simple(X_list[[1]], hypothesis = "given-trace", A = 3)


A <- matrix(rnorm(d*d), ncol = d, nrow = d)
TestCovariance_simple(X_list[[1]], hypothesis = "given-matrix", A = A )


TestCovariance_simple(X_list[[1]], hypothesis = "uncorrelated", A = A)


TestCovariance_structure(X = X_list[[1]], structure = "sphericity", method = "MC")

TestCovariance_structure(X = X_list, structure = "sphericity", method = "MC")

X_one <- as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
                           c("brainrate_temporal", "brainrate_frontal","brainrate_central",
                             "complexity_temporal","complexity_frontal", "complexity_central")])
TestCovariance_structure(X = X_one, structure = "autoregressive", method = "MC")
