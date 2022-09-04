## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(changepoints)

## -----------------------------------------------------------------------------
## parameters for simulating data
rho = 0.4 #candidate of rho (rho can not be too large so that the time series is stable)
p = 20 # dimension
sigma = 1 # standard deviation of error terms
delta1 = 10 # minimum gap for DP
delta.local = 10 # minimum gap for local.refine
n = 30 # length for each segment
v1 = 2*(seq(1,p,1)%%2) - 1
v2 = -v1
AA = matrix(0, nrow = p, ncol = p-2)
A1 = rho * cbind(v1,v2,AA) # transition matrix for the first segment
A2 = rho * cbind(v2,v1,AA) # transition matrix for the second segment
A3 = A1 # transition matrix for the third segment

set.seed(123)
# generate data
data = simu.VAR1(sigma, p, 2*n+1, A1)
data = cbind(data, simu.VAR1(sigma, p, 2*n, A2, vzero=c(data[,ncol(data)])))
data = cbind(data, simu.VAR1(sigma, p, 2*n, A3, vzero=c(data[,ncol(data)])))
cpt_true = c(2*n, 4*n)

## -----------------------------------------------------------------------------
gamma_set = c(0.1, 0.5, 1, 2, 5) # a set of tuning parameters for DP
lambda_set = c(0.1, 0.5, 1, 2, 5) # a set of tuning parameters for lasso
DP_result = CV.search.DP.VAR1(data, gamma_set, lambda_set, delta1) # grid search through cross-validation
min_idx = as.vector(arrayInd(which.min(DP_result$test_error), dim(DP_result$test_error)))# select gamma achieves the minimum validation error
cpt_DP_hat = unlist(DP_result$cpt_hat[min_idx[1], min_idx[2]]) # estimated changepoints by DP
cpt_DP_hat
Hausdorff.dist(cpt_DP_hat, cpt_true)
zeta_set = c(1, 1.5, 2) # tuning parameter for group lasso
cpt_DPlr_hat = local.refine.CV.VAR1(cpt_DP_hat, data, zeta_set, delta.local) # perform local refinement
cpt_DPlr_hat
Hausdorff.dist(cpt_DPlr_hat$cpt_hat, cpt_true)

