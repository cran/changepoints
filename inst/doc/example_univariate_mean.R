## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(changepoints)

## -----------------------------------------------------------------------------
delta = 5 # 2*delta represents the minimum gap between boundaries
sigma2 = 1 # error variance

set.seed(1234)
y = c(rep(0, 50), rep(2, 50), rep(0, 50), rep(-2, 50)) + rnorm(200, mean = 0, sd = sqrt(sigma2)) # univariate time series
cpt_true = c(50, 100, 150)
n = length(y) # sample size

## -----------------------------------------------------------------------------
gamma_set = c(0.01, 0.5, 1, 5, 10, 50) # a set of tuning parameters for DP
DP_result = CV.search.DP.univar(y, gamma_set, delta) # grid search through cross-validation
min_idx = which.min(DP_result$test_error) # select gamma achieves the minimum validation error
cpt_DP_hat = unlist(DP_result$cpt_hat[[min_idx]]) # estimated change points by DP
cpt_DP_hat
Hausdorff.dist(cpt_DP_hat, cpt_true)
cpt_DPlr_hat = local.refine.univar(cpt_DP_hat, y) # perform local refinement
cpt_DPlr_hat
Hausdorff.dist(cpt_DPlr_hat, cpt_true)

## -----------------------------------------------------------------------------
tau_BS = 4 # user specified threshold parameter for BS
BS_object = BS.univar(y, 1, n, delta)
BS_result = thresholdBS(BS_object, tau_BS)
BS_result
cpt_BS_hat = sort(BS_result$cpt_hat[,1]) # estimated change points by BS
cpt_BS_hat
Hausdorff.dist(cpt_BS_hat, cpt_true)
cpt_BSlr_hat = local.refine.univar(cpt_BS_hat, y) # perform local refinement
cpt_BSlr_hat
Hausdorff.dist(cpt_BSlr_hat, cpt_true)

## -----------------------------------------------------------------------------
BS_CPD_result = tuneBSunivar(BS_object, y)
BS_CPD_result$cpt
Hausdorff.dist(BS_CPD_result$cpt, cpt_true)
BS_CPD_LR = local.refine.univar(BS_CPD_result$cpt, y)
BS_CPD_LR
Hausdorff.dist(BS_CPD_LR, cpt_true)

## -----------------------------------------------------------------------------
tau_WBS = 4 # user specified threshold parameter for WBS
intervals = WBS.intervals(M = 300, lower = 1, upper = n)
WBS_object = WBS.univar(y, 1, n, intervals$Alpha, intervals$Beta, delta)
WBS_result = thresholdBS(WBS_object, tau_WBS)
WBS_result
cpt_WBS_hat = sort(WBS_result$cpt_hat[,1]) # estimated change points by WBS
cpt_WBS_hat
Hausdorff.dist(cpt_WBS_hat, cpt_true)
cpt_WBSlr_hat = local.refine.univar(cpt_WBS_hat, y) # perform local refinement
cpt_WBSlr_hat
Hausdorff.dist(cpt_WBSlr_hat, cpt_true)

## -----------------------------------------------------------------------------
WBS_CPD_result = tuneBSunivar(WBS_object, y)
WBS_CPD_result$cpt
Hausdorff.dist(WBS_CPD_result$cpt, cpt_true)
WBS_CPD_LR = local.refine.univar(WBS_CPD_result$cpt, y)
WBS_CPD_LR
Hausdorff.dist(WBS_CPD_LR, cpt_true)

