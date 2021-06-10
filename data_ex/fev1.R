library(foreign)
library(lme4)
library(lmmstest)
library(merDeriv)

fev1_dat <- read.dta("http://www.hsph.harvard.edu/fitzmaur/ala2e/fev1.dta")

# Fit model
fit <- lmer(exp(logfev1) ~ age + ht + baseage + baseht + (1|id) +
                    (0 + age|id), data = fev1_dat, REML = FALSE)

# CI
confint(fit)
confint(fit, method = "Wald")
obs_inf <- vcov(fit, full = TRUE, information = "observed")
est <- as.data.frame(VarCorr(fit))[1:3, 4]
se <- sqrt(diag(obs_inf)[6:8])
sqrt(est + qnorm(0.975) * se)
sqrt(est - qnorm(0.975) * se)

# Scale parameters CI
# Use data resulting from running the code in the section "Get confidence
# regions" below
ex_dat <- readRDS("~/GitHub/int-est/data_ex/data_ex.Rds")

# Confidence interval for random intercept standard deviation
lam_val <- ex_dat[[1]]
test_val <- ex_dat[[4]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[30:31]
c(0, mean(lam_val[30:31]))

# Confidence interval for random slope standard deviation
lam_val <- ex_dat[[2]]
test_val <- ex_dat[[5]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[c(8:9, 39:40)]
c(mean(lam_val[8:9]), mean(lam_val[39:40]))

# Confidence interval for error term standard deviation
sig_val <- ex_dat[[3]]
test_val <- ex_dat[[7]] - qchisq(0.95, 1)
test_val
# Take mean where sign switches
test_val[c(5:6, 30:31)]
c(mean(sig_val[5:6]), mean(sig_val[30:31]))

# # Get confidence regions ------------------------------------------------------
# 
# X <- getME(fit, "X")
# Z <- as.matrix(getME(fit, "Z"))
# y <- getME(fit, "y")
# Beta_hat <- unname(fixef(fit))
# VC <- as.data.frame(VarCorr(fit))
# sig_hat <- VC[3, 5]
# lam_hat <- VC[1:2, 5]
# theta_hat <- c(Beta_hat, sig_hat, lam_hat)
# 
# # lambda_1
# get_score_stat_1 <- function(lam_null){
#   lmmstest::score_test(y = y,
#                        X = X,
#                        Z = Z,
#                        Beta = Beta_hat,
#                        sigma = sig_hat,
#                        lambda = c(lam_null, lam_hat[2]),
#                        lam_idx = rep(1:2, each = 300),
#                        test_idx = 7)["chi_sq"]
# }
# 
# lam_vec_1 <- seq(0, 0.08, length.out = 50)
# 
# test_vec_1 <- rep(0, length(lam_vec_1))
# for(ii in 1:length(lam_vec_1)){
#   test_vec_1[ii] <- get_score_stat_1(lam_vec_1[ii])
# }
# 
# # lambda_2
# get_score_stat_2 <- function(lam_null){
#   lmmstest::score_test(y = y,
#                        X = X,
#                        Z = Z,
#                        Beta = Beta_hat,
#                        sigma = sig_hat,
#                        lambda = c(lam_hat[1], lam_null),
#                        lam_idx = rep(1:2, each = 300),
#                        test_idx = 8)["chi_sq"]
# }
# 
# lam_vec_2 <- seq(0.018, 0.023, length.out = 50)
# test_vec_2 <- rep(0, length(lam_vec_2))
# for(ii in 1:length(lam_vec_2)){
#   test_vec_2[ii] <- get_score_stat_2(lam_vec_2[ii])
# }
# 
# 
# # Joint lambda
# library(foreach)
# library(doParallel)
# library(parallel)
# 
# get_score_stat_3 <- function(lam_null){
#   lmmstest::score_test(y = y,
#                        X = X,
#                        Z = Z,
#                        Beta = Beta_hat,
#                        sigma = sig_hat,
#                        lambda = lam_null,
#                        lam_idx = rep(1:2, each = 300),
#                        test_idx = 7:8)["chi_sq"]
# }
# 
# test_mat <- matrix(NA, nrow = length(lam_vec_1), ncol = length(lam_vec_2))
# 
# numCores <- detectCores() - 2
# cl <- makeCluster(numCores)
# registerDoParallel(cl)
# 
# for(ii in 1:length(lam_vec_1)){
#   cat("ii = ", ii, "\n")
#   out <- foreach(jj = 1:length(lam_vec_2), .combine = c) %dopar% {
#     res <- try(get_score_stat_3(c(lam_vec_1[ii], lam_vec_2[jj])), TRUE)
#     if(class(res) == "try-error"){
#       res <- NA
#     }
#     res
#   }
#   test_mat[ii, ] <- out
# }
# stopCluster(cl)
# 
# # sigma
# get_score_stat_4 <- function(sig_null){
#   lmmstest::score_test(y = y,
#                        X = X,
#                        Z = Z,
#                        Beta = Beta_hat,
#                        sigma = sig_null,
#                        lambda = lam_hat,
#                        lam_idx = rep(1:2, each = 300),
#                        test_idx = 6)["chi_sq"]
# }
# 
# sig_vec <- seq(0.15, 0.17, length.out = 50)
# test_vec_3 <- rep(0, length(sig_vec))
# for(ii in 1:length(sig_vec)){
#   test_vec_3[ii] <- get_score_stat_4(sig_vec[ii])
# }
# 
# data_ex_conf <- list(lam_vec_1, lam_vec_2, sig_vec, test_vec_1, test_vec_2,
#                      test_mat, test_vec_3, fit)
# saveRDS(data_ex_conf, "~/GitHub/int-est/data_ex/data_ex.Rds")
