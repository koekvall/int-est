library(lme4)
# Settings are (n, T, sigma_e, sigma_u1, sigma_u2)
settings <- c(40, 10, 1.5, 0.1, 0.5)

# Settings
n_i <- settings[1]
n_t <- settings[2]
sigma_e <- settings[3]
sigma_u1 <-settings[4]
sigma_u2 <- settings[5]

X <- matrix(runif(n_i * n_t, -1, 2) , nrow = n_i, ncol = n_t)
U1 <- rnorm(n_i, sd = sigma_u1)
U2 <- rnorm(n_i, sd = sigma_u2)
E <- matrix(rnorm(n_i * n_t, sd = sigma_e), nrow = n_i, ncol = n_t)
Y <- X + 1 + E
for(ii in 1:n_i){
  Y[ii, ] <- Y[ii, ] + U1[ii] + U2[ii] * X[ii, ]
}
D <- data.frame(c(t(Y)), c(t(X)))
D$time <- rep(1:n_t, n_i)
D$obs <- as.factor(rep(1:n_i, each = n_t))
names(D)[1:2] <- c("y", "x")
fit <- lmer(y ~ x + (1|obs) + (0 + x|obs), data = D, REML = FALSE)


# Test correctness of output
X_tall <- getME(fit, "X")
Z_tall <- as.matrix(getME(fit,"Z"))
n <- nrow(X_tall)
Beta_test <- c(1, 1)
sigma_test <- sigma_e
lambda_test <- c(sigma_u1, sigma_u2)
Lam_test <- diag(rep(lambda_test, each = n_i), n_i * 2)
Sigma_test <- diag(sigma_test^2, n_i * n_t) + Z_tall %*% Lam_test^2 %*% t(Z_tall)

out <- lmmstest::log_lik(y = D$y,
                        X = X_tall,
                        Z = Z_tall,
                        Beta = Beta_test,
                        sigma = sigma_test,
                        lambda = lambda_test,
                        lam_idx = rep(1:2, each = n_i),
                        diffs = 2)
# Value test
out[[1]] - mvtnorm::dmvnorm(x = matrix(D$y - X_tall %*% Beta_test, nrow = 1),
                            sigma = Sigma_test, log = T)
# Score test
test_fun <- function(theta){
  lmmstest::log_lik(y = D$y,
                   X = X_tall,
                   Z = Z_tall,
                   Beta = theta[1:2],
                   sigma = theta[3],
                   lambda = theta[4:5],
                   lam_idx = rep(1:2, each = n_i),
                   diffs = 0)[[1]]
}

numDeriv::grad(test_fun, x = c(Beta_test, sigma_test, lambda_test)) - out[[2]] * c(1, 1, sigma_test, lambda_test)
numDeriv::grad(test_fun, x = c(Beta_test, sigma_test, lambda_test)) - 
  lmmstest::score(y = D$y,
       X = X_tall,
       Z = Z_tall,
       Beta = Beta_test,
       sigma = sigma_test,
       lambda = lambda_test,
       lam_idx = rep(1:2, each = n_i))

# Simulation to test Fisher information
library(lme4)
# Settings are (n, T, sigma_u1, sigma_u2, sigma_e)
do_one_sim <- function(seed, settings){
  set.seed(seed)
  
  n_i <- settings[1]
  n_t <- settings[2]
  sigma_e <- settings[3]
  sigma_u1 <-settings[4]
  sigma_u2 <- settings[5]
  
  # Generate data
  X <- matrix(runif(n_i * n_t, -1, 2) , nrow = n_i, ncol = n_t)
  U1 <- rnorm(n_i, sd = sigma_u1)
  U2 <- rnorm(n_i, sd = sigma_u2)
  E <- matrix(rnorm(n_i * n_t, sd = sigma_e), nrow = n_i, ncol = n_t)
  Y <- X + 1 + E
  for(ii in 1:n_i){
    Y[ii, ] <- Y[ii, ] + U1[ii] + U2[ii] * X[ii, ]
  }
  D <- data.frame(c(t(Y)), c(t(X)))
  D$time <- rep(1:n_t, n_i)
  D$obs <- as.factor(rep(1:n_i, each = n_t))
  names(D)[1:2] <- c("y", "x")
  
  # Fit model
  fit <- lme4::lmer(y ~ x + (1|obs) + (0 + x|obs), data = D, REML = FALSE)
  comp <- lmmstest::log_lik(
    y = D$y,
    X = lme4::getME(fit, "X"),
    Z = as.matrix(lme4::getME(fit, "Z")),
    Beta = c(1, 1),
    sigma = sigma_e,
    lambda = c(sigma_u1, sigma_u2),
    lam_idx = rep(1:2, each = n_i),
    diffs = 2
  )
  cbind(comp[[2]], comp[[3]])
}

# Do simulation
library(doParallel)
library(doRNG)
cl <- makeCluster(6)
registerDoParallel(cl)
# Settings are (n, T, sigma_e, sigma_u1, sigma_u2)
settings <- c(40, 10, 1, 0.1, 0.5)

res <- foreach(ii = 1:1e4,
                   .errorhandling = "remove",
                   .packages = c("lmmstest", "lme4")) %dorng%{
                     do_one_sim(ii, settings)
                   }
# Check mean of scaled score
colMeans(do.call(rbind, lapply(res, function(x)x[, 1])))

# Check covariance of Score and Fisher information
emp_score_cov <- cov(do.call(rbind, lapply(res, function(x)x[, 1])))
emp_mean_inf <- Reduce("+", lapply(res, function(x)x[, 2:6])) / length(res)
(emp_score_cov - emp_mean_inf) / emp_score_cov
