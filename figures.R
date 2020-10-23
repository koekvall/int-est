out_dir <- "~/GitHub/int_est/" # end in "/"
PDF <- TRUE

ex_1_ll <- function(y, lam)
{
  sum(- 0.5 * log(1 + lam^2) - 0.5 * y^2 / (1 + lam^2))
}

ex_1_s <- function(y, lam)
{
  sum(-lam / (1 + lam^2) + y^2 * lam / (1 + lam^2)^2)
}

ex_1_h <- function(y, lam)
{
  sum(y^2 * (1 - 3 * lam^2) / (1 + lam^2)^3 - (1 - lam^2) / (1 + lam^2)^2)
}

ex_1_stat <- function(y, lam)
{
  sum(-1 + y^2 / (1 + lam^2))^2 / (2 * length(y))
}

ex_1_obs_stat <- function(y, lam)
{
  ex_1_s(y, lam)^2 / (-ex_1_h(y, lam))
}

set.seed(3)

# 1
y1 <- rnorm(100, sd = 1.1)
y2 <- rnorm(100, sd = 1.1)
y3 <- rnorm(100, sd = 1.1)
lam_vals <- seq(-0.1, 1.5, length.out = 10000)
ll_vals1 <- sapply(lam_vals, function(x)ex_1_ll(y1, x))
ll_vals2 <- sapply(lam_vals, function(x)ex_1_ll(y2, x))
# ll_vals3 <- sapply(lam_vals, function(x)ex_1_ll(y3, x))
all_vals <- c(ll_vals1, ll_vals2)


if(PDF) pdf(paste(out_dir, "fig_1.pdf", sep = ""), width = 12.5, height = 5)

par(cex.axis = 1.3, cex.lab = 1.3, cex = 1.3)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 3))

# Likelihood
plot(lam_vals, ll_vals1, type = "l", lwd = 3, ylim = c(min(all_vals), max(all_vals)),
     ylab = "log-likelihood", xlab = expression(lambda))
lines(lam_vals, ll_vals2, type = "l", lwd = 3, lty = 1)



# Test statistics 2
our_vals <- sapply(lam_vals, function(x)ex_1_stat(y1, x))
obs_vals <-  sapply(lam_vals, function(x)ex_1_obs_stat(y1, x))
plot(lam_vals, our_vals, type = "l", lwd = 3,
     ylab = "value of test statistic", xlab = expression(lambda),
     ylim = c(0, 15), lty = 3)
lines(lam_vals, obs_vals, type = "l", lwd = 3, lty = 2)

# Test statistics 1
our_vals <- sapply(lam_vals, function(x)ex_1_stat(y2, x))
obs_vals <-  sapply(lam_vals, function(x)ex_1_obs_stat(y2, x))
plot(lam_vals, our_vals, type = "l", lwd = 3,
     ylab = "value of test statistic", xlab = expression(lambda),
     ylim = c(-15, 15), lty = 3)
lines(lam_vals, obs_vals, type = "l", lwd = 3, lty = 2)
if(PDF) dev.off()

# Test derivatives
numDeriv::grad(function(x)ex_1_ll(y2, x), 0.2)
numDeriv::hessian(function(x)ex_1_ll(y2, x), 0.2)
ex_1_s(y2, 0.2)
ex_1_h(y2, 0.2)
