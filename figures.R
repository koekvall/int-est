out_dir <- "/Users/karekv/Dropbox/Apps/Overleaf/int_est/" # end in "/"
PDF <- TRUE


### Introduction ###
get_statistics <- function(y, lam)
{
  n <- length(y)
  ll_fun <- function(y, lam){
    sum(- 0.5 * log(1 + lam^2) - 0.5 * y^2 / (1 + lam^2))
  }
  ll <- ll_fun(y, lam)
  s <- sum(-lam / (1 + lam^2) + y^2 * lam / (1 + lam^2)^2)
  h <- sum(y^2 * (1 - 3 * lam^2) / (1 + lam^2)^3 - (1 - lam^2) / (1 + lam^2)^2)
  finf <- n * 2 * lam^2 / (1 + lam^2)^2
  our_stat <- sum(-1 + y^2 / (1 + lam^2))^2 / (2 * n)

  score_obs_stat <- s^2 / (-h)

  mle <- sqrt(max(c(mean(y^2) - 1, 0)))
  lrt_stat <- 2 * (ll_fun(y, mle) - ll)
  wald_stat <- (mle - lam)^2 * finf
  wald_obs_stat <-  (mle - lam)^2 * (-h)
  return(c(ll = ll, s = s, h = h, finf = finf, our_stat = our_stat,
              score_obs_stat = score_obs_stat, lrt_stat = lrt_stat,
              wald_stat = wald_stat, wald_obs_stat = wald_obs_stat))
}

set.seed(3)
y1 <- rnorm(100, sd = 1.1)
y2 <- rnorm(100, sd = 1.1)
lam_vals <- seq(-0.1, 1, length.out = 1000)
dat1 <- t(sapply(lam_vals, function(x)get_statistics(y1, x)))
dat2 <-  t(sapply(lam_vals, function(x)get_statistics(y2, x)))


# Create figure
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if(PDF) pdf(paste(out_dir, "fig_1.pdf", sep = ""), width = 6, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))

# Likelihood
plot(lam_vals, dat1[, "ll"], type = "l", lwd = 2,
     ylab = "log-likelihood", xlab = expression(lambda) ,
     ylim = c(min(dat2[, "ll"]), max(dat1[, "ll"])),
     col = cbbPalette[1])
lines(lam_vals, dat2[, "ll"], type = "l", lwd = 2, lty = 1,
      col = cbbPalette[1])

if(PDF) dev.off()

### Simulation Figures ####
sim_data <- readRDS("~/GitHub/int-est/sims/sims_big.Rds")
n_sims <- nrow(sim_data[[1]])
# Verify this using names(sim_data)
gamma <- c(0, 0.01, 0.05, seq(0.1, 0.5, length.out = 5))
n_vec <- c(20, 40, 80)

# Cover curves
cover <- matrix(0, nrow = length(sim_data), ncol = 8)
for(ii in 1:nrow(cover)){
  cover[ii, 1:3] <- 1 - colMeans(sim_data[[ii]][, c("p_val", "lrt_p_val",
                                                    "wald_p_val")] < 0.05)
  cover[ii, 4:6] <- apply(sim_data[[ii]][, c("p_val", "lrt_p_val",
                                             "wald_p_val")] < 0.05,
                          2, sd) / sqrt(n_sims)
}
cover[, 7] <- rep(n_vec, each = length(gamma))
cover[, 8] <- rep(gamma, length(n_vec))
colnames(cover) <- c("our", "lrt", "wald", "se_our", "se_lrt", "se_wald",
                     "n", "gamma")

if(PDF) pdf(paste(out_dir, "fig_2.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

for(jj in c(20, 80)){
  plot_dat <- cover[cover[, "n"] == jj, ]
  plot(x = gamma,
       y = plot_dat[, "our"],
       type = "l",
       ylim = c(0.82, 1),
       lwd = 2,
       ylab = "estimated cover probability",
       xlab = expression(lambda),
       main = paste0("n = ", jj),
       col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] + 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  arrows(gamma,
         plot_dat[, "our"],
         gamma,
         plot_dat[, "our"] - 2 * plot_dat[, "se_our"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[1])
  abline(h = 0.95, lwd = 1)
  lines(x = gamma,
        y = plot_dat[, "lrt"],
        lty = 2,
        lwd = 2,
        col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] + 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  arrows(gamma,
         plot_dat[, "lrt"],
         gamma,
         plot_dat[, "lrt"] - 2 * plot_dat[, "se_lrt"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[2])
  lines(x = gamma,
        y = plot_dat[, "wald"],
        lty = 3,
        lwd = 2,
        col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] + 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
  arrows(gamma,
         plot_dat[, "wald"],
         gamma,
         plot_dat[, "wald"] - 2 * plot_dat[, "se_wald"],
         length = 0.05,
         angle = 90,
         col = cbbPalette[3])
}
if(PDF) dev.off()
# Distribution
pp <- ppoints(n_sims)
if(PDF) pdf(paste(out_dir, "fig_3.pdf", sep = ""), width = 10, height = 4)
par(cex.axis = 1.1, cex.lab = 1.1, cex = 1.1)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.5,1.1))
par(mfrow = c(1, 2))

# Small but non-zero scale parameters
plot_dat <- sim_data[[18]]
plot(x = qchisq(pp, df = 3),
     y = quantile(plot_dat[, "chi_sq"], pp),
     xlab = "theoretical quantiles",
     ylab = "sample quantiles",
     main = "small scale parameters",
     col = cbbPalette[1])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "lrt_chi_sq"], pp),
       pch = 2,
       col = cbbPalette[2])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "wald_chi_sq"], pp),
       pch = 3,
       col = cbbPalette[3])
abline(a = 0, b = 1, col = "red")

# Large scale parameters
plot_dat <- sim_data[[24]]
plot(x = qchisq(pp, df = 3),
     y = quantile(plot_dat[, "chi_sq"], pp),
     xlab = "theoretical quantiles",
     ylab = "sample quantiles",
     main = "large scale parameters",
     col = cbbPalette[1])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "lrt_chi_sq"], pp),
       pch = 2,
       col = cbbPalette[2])
points(x = qchisq(pp, df = 3),
       y = quantile(plot_dat[, "wald_chi_sq"], pp),
       pch = 3,
       col = cbbPalette[3])
abline(a = 0, b = 1, col = "red")
if(PDF) dev.off()
