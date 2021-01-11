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
# par(bty="n")
#par(mfrow = c(1, 3))

# Likelihood
plot(lam_vals, dat1[, "ll"], type = "l", lwd = 2,
     ylab = "log-likelihood", xlab = expression(lambda) ,
     ylim = c(min(dat2[, "ll"]), max(dat1[, "ll"])),
     col = cbbPalette[1])
lines(lam_vals, dat2[, "ll"], type = "l", lwd = 2, lty = 1,
      col = cbbPalette[1])

if(PDF) dev.off()

### Simulations ####



