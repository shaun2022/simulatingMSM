# This R code generates Figure 2

gauscopula <- function(g, uh, rho)
  return( pnorm( ( qnorm(g) - rho * qnorm(uh) ) / sqrt(1 - rho^2) ) )

uhval <- seq(0.00001, 0.99999, 0.00001)

marg.prob.vals <- c(0.05, 0.1, 0.5, 0.9)
rho <- - c(0, 0.1, 0.5, 0.9, 0.99)
nrho <- length(rho)
prob.rho <- matrix(0, nrow=length(uhval), ncol=length(rho))

pdf("copulaplots.pdf")

par(mfrow=c(2,2))
par(mgp=c(1.7, 0.5, 0), mar=c(3, 3, 3, 0.5))

for (marg.prob in marg.prob.vals) {
  for (j in 1:nrho)
    prob.rho[,j] <- gauscopula(marg.prob, uhval, rho[j])
  plot(uhval, prob.rho[,1], type="l", lwd=2, xlim=c(0,1), ylim=c(0,1),
    xlab="risk quantile", ylab="P(failure | risk quantile)", cex.lab=1.4)
  for (j in 2:nrho)
    lines(uhval, prob.rho[,j], col=j, lwd=2)
  title(paste0("marginal P(failure)=", marg.prob))
  if (marg.prob==0.05)
    legend("topleft", lty=rep(1, nrho), lwd=2, col=c(1:nrho), legend=rho)
}

dev.off()
