# This program is called by the file `simulateMSM_runall.R'.

library(sandwich) # For calculating sandwich SEs
expit <- function(x) return( exp(x) / (1+exp(x)) ) # The inverse logistic function.

OVERALL.SEED <- 20 # This seed will be used to generate multiple seeds.

ncopy.passes <- c(5000, 100000) # Number of matches plus 1 for 1st and 2nd passes.
min.copies <- c(500, 0) # Minimum allowed number of surviving matches for 1st and 2nd passes.
nmarg <- 2 # Number of baseline variables in marginal model.
nbase <- 2 # Number of baseline variables not in marginal model.
nvary <- 2 # Number of time-dependent confounders.
ntimes <- 10 # Number of visits.
copula.cor <- 0.9 # Negative of the correlation parameter of Gaussian copula.

# If MARGINAL.RISK is, respectively, 1/2/3, then 10%/50%/90% of individuals fail before time 11.
# haz.base is baseline (discrete-time) hazard at each time point.
# It is -4.1 for 10% failing, -2.5 for 50% failing, and -1.2 for 90% failing.
if (MARGINAL.RISK==1)
  haz.base <- rep(-4.1, ntimes)
if (MARGINAL.RISK==2)
  haz.base <- rep(-2.5, ntimes)
if (MARGINAL.RISK==3)
  haz.base <- rep(-1.2, ntimes)

set.seed(OVERALL.SEED)
# Generate 1+n/1000 random number seeds.
seeds <- floor( runif(1 + ceiling(n/1000)) * 1e8 )

# Set up matrices, arrays and vectors to store all the simulated data.
# n is the number of (original) individuals to generate data for.
m.all <- matrix(0, n, nmarg) # Baseline covariates in MSM (called "X" in the article).
b.all <- matrix(0, n, nbase) # Baseline variables not in MSM.
x.all <- array(0, c(n, nvary, ntimes)) # Time-dependent confounders (called "L" in the article).
a.all <- matrix(0, n, ntimes) # Time-dependent treatment.
tim.all <- rep(ntimes, n) # Failure times (will be ntimes if does not fail).
fail.all <- rep(0, n) # Indicator of failure.
nsurv.copy.all <- rep(0, n) # Number of surviving matches. 
last.enough.all <- rep(0, n) # Last time point where there are at least min.copies matches.

set.seed(seeds[1])

# Sample baseline covariates in MSM for the n originals
m.all[,1] <- rnorm(n)
m.all[,2] <- runif(n) < 0.5

for (whichpass in 1:2) {
  ncopy <- ncopy.passes[whichpass] # Number of matches plus the one original.
  # If whichpass=1, nsurv.copies.all equals ncopy.passes[1] for everyone.
  # If whichpass=2, nsurv.copies.all equals ncopy.passes[2] for everyone who has not yet failed.
  nsurv.copy.all[last.enough.all < ntimes & !fail.all] <- ncopy.passes[whichpass]
  
  # The following matrices will be needed later.
  b <- matrix(0, ncopy, nbase)
  x <- matrix(0, ncopy, nvary)
  
  # Parameters of data-generating model for B given X.
  b.coeff <- matrix(0, nbase, 1+nmarg) 
  b.coeff[1,] <- c(-0.2, 0, 0.4)
  b.coeff[2,] <- c(0, 0.2, 0)
  
  for (i in (1:n)[last.enough.all < ntimes & !fail.all]) {
      # Reset the seed every 1000 individuals to make debugging easier.
    if ( whichpass==1 & (i+999)%%1000==0 ) {
      print(paste0("Generating data for individual ", i))
      set.seed(seeds[(i+999)%/%1000 + 1])
    }
    
    if (whichpass==2)
      print(paste0("Refreshing matches for individual ", i))
    
    fail <- rep(F, ncopy)  
    # iden will track which matches `donated' themselves to failed matches (called "I" in article).
    iden <- 1:ncopy
    
    # Sample baseline variables not in MSM for original and matches
    for (j in 1:nbase) {
      b[,j] <- as.vector( b.coeff[j,1] + m.all[i,] %*% b.coeff[j,-1] ) + rnorm(ncopy)
    }
    if (whichpass==1) {
      # Store the generated value of b for the original
      b.all[i,] <- b[1,]
    } else {
      # Keep the previously generated value of b for the original if it has already been generated.
      b[1,] <- b.all[i,]
    }
    
    timep <- 1 # Counter used to denote visit.
    
    while (timep <= ntimes & !fail[1] & nsurv.copy.all[i] >= min.copies[whichpass]) {
      if (timep==1) {
        x[,1] <- 0.2 * m.all[i, 1] + rnorm(ncopy)
        x[,2] <- runif(ncopy) < expit( -0.2 + 0.4 * m.all[i, 2] )
      } else {
        x.prev <- x
        x[,1] <- 0.3 + 0.4 * b[,2] + 0.7 * x.prev[,1] - 0.6 * a.all[i, timep-1] + rnorm(ncopy)
        x[,2] <- runif(ncopy) < expit(-0.2 + 0.4 * b[, 2] + 1 * x.prev[,2] - 0.6 * a.all[i, timep-1])
      }
      
      # Unless timep is after last.enough, keep value of x for the original.
      # If timep is after last.enough, sample treatment for the original.
      if (timep > last.enough.all[i]) {
        x.all[i,, timep] <- x[1,]
        logitp <- -1 + m.all[i,] %*% c(0.2, 0.3) + b.all[i,] %*% c(0.2, 0) + x[1,] %*% c(0.6, 0.6)
        if (timep>1)
          logitp <- logitp + a.all[i, timep-1]
        a.all[i, timep] <- runif(1) < expit(logitp)
      } else {
        x[1,] <- x.all[i,, timep]
      }
      
      # Calculate the negative of H, the `risk score'.
      h <- as.vector( b %*% c(-0.3, -0.5) + x %*% c(-1, -1) )
      
      # Use the copula to get from H to u.y (which determines who fails).
      u.h <- (rank(h) - runif(ncopy)) / ncopy
      z.h <- qnorm(u.h)
      z.y <- copula.cor * z.h + rnorm(ncopy, sd=sqrt(1 - copula.cor^2))
      u.y <- pnorm(z.y)

      # Calculate probability of failure given m.all and a.all.
      prob.fail <- as.vector( expit( haz.base[timep] + m.all[i,] %*% c(0.5, 0.5)
                                    + a.all[i, timep] * (-1) ) )
      
      # Establish which originals and matches fail at this time point.
      ind <- u.y < prob.fail
      # But don't kill the original unless timep is after last.enough.
      if (timep <= last.enough.all[i])
        ind[1] <- F
      fail[ind==T] <- T
        
      if (fail[1]) {
        # If the original has failed, store this information.
        fail.all[i] <- T
        tim.all[i] <- timep
      }

      if (!fail[1] & timep < ntimes) {
        # Replace the set of matches who fail with a random sample of matches who did not fail.
        fail.m1 <- fail[-1] # Failure indicators of only the matches.
        if (timep==1) {
          donors <- rep((2:ncopy)[!fail.m1], length.out=sum(fail.m1))
        } else {
          donors <- sample(x=(2:ncopy)[!fail.m1], size=sum(fail.m1), replace=T)
        }
        b[fail,] <- b[donors,]
        x[fail,] <- x[donors,]
        iden[fail] <- iden[donors]
        fail <- rep(F, ncopy)
      } # End of if (!fail[1] & timep < ntimes) condition.

      timep <- timep + 1
      nsurv.copy.all[i] <- length(unique(iden[-1]))
    } # End of while (timep <= ntimes & !fail[1] & nsurv.copy.all[i] >= min.copies[whichpass]) loop.
    
    last.enough.all[i] <- timep - 1
  } # End of for (i in (1:n)[last.enough.all < ntimes & !fail.all]) loop.
} # End of for (whichpass in 1:2) loop

# The data have now been simulated.  Now analyse the simulated data.


# Calculate stabilised weights for the MSM.
wei <- matrix(1, length(tim.all), ntimes)
p.denom <- rep(0, ntimes)
p.numer <- rep(0, ntimes)

for (timep in 1:ntimes) {
  use <- tim.all >= timep
  if (timep==1) {
    # Fit models for the denominator and numerator of the stabilised weights.
    fit.denom <- glm(a.all[, timep] ~ m.all + b.all[,1] + x.all[,, timep], family=binomial)
    fit.numer <- glm(a.all[, timep] ~ m.all, family=binomial)
  }
  if (timep>1) {
    fit.denom <- glm(a.all[, timep] ~ m.all + b.all[,1] + x.all[,, timep] + a.all[, timep-1],
                     family=binomial, subset=use==T)
    fit.numer <- glm(a.all[, timep] ~ m.all + a.all[, timep-1], family=binomial, subset=use==T)
  }
  p.denom[use==T] <- a.all[use==T, timep] * fit.denom$fitted.values +
    (1-a.all[use==T, timep]) * (1-fit.denom$fitted.values)
  p.numer[use==T] <- a.all[use==T, timep] * fit.numer$fitted.values +
    (1-a.all[use==T, timep]) * (1-fit.numer$fitted.values)
  if (timep==1) {
    wei[, 1] <- p.numer / p.denom
  } else {
    wei[, timep] <- wei[, timep-1] * p.numer / p.denom
  }
  wei[use==F, timep] <- 0
}

# Truncate the weights to deal with the possibility that a small number of them might be enormous.
# The truncation is very `generous': weights of up to 1000 are allowed.
wei[wei>1000] <- 1000

# Set up vectors to store results of fitting MSM using with and without weights.
est.nowei <- rep(0, ntimes+nmarg*2+2)
est.wei <- rep(0, ntimes+nmarg*2+2)
se.sand.wei <- rep(0, ntimes+nmarg*2+2)
se.model.nowei <- rep(0, ntimes+nmarg*2+2)
se.model.wei <- rep(0, ntimes+nmarg*2+2)

# Fit the MSM using no weights and using inverse probability of treatment weights.
# Need to put data in `long' format for this.
# The following vectors will contain the y, m1, m2, a, t and weight values for each individual and
# each visit (until individual fails).
y.vec <- NULL
m1.vec <- NULL
m2.vec <- NULL
a.vec <- NULL
t.vec <- NULL
wei.vec <- NULL

for (j in 1:ntimes) {
  stillatrisk <- tim.all >= j
  y.vec <- c(y.vec, (tim.all==j & fail.all)[stillatrisk==T] * 1)
  m1.vec <- c(m1.vec, m.all[stillatrisk==T, 1])
  m2.vec <- c(m2.vec, m.all[stillatrisk==T, 2])
  a.vec <- c(a.vec, a.all[stillatrisk==T, j])
  t.vec <- c(t.vec, rep(j-1, sum(stillatrisk==T)))
  wei.vec <- c(wei.vec, wei[stillatrisk==T, j])
}

# Fit model without weights.
fit.nowei <- glm(y.vec ~ factor(t.vec) + m1.vec + m2.vec + a.vec + m1.vec:t.vec + m2.vec:t.vec +
                   a.vec:t.vec - 1, family=binomial)
# Store parameter estimates and model-based estimates of SEs.
est.nowei <- fit.nowei$coefficients
se.model.nowei <- sqrt(diag(summary(fit.nowei)$cov.scaled))

# Fit model with weights.
fit.wei <- glm(y.vec ~ factor(t.vec) + m1.vec + m2.vec + a.vec + m1.vec:t.vec + m2.vec:t.vec +
                 a.vec:t.vec - 1, family=binomial, weights=wei.vec)
# Store parameter estimates and sandwich estimates of SEs.
est.wei <- fit.wei$coefficients
sand.wei <- sandwich(fit.wei, bread=bread(fit.wei), meat=meat(fit.wei))
se.sand.wei <- sqrt(diag(sand.wei))

# Output results to file.
whichtable <- switch(MARGINAL.RISK, "A2", "A1", "A3")
sink(paste0("results/Table_", whichtable, "_samplesize", n, ".txt"))
results.mat <- cbind(c(haz.base, c(0.5, 0.5, -1, 0, 0, 0)), est.nowei, se.model.nowei, est.wei,
                     se.sand.wei)
colnames(results.mat) <- c("true", "unwei.est", "unwei.se", "wei.est", "wei.se")
rownames(results.mat) <- c(paste0("beta", 0:(ntimes-1), "0"), paste0("beta", 1:6))
print(round(results.mat, 3))
sink()
