# This program is called by the file `simulateMSM_sensitivity_runall.R'.

expit <- function(x) return( exp(x) / (1+exp(x)) ) # The inverse logistic function.

ORIGINAL.SEED <- 20
# This seed is for generating the data for the originals (rather than for the matches).

nmarg <- 2 # Number of baseline variables in marginal model.
nbase <- 2 # Number of baseline variables not in marginal model.
nvary <- 2 # Number of time-varying confounders.
ntimes <- 10 # Number of visits.

# haz.base is baseline (discrete-time) hazard at each time point.
# It is -4.1 for 10% failing, -2.5 for 50% failing, and -1.2 for 90% failing
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
m.all <- matrix(0, n, nmarg) # Baseline covariates in MSM (called "X" in the article).
b.all <- matrix(0, n, nbase) # Baseline variables not in MSM.
x.all <- array(0, c(n, nvary, ntimes)) # Time-dependent confounders (called "L" in the article).
a.all <- matrix(0, n, ntimes) # Time-dependent treatment.
tim.all <- rep(ntimes, n) # Failure times (will be ntimes if does not fail).
fail.all <- rep(0, n) # Indicator of failure.
nsurv.copy.all <- rep(0, n) # Number of surviving matches. 
last.enough.all <- rep(0, n) # Last time point where at least min.copies matches.

set.seed(ORIGINAL.SEED)
  
# Sample baseline covariates in MSM for the n originals.
m.all[,1] <- rnorm(n)
m.all[,2] <- runif(n) < 0.5

# Parameters of data-generating model for B given X.
b.coeff <- matrix(0, nbase, 1+nmarg) 
b.coeff[1,] <- c(-0.2, 0, 0.4)
b.coeff[2,] <- c(0, 0.2, 0)

# Sample baseline variables not in MSM for the n originals.
for (j in 1:nbase)
  b.all[,j] <- as.vector( b.coeff[j,1] + m.all %*% b.coeff[j,-1] ) + rnorm(n)

# Sample time-dependent variables (including treatment) for the n originals.
for (timep in 1:ntimes) {
  if (timep==1) {
    x.all[, 1, timep] <- 0.2 * m.all[, 1] + rnorm(n)
    x.all[, 2, timep] <- runif(n) < expit( -0.2 + 0.4 * m.all[, 2] )
  } else {
    x.all[, 1, timep] <- 0.3 + 0.4 * b.all[,2] + 0.7 * x.all[, 1, timep-1] -
                         0.6 * a.all[, timep-1] + rnorm(n)
    x.all[, 2, timep] <- runif(n) < expit( -0.2 + 0.4 * b.all[, 2] +
                                             1 * x.all[, 2, timep-1] - 0.6 * a.all[, timep-1] )
  }
  
  logitp <- -1 + m.all %*% c(0.2, 0.3) + b.all %*% c(0.2, 0) + x.all[,, timep] %*% c(0.6, 0.6)
  if (timep>1)
    logitp <- logitp + a.all[, timep-1]
  a.all[, timep] <- runif(n) < expit(logitp)
} # end of for (timep in 1:ntimes) loop.

# Generate random numbers that will be used later when calculating which of the n originals
# fails and when.
u.h.rand <- matrix( runif(n*ntimes), n, ntimes )
z.y.rand <- matrix( rnorm(n*ntimes, sd=sqrt(1 - copula.cor^2)), n, ntimes )

for (whichpass in 1:2) {
  ncopy <- ncopy.passes[whichpass] # Number of matches plus one.
  # If whichpass=1, nsurv.copies.all equals ncopy.passes[1] for everyone.
  # If whichpass=2, nsurv.copies.all equals ncopy.passes[2] for everyone who has not yet failed.
  nsurv.copy.all[last.enough.all < ntimes & !fail.all] <- ncopy.passes[whichpass]
  
  # The following matrices will be needed later.
  b <- matrix(0, ncopy, nbase)
  x <- matrix(0, ncopy, nvary)
  
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
    
    # Sample baseline variables not in MSM for original and matches.
    for (j in 1:nbase)
      b[,j] <- as.vector( b.coeff[j,1] + m.all[i,] %*% b.coeff[j,-1] ) + rnorm(ncopy)
    # Keep the value of b for the original.
    b[1,] <- b.all[i,]
    
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

      # Keep previously generated value of x for the original.
      x[1,] <- x.all[i,, timep]
      
      # Calculate the negative of H, the `risk score'.
      h <- as.vector( b %*% c(-0.3, -0.5) + x %*% c(-1, -1) )
      
      # Use the copula to get from H to u.y (which determines who fails).
      rand <- runif(ncopy)
      rand[1] <- u.h.rand[i]
      u.h <- (rank(h) - rand) / ncopy
      z.h <- qnorm(u.h)
      rand <- rnorm(ncopy, sd=sqrt(1 - copula.cor^2))
      # For the original, use the previously generated random number.
      rand[1] <- z.y.rand[i, timep]
      z.y <- copula.cor * z.h + rand
      u.y <- pnorm(z.y)
      
      # Calculate probability of failure given m.all and a.all.
      prob.fail <- as.vector( expit( haz.base[timep] + m.all[i,] %*% c(0.5, 0.5) +
                                       a.all[i, timep] * (-1) ) )
        
      # Establish which originals and matches fail at this time point.
      ind <- u.y < prob.fail
      # But don't kill the original unless timep is after last.enough.
      if (timep <= last.enough.all[i])
        ind[1] <- F
      fail[ind==T] <- T
        
      if (fail[1]) {
        # The original has failed.
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
      } # End of if (!fail[1] & timep < ntimes) condition
        
      timep <- timep + 1
      nsurv.copy.all[i] <- length(unique(iden[-1]))
    } # End of while (timep <= ntimes & !fail[1] & nsurv.copy.all[i] >= min.copies[whichpass]) loop.
      
    last.enough.all[i] <- timep - 1
  } # End of for (i in (1:n)[last.enough.all < ntimes & !fail.all]) loop.
} # End of for (whichpass in 1:2) loop.
