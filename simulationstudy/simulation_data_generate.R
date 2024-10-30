##############################################################
#This file simulates a longitudinal data set as described in Section 6.1.
#This file is called by `simulation_analysis.R'.
##############################################################

#--------------
# Data generation settings used in the paper, but which could be modified to generate other scenarios
#--------------

ncopy.passes <- c(5000, 100000) # number of matches plus 1 for 1st and 2nd passes
min.copies <- c(500, 0) # minimum allowed number of surviving matches for 1st and 2nd passes
nmarg <- 2 # number of baseline variables in marginal model
nbase <- 2 # number of baseline variables not in marginal model
nvary <- 2 # number of time-varying covariates
ntimes <- 10 # number of visits (K=9)

haz.base <- rep(-3.3, 10) # baseline hazard for event times
lambda.cens <-exp(-3.6) # hazard for censoring times (should give 50% censored before time 10)

#--------------
# Set up matrices, arrays and vectors to store all the simulated data
#--------------

m.all <- matrix(0, n, nmarg) # baseline covariates in MSM (called "X" in the article)
b.all <- matrix(0, n, nbase) # baseline variables not in MSM
x.all <- array(0, c(n, nvary, ntimes)) # time-dependent confounders (called "L" in the article)
a.all <- matrix(0, n, ntimes) # time-dependent treatment
tim.all <- rep(ntimes, n) # failure times (will be ntimes if does not fail)
fail.all <- rep(0, n) # indicator of failure
nsurv.copy.all <- rep(0, n) # number of surviving matches 
last.enough.all <- rep(0, n) # last time point where at least min.copies matches

#--------------
#Start of simulation algorithm
#--------------
  
# sample baseline covariates in MSM for the n originals
m.all[,1] <- rnorm(n)
m.all[,2] <- runif(n) < 0.5

for (whichpass in 1:2) {
  ncopy <- ncopy.passes[whichpass] # number of matches plus the one original
  # If whichpass=1, nsurv.copies.all equals ncopy.passes[1] for everyone.
  # If whichpass=2, nsurv.copies.all equals ncopy.passes[2] for everyone who has not yet failed.
  nsurv.copy.all[last.enough.all < ntimes & !fail.all] <- ncopy.passes[whichpass]
  
  # The following matrices will be needed later
  b <- matrix(0, ncopy, nbase)
  x <- matrix(0, ncopy, nvary)
  
  # parameters of data-generating model for B given X
  b.coeff <- matrix(0, nbase, 1+nmarg) 
  b.coeff[1,] <- c(-0.2, 0, 0.4)
  b.coeff[2,] <- c(0, 0.2, 0)
  
  for (i in (1:n)[last.enough.all < ntimes & !fail.all]) {
    # Refreshing matches for individual i    
    fail <- rep(F, ncopy)
    # iden will track which matches `donated' themselves to failed matches (called "I" in article).
    iden <- 1:ncopy
    
    # Sample baseline variables not in MSM for original and matches
    for (j in 1:nbase) {
      b[,j] <- as.vector( b.coeff[j,1] + m.all[i,] %*% b.coeff[j,-1] ) + rnorm(ncopy)
    }
    # Keep the value of b for the original if it has already been generated
    if (whichpass==1) {
      # Store the generated value of b for the original
      b.all[i,] <- b[1,]
    } else {
      # Keep the previously generated value of b for the original if it has already been generated.
      b[1,] <- b.all[i,]
    }
    
    timep <- 1 # counter used to denote visit
    
    while (timep <= ntimes & !fail[1] & nsurv.copy.all[i] >= min.copies[whichpass]) {
      if (timep==1) {
        x[,1] <- 0.2 * m.all[i, 1] + rnorm(ncopy)
        x[,2] <- runif(ncopy) < expit( -0.2 + 0.4 * m.all[i, 2] )
      } else {
        x.prev <- x
        x[,1] <- 0.3 + 0.4 * b[,2] + 0.7 * x.prev[,1] - 0.6 * a.all[i, timep-1] + rnorm(ncopy)
        x[,2] <- runif(ncopy) < expit( -0.2 + 0.4 * b[, 2] + 1 * x.prev[,2] - 0.6 * a.all[i, timep-1] )
      }
      
      # Unless timep is after last.enough, keep value of x for the original.
      # If timep is after last.enough, sample treatment for the original.
      if (timep > last.enough.all[i]) {
        x.all[i,, timep] <- x[1,]
        logitp <- -1 + propscore.mult*(m.all[i,] %*% c(0.2, 0.3) + b.all[i,] %*% c(0.2, 0) +
                                         x[1,] %*% c(0.6, 0.6))
        if (timep>1)
          logitp <- logitp + a.all[i, timep-1]
        a.all[i, timep] <- runif(1) < expit(logitp)
      } else {
        x[1,] <- x.all[i,, timep]
      } # end of if (timep > last.enough.all[i]) condition and its else statement
      
      # Calculate the negative of H, the `risk score'.
      h <- as.vector( b %*% c(-0.3, -0.5) + x %*% c(-1, -1) )
      
      # Use the copula to get from H to u.y (which determines who fails)
      u.h <- (rank(h) - runif(ncopy)) / ncopy
      z.h <- qnorm(u.h)
      z.y <- copula.cor * z.h + rnorm(ncopy, sd=sqrt(1 - copula.cor^2))
      u.y <- pnorm(z.y)
      
      # Calculate probability of failure given m.all and a.all
      prob.fail <- as.vector( loglog(haz.base[timep] + m.all[i,] %*% c(0.5, 0.5) +
                                     a.all[i, timep] * (-1) + m.all[i,1]*a.all[i, timep] * (-0.4)) )
      
      # Establish which originals and matches fail at this time point
      ind <- u.y < prob.fail
      # But don't kill the original unless timep is after last.enough
      if (timep <= last.enough.all[i])
        ind[1] <- F
      fail[ind==T] <- T
      
      if (fail[1]) {
        # If the original has failed, store this information.
        fail.all[i] <- T
        tim.all[i] <- timep - 1 + log(1-u.y)[1] / log(1-prob.fail)
      }
      
      if (!fail[1] & timep < ntimes) {
        # Replace the set of matches who fail with a random sample of matches who did not fail.
        fail.m1 <- fail[-1] # failure indicators of only the matches
        if (timep==1) {
          donors <- rep((2:ncopy)[!fail.m1], length.out=sum(fail.m1))
        } else {
          donors <- sample(x=(2:ncopy)[!fail.m1], size=sum(fail.m1), replace=T)
        }
        b[fail,] <- b[donors,]
        x[fail,] <- x[donors,]
        iden[fail] <- iden[donors]
        fail <- rep(F, ncopy)
      } # end of if (!fail[1] & timep < ntimes) condition
      
      timep <- timep + 1
      nsurv.copy.all[i] <- length(unique(iden[-1]))
    } # end of while (timep <= ntimes & !fail[1] & nsurv.copy.all[i] >= min.copies[whichpass]) loop
    
    last.enough.all[i] <- timep - 1
  } # end of for (i in (1:n)[last.enough.all < ntimes & !fail.all]) loop
} # end of for (whichpass in 1:2) loop
  

#--------------
#GENERATE RANDOM CENSORING USING EXPONENTIAL HAZARD
#--------------

cens.time<-(-log(runif(n,0,1))/lambda.cens)^(1/1)

obs.time<-ifelse(tim.all<=cens.time,tim.all,cens.time)

obs.fail<-ifelse(tim.all<=cens.time,fail.all,0)

#this added an extra indicator for administrative censoring for the purposes of getting percentages with each end type
obs.fail<-ifelse(obs.time==10 & obs.fail==0,2,obs.fail)
prop.table(table(obs.fail))
  
#--------------
#CREATE A DATA FRAME 
#--------------

#---
#wide format
dat<-data.frame(id=1:n,time=tim.all,fail=fail.all,
                m.all,b.all,a.all,x.all[,1,],x.all[,2,])
names(dat)=c("id","time","fail","M1","M2","B1","B2",paste0("A.",1:10),paste0("X1.",1:10),paste0("X2.",1:10))

#---
#long format

dat.long<-reshape(dat,
                  varying=c(paste0("A.",1:10),paste0("X1.",1:10),paste0("X2.",1:10)),
                  timevar="visit.time",idvar="id",direction="long")
dat.long$visit.time<-dat.long$visit.time-1

dat.long<-dat.long[order(dat.long$id,dat.long$visit.time),]

#delete rows after the person's event/censoring time
dat.long<-dat.long[dat.long$time>=dat.long$visit.time,]

#generate lagged treatment
dat.long<-dat.long%>%group_by(id)%>%mutate(A_lag=lag(A,default=0))

#generate row number and total number of rows
dat.long<-dat.long%>%group_by(id)%>%mutate(rownum=row_number())%>%mutate(maxrow=max(rownum))

#generate start and stop times in each row, and corresppnding time-updated fail indicator
dat.long$tstart<-dat.long$visit.time

dat.long$tstop<-ifelse(dat.long$rownum!=dat.long$maxrow,dat.long$visit.time+1,dat.long$time)

dat.long$fail.updated<-ifelse(dat.long$rownum!=dat.long$maxrow,0,dat.long$fail)
