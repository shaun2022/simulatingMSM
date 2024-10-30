##############################################################
#Master file for generating Tables 1, 2, 3 and Table A7
##############################################################

include.bootstrap<-1 #change this to 0 to exclude the bootstrap analyses.

#--------------------------
#Packages and functions
#--------------------------

library(survival) # For the coxph function when using Cox MSM
library(sandwich) # For calculating sandwich SEs
library(dplyr) #For general data manipulation

expit <- function(x) return( exp(x) / (1+exp(x)) )
loglog <- function(x) return( (1-exp(-exp(x))) ) # analogous to expit

#--------------------------
#Simulation settings
#--------------------------

#number of simulated data sets
nsim<-1000

#number of bootstrap samples to be used in bootstrap approach
nboot<-1000

#true coefficients (log HRs) in the Cox MSM (used in obtaining coverage estimates)
true.coefs<-c(0.5,0.5,-1,-0.4,0)

# n is number of individuals in each simulated data set.

# copula.cor is correlation parameter of Gaussian copula (note this is -rho in
# the notation of the paper).  This is 0.5 for rho.low and 0.9 for rho.high.

# propscore.mult is multiplier for coefficients in the propensity score model.
# It is 0.5 for delta.low and 1 for delta.high.

for (n in c(250, 500, 1000))
  for (copula.cor in c(0.5, 0.9))
    for (propscore.mult in c(0.5, 1)) {
      print(paste0("n=", n, ", copula.cor=", copula.cor, ", propscore.mult=", propscore.mult))
      
      # Run simulation
      source("simulation_analysis.R")

      # Generate intermediate results for Tables 1, 2, 3 and Table A7
      source("simulation_saveresults.R")
    }

# Create the final tables
source("simulation_createtables.R")
