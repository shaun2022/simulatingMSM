# This program reproduces Tables A5 and A6 when n=100000.

n <- 100000 # Total number of individuals to generate.
# If you want to reduce the running time, use n=5000 instead.
overall.seed.vals <- c(10, 20)
# These are the seeds that will be used to generate multiple seeds for the copies.
copula.cor.vals <- c(0.9, 0.5) # The (negative of the) correlation parameter of Gaussian copula.

# Five values for the number of matches (plus the original) at each pass.
ncopy.passes.vals <- matrix(0, 5, 2)
ncopy.passes.vals[, 1] <- c(10, 100, 500, 1000, 5000)
ncopy.passes.vals[, 2] <- c(2000, 2000, 10000, 20000, 100000)
# Five values for minimum allowed number of surviving matches for 1st and 2nd passes.
min.copies.vals <- matrix(0, 5, 2)
min.copies.vals[, 1] <- c(0, 10, 50, 100, 500)
num.ncopy <- nrow(ncopy.passes.vals)

# Arrays to store the failure times and the failure indicators for the n originals
# in each of the scenarios.  The indices are: 1) individual; 2) marginal risk;
# 3) number of copies; 4) random seed for the matches.
tim.sens <- array(NA, c(n, 3, num.ncopy, 2))
fail.sens <- array(NA, c(n, 3, num.ncopy, 2))

for (which.copula.cor in 1:2) {
  # All of the following code will be run for each of the two values of the Gaussian copula.
  copula.cor <- copula.cor.vals[which.copula.cor]
  for (MARGINAL.RISK in 1:3)
    for (which.ncopy in 1:num.ncopy)
      for (which.seed in 1:2) {
        # The following code will be run for each marginal risk, number of copies and (2) seeds.
        ncopy.passes <- ncopy.passes.vals[which.ncopy,]
        min.copies <- min.copies.vals[which.ncopy,]
        OVERALL.SEED <- overall.seed.vals[which.seed]
        # This seed is for generating data for the matches.
        
        print(paste0("MARGINAL.RISK=", MARGINAL.RISK, ", which.ncopy=", which.ncopy, ", seed=",
                     OVERALL.SEED, ", copula.cor=", copula.cor))
        
        try( {
          # Using try() because may fail when only 9 copies are used and marginal risk is high.
          # Generate data for the n originals.
          source("simulateMSM_sensitivity.R")

          # Put the data on failure times into the arrays for storing the failure times
          # for all the scenarios.
          tim.sens[, MARGINAL.RISK, which.ncopy, which.seed] <- tim.all
          fail.sens[, MARGINAL.RISK, which.ncopy, which.seed] <- fail.all
        } , silent=T )
      }

  # Set failure time equal to ntimes+1 for individuals who do not fail
  tim.sens[!fail.sens] <- ntimes+1
  
  # Set up array to store the numbers of individuals whose failure times differ by
  # 0, 1, 2, ... when different number of matches are used.
  tim.diff <- array(0, c(3, num.ncopy-1, ntimes+1))

  for (MARGINAL.RISK in 1:3)
    for (which.ncopy in 1:(num.ncopy-1))
      for (which.seed in 1:2) {
        absdiff <- abs( tim.sens[, MARGINAL.RISK, which.ncopy, which.seed] -
                          tim.sens[, MARGINAL.RISK, num.ncopy, which.seed] )
        # absdiff[i, MARGINAL.RISK, which.ncopy, which.seed] is absolute difference between sampled
        # individual i's failure time when <4999 matches are used and when 4999 copies are used.
        tim.diff[MARGINAL.RISK, which.ncopy,] <- tim.diff[MARGINAL.RISK, which.ncopy,] +
                                                 tabulate(absdiff+1, nbin=ntimes+1)
        # tim.diff[MARGINAL.RISK, which.ncopy, j] is number of sampled individuals whose failure
        # times differ by j-1.  In particular, tim.diff[MARGINAL.RISK, which.ncopy, 1] is the number
        # whose failure times do not differ.
      }
  
  # Set up array to store the numbers of individuals whose failure times differ by 0, 1, 2, ...
  # when different random seeds are used.
  tim.seed.diff <- array(0, c(3, num.ncopy, ntimes+1))
  
  for (MARGINAL.RISK in 1:3)
    for (which.ncopy in 1:num.ncopy) {
      absdiff <- abs( tim.sens[, MARGINAL.RISK, which.ncopy, 1] -
                        tim.sens[, MARGINAL.RISK, which.ncopy, 2] )
      tim.seed.diff[MARGINAL.RISK, which.ncopy,] <- tabulate(absdiff+1, nbin=ntimes+1)
    }
  
  # Output results to a file.
  sink(paste0("results/Table_A", which.copula.cor+4, "_samplesize", n, ".txt"))
  
  mytable <- cbind( rbind( t(tim.diff[,,1]), c(n*2, n*2, n*2) ) / 2, t(tim.seed.diff[,,1]) ) / n
  rownames(mytable) <- ncopy.passes.vals[, 1]
  colnames(mytable) <- c("vs.5000.low", "vs.5000.med", "vs.5000.high", "vs.same.low",
                         "vs.same.med", "vs.same.high")
  print( round(mytable*100, 1) )

  sink()
} # End of for (which.copula.cor in 1:2) loop.
