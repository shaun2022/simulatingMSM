# This program reproduces Tables A1, A2 and A3.

n <- 1e6 # Sample size
# If you want to reduce the running time, use n=50000 instead.

for (MARGINAL.RISK in 1:3) {
  # When MARGINAL.RISK=1/2/3, 10%/50%/90% fail before time 11
  print(paste0("MARGINAL.RISK"=MARGINAL.RISK))
  source("simulateMSM.R")
}
