betanames <- paste0("beta", 1:5)

# Set up tables into which to insert the intermediate results

# Table 1
table1 <- matrix(0, 20, 6)
rownames(table1) <- c( paste0(betanames, ".lowlow"), paste0(betanames, ".highlow"), paste0(betanames, ".lowhigh"), paste0(betanames, ".highhigh"))
colnames(table1) <- c("naive.1000", "IPTW.1000", "naive.500", "IPTW.500", "naive.250", "IPTW.250")

# Table 2
table2 <- matrix(0, 20, 9)
rownames(table2) <- rownames(table1)
colnames(table2) <- c("naive.1000", "sand.1000", "boot.1000", "naive.500", "sand.500", "boot.500", "naive.250", "sand.250", "boot.250")

# Table 3
table3 <- matrix(0, 4, 9)
rownames(table3) <- c("delta.low, rho.low", "delta.high, rho.low", "delta.low, rho.high", "delta.high, rho.high")
colnames(table3) <- colnames(table2)

# Table A7
tableA7 <- table2

# Values of n, rho and delta for the 12 scenarios
n.vals <- c(1000, 500, 250)
copula.cor.vals <- c(0.5, 0.5, 0.9, 0.9)
propscore.mult.vals <- c(0.5, 1, 0.5, 1)

# Cycle through the 12 scenarios, inserting the intermediate results into the
# appropriate places in the tables
for (i in 1:3)
  for (j in 1:4) {
    n <- n.vals[i]
    copula.cor <- copula.cor.vals[j]
    propscore.mult <- propscore.mult.vals[j]

    # Appropriate rows and columns of the tables into which to insert
    # intermediate results for this scenario
    table12A7.rows <- (j-1)*5 + 1:5
    table3.row <- j
    table1.cols <- (i-1)*2 + 1:2
    table23A7.cols <- (i-1)*3 + 1:3

    filename <- paste0("nsim", nsim, "_n", n, "_nboot", nboot, "_copula.cor", copula.cor, "_propscore.mult", propscore.mult, ".out")

    # Table 1
    intermediate.result <- as.matrix(read.table(paste0("intermediate_results/table1_", filename)))
    table1[table12A7.rows, table1.cols] <- intermediate.result[, 1:2]

    if (include.bootstrap) {
      # Can only create Tables 2, 3 and A7 if bootstrap has been used
      
      # Table 2
      intermediate.result <- as.matrix(read.table(paste0("intermediate_results/table2_", filename)))
      table2[table12A7.rows, table23A7.cols] <- intermediate.result[, 1:3]
      
      # Table 3
      intermediate.result <- as.vector(read.table(paste0("intermediate_results/table3_", filename))$x)
      table3[table3.row, table23A7.cols] <- intermediate.result[1:3]
      
      # Table A7
      intermediate.result <- as.matrix(read.table(paste0("intermediate_results/tableA7_", filename)))
      tableA7[table12A7.rows, table23A7.cols] <- intermediate.result[, 1:3]
    }
  }

# Create files containing the four tables.
options("width"=120)

sink("results/Table1.txt")
print(round(table1, 3))
sink()

if (include.bootstrap) {
  sink("results/Table2.txt")
  print(round(table2, 3))
  sink()
  
  sink("results/Table3.txt")
  print(round(table3, 3))
  sink()
  
  sink("results/TableA7.txt")
  print(round(tableA7, 3))
  sink()
}
