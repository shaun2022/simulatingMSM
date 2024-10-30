##############################################################
#This file creates Tables 1,2,3 and A7 using the results generated in `simulation_analysis.R'.
#This file is called by `simulation_master.R'.
##############################################################

#--------------------------
#functions used to calculate and display the results
#--------------------------

#function for bias, and corresponding MC standard error
bias=function(x,true){mean(x-true)}#x: estimated logHR
bias.mc=function(x){sqrt(sum((x-mean(x))^2)/(nsim*(nsim-1)))}#x: estimated logHR

#function for coverage, and corresponding MC standard error
coverage=function(x){mean(x)}#x: 0/1 coverage
coverage.mc=function(x){sqrt(coverage(x)*(1-coverage(x))/nsim)}

#function for test power, and corresponding MC standard error
test<-function(x){mean(1-x)}#x: 0/1 result of hypothesis test
test.mc=function(x){sqrt(test(1-x)*(1-test(1-x))/nsim)}

#function for empirical SE, and corresponding MC standard error
empse=function(x){sqrt(var(x))}#x: estimated logHR
empse.mc=function(x){sqrt(var(x))/sqrt(2*(nsim-1))}

#function for mean of the square root of the estimated variance, and corresponding MC standard error
estimatedse=function(x){sqrt(mean(x))}#x: estimated var of the estimated logHR
estimatedse.mc=function(x){sqrt((sum((x-mean(x)/nsim)^2)/(nsim-1))/(4*nsim*mean(x)))}#x: estimated var of the estimated logHR

#which file to put intermediate results in
filename <- paste0("nsim", nsim, "_n", n, "_nboot", nboot, "_copula.cor", copula.cor, "_propscore.mult", propscore.mult, ".out")

#--------------------------
#Table 1
#--------------------------

bias.naive<-sapply(1:5,function(x){bias(coef.naive[,x],true.coefs[x])})
bias.iptw<-sapply(1:5,function(x){bias(coef.iptw[,x],true.coefs[x])})

biasmc.naive<-sapply(1:5,function(x){bias.mc(coef.naive[,x])})
biasmc.iptw<-sapply(1:5,function(x){bias.mc(coef.iptw[,x])})

tab <- round(cbind(bias.naive, bias.iptw, biasmc.naive, biasmc.iptw), 3)
rownames(tab) <- paste0("beta", 1:5)
colnames(tab) <- c("naive", "IPTW", "naive.mcse", "IPTW.mcse")
write.table(file=paste0("intermediate_results/table1_", filename), x=tab)

#--------------------------
#Table 2
#--------------------------

cov.naive<-sapply(1:5,function(x){coverage(coverage.naive[,x])})
cov.sand<-sapply(1:5,function(x){coverage(coverage.sand[,x])})
cov.boot<-sapply(1:5,function(x){coverage(coverage.boot[,x])})

covmc.naive<-sapply(1:5,function(x){coverage.mc(coverage.naive[,x])})
covmc.sand<-sapply(1:5,function(x){coverage.mc(coverage.sand[,x])})
covmc.boot<-sapply(1:5,function(x){coverage.mc(coverage.boot[,x])})

tab <- round(cbind(cov.naive, cov.sand, cov.boot, covmc.naive, covmc.sand, covmc.boot), 3)
rownames(tab) <- paste0("beta", 1:5)
colnames(tab) <- c("naive", "sand", "boot", "naive.mcse", "sand.mcse", "boot.mcse")
write.table(file=paste0("intermediate_results/table2_", filename), x=tab)

#--------------------------
#Table 3
#--------------------------

testpower.naive<-test(test.naive)
testpower.sand<-test(test.sand)
testpower.boot<-test(test.boot)

testpowermc.naive<-test.mc(test.naive)
testpowermc.sand<-test.mc(test.sand)
testpowermc.boot<-test.mc(test.boot)

tab <- round(c(testpower.naive, testpower.sand, testpower.boot, testpowermc.naive, testpowermc.sand, testpowermc.boot), 3)
names(tab) <- c("naive", "sand", "boot", "naive.mcse", "sand.mcse", "boot.mcse")
write.table(file=paste0("intermediate_results/table3_", filename), x=tab)

#--------------------------
#Table A7
#--------------------------

empse.iptw<-sapply(1:5,function(x){empse(coef.iptw[,x])})

empsemc.iptw<-sapply(1:5,function(x){empse.mc(coef.iptw[,x])})

estimatedse.sand<-sapply(1:5,function(x){estimatedse(var.sand[,x])})
estimatedse.boot<-sapply(1:5,function(x){estimatedse(var.boot[,x])})

estimatedsemc.sand<-sapply(1:5,function(x){estimatedse.mc(var.sand[,x])})
estimatedsemc.boot<-sapply(1:5,function(x){estimatedse.mc(var.boot[,x])})

tab <- round(cbind(empse.iptw, estimatedse.sand, estimatedse.boot, empsemc.iptw, estimatedsemc.sand, estimatedsemc.boot), 3)
rownames(tab) <- paste0("beta", 1:5)
colnames(tab) <- c("naive", "sand", "boot", "naive.mcse", "sand.mcse", "boot.mcse")
write.table(file=paste0("intermediate_results/tableA7_", filename), x=tab)
