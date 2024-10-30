##############################################################
#This file applies the analyses described in Section 6.1 to each simulated data set and stores the results.
#This file is called by `simulation_master.R'.
##############################################################

#---
#storage for simulation results

coef.naive<-matrix(nrow=nsim,ncol=5)
coef.iptw<-matrix(nrow=nsim,ncol=5)

var.sand<-matrix(nrow=nsim,ncol=5)
var.boot<-matrix(nrow=nsim,ncol=5)

coverage.naive<-matrix(nrow=nsim,ncol=5)
coverage.sand<-matrix(nrow=nsim,ncol=5)
coverage.boot<-matrix(nrow=nsim,ncol=5)

test.naive<-rep(NA,nsim)
test.sand<-rep(NA,nsim)
test.boot<-rep(NA,nsim)

#---
#generating random seeds, enabling us to reproduce the results in the paper
#Note that we do not just use a single seed because the analyses for the paper 
# were performed using a high performance cluster in which simulations were run in parallel rather than sequentially

set.seed(2611)
seed.vector=sample(1:999999,1000) 

#---
#start of simulation loop

for (s in 1:nsim){

  print(paste0("simulation ",s))
  
  set.seed(seed.vector[s])
  
  #--------------
  #generate data
  #--------------
  
  source("simulation_data_generate.R")
  
  #--------------
  #Estimate stabilised weights
  #The weights models are fitted using data pooled across visits
  #--------------
  
  #fit models: uses logistic regression and is pooled across time points
  wtmod.denom<-glm(A~as.factor(visit.time)+A_lag+M1+M2+B1+X1+X2,data=dat.long,family="binomial")
  wtmod.numer<-glm(A~as.factor(visit.time)+A_lag+M1+M2,data=dat.long,family="binomial")
  
  #get fitted values from weights models
  wtmod.denom.fitted<-predict(wtmod.denom,type="response")
  wtmod.numer.fitted<-predict(wtmod.numer,type="response")
  
  #obtain weights at each time point
  dat.long$wt.time<-ifelse(dat.long$A==1,
                           wtmod.numer.fitted/wtmod.denom.fitted,
                           (1-wtmod.numer.fitted)/(1-wtmod.denom.fitted))
  
  #take cumulative product over visit.time to give the IPTW
  dat.long<-dat.long%>%group_by(id)%>%mutate(iptw=cumprod(wt.time))
  
  #--------------
  #Fit Cox model using stabilised weights
  #--------------
  
  cox.iptw<-coxph(Surv(tstart,tstop,fail.updated)~M1+M2+A+M1*A+M2*A,data=dat.long,cluster=id,weights = iptw)    
  
  #--------------
  #Fit Cox model without any weights (naive analysis)
  #--------------
  
  cox.naive<-coxph(Surv(tstart,tstop,fail.updated)~M1+M2+A+M1*A+M2*A,data=dat.long,cluster=id)    
  
  #--------------
  #store results
  #--------------
  
  #coefficients
  coef.naive[s,]<-cox.naive$coef
  coef.iptw[s,]<-cox.iptw$coef
  
  #sandwich variances
  var.sand[s,]<-diag(cox.iptw$var)

  #indicator of whether each 95% CI contains the true value (for 95% coverage)
  coverage.naive[s,]<-sapply(1:5,FUN=function(x){true.coefs[x]>=confint(cox.naive)[x,1] & true.coefs[x]<=confint(cox.naive)[x,2]})
  coverage.sand[s,]<-sapply(1:5,FUN=function(x){true.coefs[x]>=confint(cox.iptw)[x,1] & true.coefs[x]<=confint(cox.iptw)[x,2]})

  #indicator of whether 95% CI for M1*A coefficient contains 0 (for type II errors)
  test.naive[s]<-(0>=confint(cox.naive)[4,1] & 0<=confint(cox.naive)[4,2])
  test.sand[s]<-(0>=confint(cox.iptw)[4,1] & 0<=confint(cox.iptw)[4,2])

  #--------------
  #Now we perform the analysis across bootstrap samples
  #--------------
  
  if(include.bootstrap==1){
  
    coef.boot<-matrix(NA,nrow=nboot,ncol=5) #storage for bootstrap estimates
    
    for(b in 1:nboot){
      
      if (b==1|(b%%100==0)) {
        print(paste0("bootstrap ",b))}
      #---
      #Take bootstrap sample of individuals in dat, then recreate dat.long
      #Note that we need to generate new ID number when people are repeated in the bootstrap sample
      #---
      
      boot.samp<-sample(dat$id,n,replace=T)
      dat.boot<-dat[boot.samp,]
      dat.boot$id<-1:n
      
      dat.boot.long<-reshape(dat.boot,
                             varying=c(paste0("A.",1:10),paste0("X1.",1:10),paste0("X2.",1:10)),
                             timevar="visit.time",idvar="id",direction="long")
      dat.boot.long$visit.time<-dat.boot.long$visit.time-1
      
      dat.boot.long<-dat.boot.long[order(dat.boot.long$id,dat.boot.long$visit.time),]
      
      #delete rows after the person's event/censoring time
      dat.boot.long<-dat.boot.long[dat.boot.long$time>=dat.boot.long$visit.time,]
      
      #generate lagged treatment
      dat.boot.long<-dat.boot.long%>%group_by(id)%>%mutate(A_lag=lag(A,default=0))
      
      #generate row number and total number of rows
      dat.boot.long<-dat.boot.long%>%group_by(id)%>%mutate(rownum=row_number())%>%mutate(maxrow=max(rownum))
      
      #generate start and stop times in each row, and corresponding time-updated failure indicator
      dat.boot.long$tstart<-dat.boot.long$visit.time
      
      dat.boot.long$tstop<-ifelse(dat.boot.long$rownum!=dat.boot.long$maxrow,dat.boot.long$visit.time+1,dat.boot.long$time)
      
      dat.boot.long$fail.updated<-ifelse(dat.boot.long$rownum!=dat.boot.long$maxrow,0,dat.boot.long$fail)
      
      #--------------
      #Estimate stabilised weights in bootstrap sample
      #The weights models are fitted using data pooled across visits
      #--------------
      
      #fit models: uses logistic regression and is pooled across time points
      wtmod.denom.boot<-glm(A~as.factor(visit.time)+A_lag+M1+M2+B1+X1+X2,data=dat.boot.long,family="binomial")
      wtmod.numer.boot<-glm(A~as.factor(visit.time)+A_lag+M1+M2,data=dat.boot.long,family="binomial")
      
      #get fitted values from weights models
      wtmod.denom.boot.fitted<-predict(wtmod.denom.boot,type="response")
      wtmod.numer.boot.fitted<-predict(wtmod.numer.boot,type="response")
      
      #obtain weights at each time point
      dat.boot.long$wt.boot.time<-ifelse(dat.boot.long$A==1,
                                         wtmod.numer.boot.fitted/wtmod.denom.boot.fitted,
                                         (1-wtmod.numer.boot.fitted)/(1-wtmod.denom.boot.fitted))
      
      #take cumulative product over visit.time to give the IPTW
      dat.boot.long<-dat.boot.long%>%group_by(id)%>%mutate(iptw.ReEst=cumprod(wt.boot.time))
      
      #--------------
      #obtain stabilised weights based on weights models fitted in the original complete data (i.e. not in the bootstrap sample)
      #The weights models are fitted using data pooled across visits
      #--------------
      
      #get fitted values from weights models
      wtmod.denom.fitted<-predict(wtmod.denom,newdata=dat.boot.long,type="response")
      wtmod.numer.fitted<-predict(wtmod.numer,newdata=dat.boot.long,type="response")
      
      #obtain weights at each time point
      dat.boot.long$wt.time<-ifelse(dat.boot.long$A==1,
                                    wtmod.numer.fitted/wtmod.denom.fitted,
                                    (1-wtmod.numer.fitted)/(1-wtmod.denom.fitted))
      
      #take cumulative product over visit.time to give the IPTW
      dat.boot.long<-dat.boot.long%>%group_by(id)%>%mutate(iptw=cumprod(wt.time))
      
      #--------------
      #Fit Cox model using stabilised weights, with weights models fitted in bootstrap sample
      #--------------
      
      cox.iptw.boot<-coxph(Surv(tstart,tstop,fail.updated)~M1+M2+A+M1*A+M2*A,data=dat.boot.long,weights = iptw.ReEst)  
      
      #--------------
      #store results from this bootstrap sample
      #--------------
      
      coef.boot[b,]<-cox.iptw.boot$coefficients
    }
    
    #--------------
    #obtain results from bootstrapping
    #--------------
    
    #variances
    var.boot[s,]<-sapply(1:5,FUN=function(x){var(coef.boot[,x])})

    #indicator of whether each 95% CI contains the true value (for 95% coverage)
    coverage.boot[s,]<-sapply(1:5,FUN=function(x){true.coefs[x]>=quantile(coef.boot[,x],0.025) & true.coefs[x]<=quantile(coef.boot[,x],0.975)})
    
    #indicator of whether 95% CI for M1*A coefficient contains 0 (for type II errors)
    test.boot[s]<-(0>=quantile(coef.boot[,4],0.025) & 0<=quantile(coef.boot[,4],0.975))
    
  }
  
}
