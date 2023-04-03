rm(list = ls())

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(mnormt)
library(matrixcalc)
library(cmprsk)
source("functions_competing_risks.R")

init.value.theta_1=1
init.value.theta_2=1
init.value.theta_3=1
parN = list(beta1=c(4,1,1.8,2),beta2=c(4.65,-1,-1.2,1.3),beta3=c(4.15,1.5,-1.5,-1.8),sd=c(1.5,1,1.2,-0.80,0.65,-0.75,1.5,0.5, 1),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 3*parl
parlgamma = (parl-1)
namescoef =  c("beta_{01}","beta_{11}","alpha_1","lambda_1","beta_{02}","beta_{12}","alpha_2","lambda_2","beta_{03}","beta_{13}","alpha_3","lambda_3","sigma_1","sigma_2","sigma_3","rho_12","rho_13","rho_23","theta_1","theta_2","theta_3")

samsize= c(1000)


# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

nsim = 500
myseed = 750751

for(l in samsize)
{
  message("sample size = ",l)
  results=SimulationCI22_Sara(l,nsim,myseed,init.value.theta_1, init.value.theta_2, init.value.theta_3) 
}

########### Results of parametric estimation

results.mean = results[[1]]
results.RMSE = results[[2]]

########### Nonparametric estimation

Time <- seq(from=1,to=200,by=0.5)


# Expected CIR (necessary to calculate RMSE)
C1R.11<-c()
C2R.11<-c()
C1R.01<-c()
C2R.01<-c()
C1R.10<-c()
C2R.10<-c()
C1R.00<-c()
C2R.00<-c()

XandW.11<-c(1,1,1)
Z.11<-1
V.true.11 <- (1-Z.11)*((1+exp(XandW.11%*%parN[[5]]))*log(1+exp(XandW.11%*%parN[[5]]))-(XandW.11%*%parN[[5]])*exp(XandW.11%*%parN[[5]]))-Z.11*((1+exp(-(XandW.11%*%parN[[5]])))*log(1+exp(-(XandW.11%*%parN[[5]])))+(XandW.11%*%parN[[5]])*exp(-(XandW.11%*%parN[[5]])))
M.true.11 <- c(1,1,Z.11,V.true.11)
for (i in 1:length(Time)){
  C1R.11[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.11) 
  C2R.11[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.11) 
}

XandW.10<-c(1,1,0)
Z.10<-1
V.true.10 <- (1-Z.10)*((1+exp(XandW.10%*%parN[[5]]))*log(1+exp(XandW.10%*%parN[[5]]))-(XandW.10%*%parN[[5]])*exp(XandW.10%*%parN[[5]]))-Z.10*((1+exp(-(XandW.10%*%parN[[5]])))*log(1+exp(-(XandW.10%*%parN[[5]])))+(XandW.10%*%parN[[5]])*exp(-(XandW.10%*%parN[[5]])))
M.true.10 <- c(1,1,Z.10,V.true.10)
for (i in 1:length(Time)){
  C1R.10[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.10) 
  C2R.10[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.10) 
}

XandW.01<-c(1,1,1)
Z.01<-0
V.true.01 <- (1-Z.01)*((1+exp(XandW.01%*%parN[[5]]))*log(1+exp(XandW.01%*%parN[[5]]))-(XandW.01%*%parN[[5]])*exp(XandW.01%*%parN[[5]]))-Z.01*((1+exp(-(XandW.01%*%parN[[5]])))*log(1+exp(-(XandW.01%*%parN[[5]])))+(XandW.01%*%parN[[5]])*exp(-(XandW.01%*%parN[[5]])))
M.true.01 <- c(1,1,Z.01,V.true.01)
for (i in 1:length(Time)){
  C1R.01[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.01) 
  C2R.01[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.01) 
}

XandW.00<-c(1,1,0)
Z.00<-0
V.true.00 <- (1-Z.00)*((1+exp(XandW.00%*%parN[[5]]))*log(1+exp(XandW.00%*%parN[[5]]))-(XandW.00%*%parN[[5]])*exp(XandW.00%*%parN[[5]]))-Z.00*((1+exp(-(XandW.00%*%parN[[5]])))*log(1+exp(-(XandW.00%*%parN[[5]])))+(XandW.00%*%parN[[5]])*exp(-(XandW.00%*%parN[[5]])))
M.true.00 <- c(1,1,Z.00,V.true.00)
for (i in 1:length(Time)){
  C1R.00[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.00) 
  C2R.00[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.00) 
}

results.C1R.11 = matrix(rep(unlist(C1R.11),nsim), nrow=nsim, byrow=TRUE)
results.C2R.11 = matrix(rep(unlist(C2R.11),nsim), nrow=nsim, byrow=TRUE)
results.C1R.10 = matrix(rep(unlist(C1R.10),nsim), nrow=nsim, byrow=TRUE)
results.C2R.10 = matrix(rep(unlist(C2R.10),nsim), nrow=nsim, byrow=TRUE)
results.C1R.01 = matrix(rep(unlist(C1R.01),nsim), nrow=nsim, byrow=TRUE)
results.C2R.01 = matrix(rep(unlist(C2R.01),nsim), nrow=nsim, byrow=TRUE)
results.C1R.00 = matrix(rep(unlist(C1R.00),nsim), nrow=nsim, byrow=TRUE)
results.C2R.00 = matrix(rep(unlist(C2R.00),nsim), nrow=nsim, byrow=TRUE)


# Non-parametric estimator of CIF
results.C1C.11 <-c()
results.C2C.11<-c()
results.C1C.10<-c()
results.C2C.10<-c()
results.C1C.01<-c()
results.C2C.01<-c()
results.C1C.00<-c()
results.C2C.00<-c()

for (i in 1:nsim){
  data = dat.sim.reg(samsize,parN,myseed+i,2,2)
  data1 = data[data[,7]==1,] #value of X is equal to 1
  Y = exp(data1[,1])
  d1 = data1[,2]
  d2 = data1[,3]
  d3 = data1[,4]
  da = data1[,5]
  D = d1+2*d2 #1 of C1, 2 if C2 and 0 if censored (dependent or administrative)
  Z = data1[,parl+4]
  W = data1[,parl+5]
  group=factor(Z+2*W) #0 (Z=0,W=0), 1(Z=1,W=0), 2(Z=0,W=1), 3(Z=1,W=1)
  
  fit = timepoints(cuminc(Y,D,group,cencode=0),Time) # nonparametric estimation
  
  results.C1C.11 = rbind(results.C1C.11,fit$est[4,])
  results.C2C.11 = rbind(results.C2C.11,fit$est[8,])
  results.C1C.10 = rbind(results.C1C.10,fit$est[2,])
  results.C2C.10 = rbind(results.C2C.10,fit$est[6,])
  results.C1C.01 = rbind(results.C1C.01,fit$est[3,])
  results.C2C.01 = rbind(results.C2C.01,fit$est[7,])
  results.C1C.00 = rbind(results.C1C.00,fit$est[1,])
  results.C2C.00 = rbind(results.C2C.00,fit$est[5,])
  
}

# results.mean2 -> mean estimates CIF non-parametric
results.mean2<-Time
results.mean2<-rbind(results.mean2,apply(results.C1C.11,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C2C.11,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C1C.10,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C2C.10,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C1C.01,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C2C.01,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C1C.00,2,mean))
results.mean2<-rbind(results.mean2,apply(results.C2C.00,2,mean))

rownames(results.mean2) <- c("Time","C1.11","C2.11","C1.10","C2.10","C1.01","C2.01","C1.00","C2.00")

#results.RMSE2 -> RMSE of non-parametric estimation
results.RMSE2<-Time
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C1C.11-results.C1R.11)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C2C.11-results.C2R.11)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C1C.10-results.C1R.10)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C2C.10-results.C2R.10)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C1C.01-results.C1R.01)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C2C.01-results.C2R.01)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C1C.00-results.C1R.00)^2,2,mean)))
results.RMSE2<-rbind(results.RMSE2,sqrt(apply((results.C2C.00-results.C2R.00)^2,2,mean)))

rownames(results.RMSE2) <- c("Time","C1.11","C2.11","C1.10","C2.10","C1.01","C2.01","C1.00","C2.00")


############## Plot curves

par(mfrow=c(2,2))
plot(c(0,results.mean[1,]),c(0,results.mean[2,]),type="l", lty=1, main=expression(paste(C1-Z1-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,30), ylim=c(0,0.5))
lines(c(0,results.mean2[1,]),c(0,results.mean2[2,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[6,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[8,]), lty=1,type="l", col="red")

plot(c(0,results.mean[1,]),c(0,results.mean[3,]),type="l", lty=1, main=expression(paste(C2-Z1-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,30))
lines(c(0,results.mean2[1,]),c(0,results.mean2[3,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[7,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[9,]), lty=1,type="l", col="red")

plot(c(0,results.mean[1,]),c(0,results.mean[26,]),type="l", lty=1, main=expression(paste(C1-Z0-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,100), ylim=c(0,0.8))
lines(c(0,results.mean2[1,]),c(0,results.mean2[8,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[30,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[32,]), lty=1,type="l", col="red")

plot(c(0,results.mean[1,]),c(0,results.mean[27,]),type="l", lty=1, main=expression(paste(C2-Z0-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,100), ylim=c(0,0.05))
lines(c(0,results.mean2[1,]),c(0,results.mean2[9,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[31,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[33,]), lty=1,type="l", col="red")
par(mfrow=c(1,1))



par(mfrow=c(2,2))
plot(c(0,results.mean[1,]),c(0,results.mean[18,]),type="l", lty=1, main=expression(paste(C1-Z0-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[20,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[22,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[24,]), lty=1,type="l", col="red")

plot(results.mean[1,],results.mean[19,],type="l", lty=1, main=expression(paste(C2-Z0-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[21,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[23,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[25,]), lty=1,type="l", col="red")

plot(c(0,results.mean[1,]),c(0,results.mean[10,]),type="l", lty=1, main=expression(paste(C1-Z1-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[12,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[14,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[16,]), lty=1,type="l", col="red")

plot(c(0,results.mean[1,]),c(0,results.mean[11,]),type="l", lty=1, main=expression(paste(C2-Z1-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[13,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[15,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[17,]), lty=1,type="l", col="red")
par(mfrow=c(1,1))

# look at results
results.mean[,c(19,39,59,99,199)]
results.mean2[,c(19,39,59,99,199)]
results.RMSE[,c(19,39,59,99,199)]
results.RMSE2[,c(19,39,59,99,199)]


# Global RMSE

grid.lower.bound <- 1
grid.upper.bound1 <- 30
grid.upper.bound2 <- 100

number.of.grid.cells1 <- (grid.upper.bound1-grid.lower.bound)*2+1
number.of.grid.cells2 <- (grid.upper.bound2-grid.lower.bound)*2+1
cell.width <- 0.5


int.C1.11 <- 0
int.C2.11 <- 0
int.C1E.11 <- 0
int.C2E.11 <- 0
int.C1I.11 <- 0
int.C2I.11 <- 0

# RMSE for 11
for (i in 1:number.of.grid.cells1) {
  int.C1.11 <- int.C1.11 + results.RMSE[2,i]*(cell.width)
  int.C2.11 <- int.C2.11 + results.RMSE[3,i]*(cell.width)
  int.C1E.11 <- int.C1E.11 + results.RMSE[6,i]*(cell.width)
  int.C2E.11 <- int.C2E.11 + results.RMSE[7,i]*(cell.width)
  int.C1I.11 <- int.C1I.11 + results.RMSE2[2,i]*(cell.width)
  int.C2I.11 <- int.C2I.11 + results.RMSE2[3,i]*(cell.width)
}


int.C1.00 <- 0
int.C2.00 <- 0
int.C1E.00 <- 0
int.C2E.00 <- 0
int.C1I.00 <- 0
int.C2I.00 <- 0

# RMSE for 11
for (i in 1:number.of.grid.cells2) {
  int.C1.00 <- int.C1.00 + results.RMSE[20,i]*(cell.width)
  int.C2.00 <- int.C2.00 + results.RMSE[21,i]*(cell.width)
  int.C1E.00 <- int.C1E.00 + results.RMSE[24,i]*(cell.width)
  int.C2E.00 <- int.C2E.00 + results.RMSE[25,i]*(cell.width)
  int.C1I.00 <- int.C1I.00 + results.RMSE2[8,i]*(cell.width)
  int.C2I.00 <- int.C2I.00 + results.RMSE2[9,i]*(cell.width)
}






