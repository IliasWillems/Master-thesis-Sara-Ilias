library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
source("Functions_multipleEV.R")

nEV <- 3-1
init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,-1,2,1.5),eta=c(1.8,0.9,0.5,1.5,-2.2,1.3),sd=c(1.1,1.4,0.75,1, 0.5),gamma1=c(-1,0.6,2.3,1.4),gamma2=c(1.3,2.6,-1.5,-1))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-nEV)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_{T,1}","alpha_{T,2}","lambda_{T,1}","lambda_{T,2}",
               "beta_{C,0}","beta_{C,1}","alpha_{C,1}","alpha_{C,2}","lambda_{C,1}","lambda_{C,2}",
               "sigma_T","sigma_C","rho","theta_1","theta_2")

samsize= c(250, 500, 1000)

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

for(l in samsize)
{
  nsim = 10
  myseed = 876661
  message("sample size = ",l)
  SimulationCI22_EV(l,nsim,myseed,init.value.theta_1, init.value.theta_2)
}

