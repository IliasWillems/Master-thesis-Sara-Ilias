library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)

# Make sure that the correct 'Functions_2trans.r'-script is loaded (not the one
# taking administrative censoring into account).
getwd()
source("Functions_2trans.r")

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(2.8,1.9,1.5,1.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")

samsize= c(10000) #c(250, 500, 1000)

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

for(l in samsize)
{
  nsim = 300
  myseed = 876661
  message("sample size = ",l)
  SimulationCI11_SaraIlias_simplified(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 1
  #SimulationCI12_SaraIlias(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 2
  #SimulationCI21_SaraIlias(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 3
  #SimulationCI22_SaraIlias(l,nsim,myseed,init.value.theta_1, init.value.theta_2) # Design 4
}
