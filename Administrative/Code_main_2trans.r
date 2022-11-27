# Clear workspace
rm(list = ls())

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)

# Make sure you are in the right directory (there is an equally named r-script
# in another directory).
getwd()
source("Functions_2trans.R")

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(1.8,0.9,0.5,-2.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

n <- 250
nsim <- 300
iseed <- 747852

# To split up simulation in different parts
number.of.parts <- 100
parts.to.evaluate <- 1:2
for (part.to.evaluate in parts.to.evaluate) {
  message("Evaluating part ", which(part.to.evaluate == parts.to.evaluate),
          " out of ", length(parts.to.evaluate))
  
  start.time <- Sys.time()
  SimulationCI11_SaraIlias_Simplified(n, nsim, iseed, init.value.theta_1,
                                      init.value.theta_2, part.to.evaluate,
                                      number.of.parts)
  
  diff <- start.time - Sys.time()
  message("This iteration ran for approximately ", round(diff/60), " minutes")
  message("")
}




