library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
source("Functions_ad_splitUp_Simulations.R")

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(1.8,0.9,0.5,-2.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)

namescoef = c("$\\beta_{T,0}$", "$\\beta_{T,1}$", "$\\alpha_T$", "$\\lambda_T$",
              "$\\beta_{C,0}$", "$\\beta_{C,1}$", "$\\alpha_C$", "$\\lambda_C$",
              "$\\sigma_T$", "$\\sigma_C$", "$\\rho$", "$\\theta_1$",
              "$\\theta_2$")

samsize = c(250, 500, 1000)

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

nsim <- 2500
myseed <- 750751
number.of.parts <- 250
#126-250
parts.to.evaluate <- 76:86


for (part.to.evaluate in parts.to.evaluate) {
  for (n in samsize){
  # Start message
  message("Starting iteration ", part.to.evaluate - parts.to.evaluate[1] + 1,
          " out of ", length(parts.to.evaluate), " with sample size ",n)
  
  start.time <- Sys.time()
  
  #
  # Make sure to select the correct simulation function
  #
  SimulationCI11_splitup(n, nsim, myseed, init.value.theta_1, init.value.theta_2,
                         part.to.evaluate, number.of.parts)
  SimulationCI12_splitup(n, nsim, myseed, init.value.theta_1, init.value.theta_2,
                         part.to.evaluate, number.of.parts)
  SimulationCI21_splitup(n, nsim, myseed, init.value.theta_1, init.value.theta_2,
                         part.to.evaluate, number.of.parts)
  SimulationCI22_splitup(n, nsim, myseed, init.value.theta_1, init.value.theta_2,
                         part.to.evaluate, number.of.parts)
  
  # Tell the time this iteration ran
  diff <- difftime(Sys.time(), start.time, units = "mins")
  diff <- as.integer(round(diff))
  if (diff > 60) {
    message("This iteration took ", floor(diff / 60), " hours and ", diff %% 60, " minutes.")
  } else {
    message("This iteration took ", diff, " minutes.")
  }
  message("")
  }
}

summarize_results("CI12", 1000)




