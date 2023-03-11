################################################################################
# This file should eventually be merged with the Code_main_ad.R and            #
# Functions_ad.R scripts.                                                      #
################################################################################

rm(list = ls()) # Remove this line later!

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(pbivnorm)
library(LaplacesDemon)
library(survival)
library(foreach)
library(doParallel)
source("Functions_ad.R")
source("Goodness-of-fit-test_functions.R")
source("Misspecification_functions.R")

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(1.8,0.9,0.5,-2.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")

n <- 1000
iseed <- 123
Zbin <- 1
Wbin <- 1
B <- 500
nruns <- 100

# Generate data and give it representative column names
data <- dat.sim.reg(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")

#
# Parametric bootstrap approach
#

# Get an idea of computing time of the full simulation. This function might take
# a while.
GOF_EstimateRunDuration(data, B, nruns, Zbin, Wbin, parallel = TRUE)

# Simulate the type-I error of the Goodness-Of-Fit test
typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns, B, iseed, Zbin, Wbin, parallel = TRUE)
reject90 <- typeIerror_results[["reject90"]]
reject95 <- typeIerror_results[["reject95"]]
write.csv(c(reject90, reject95), file = "typeIerror_iseed123_100runs_B500.csv")

typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns, B, iseed + 100, Zbin, Wbin, parallel = TRUE)
reject90 <- typeIerror_results[["reject90"]]
reject95 <- typeIerror_results[["reject95"]]
write.csv(c(reject90, reject95), file = "typeIerror_iseed223_100runs_B500.csv")

typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns, B, iseed + 200, Zbin, Wbin, parallel = TRUE)
reject90 <- typeIerror_results[["reject90"]]
reject95 <- typeIerror_results[["reject95"]]
write.csv(c(reject90, reject95), file = "typeIerror_iseed323_100runs_B500.csv")

typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns, B, iseed + 300, Zbin, Wbin, parallel = TRUE)
reject90 <- typeIerror_results[["reject90"]]
reject95 <- typeIerror_results[["reject95"]]
write.csv(c(reject90, reject95), file = "typeIerror_iseed423_100runs_B500.csv")

typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns, B, iseed + 400, Zbin, Wbin, parallel = TRUE)
reject90 <- typeIerror_results[["reject90"]]
reject95 <- typeIerror_results[["reject95"]]
write.csv(c(reject90, reject95), file = "typeIerror_iseed523_100runs_B500.csv")

# Plot the l'th rejected bootstrap sample together with the bootstrap T_{CM}.
l <- 3
hist(typeIerror_results[["TCMb_rejected"]][,l],
     main = "Histogram of bootstrap TCM leading to rejection",
     xlab = c(),
     xlim = c(0, max(max(typeIerror_results[["TCMb_rejected"]][l]), typeIerror_results[["TCM_rejected"]][l])))
abline(v = typeIerror_results[["TCM_rejected"]][l], col = "red")

# Combine results over 500 runs
first100 <- read.csv("typeIerror_iseed123_100runs.csv")
second100 <- read.csv("typeIerror_iseed223_100runs.csv")
third100 <- read.csv("typeIerror_iseed323_100runs.csv")
last200 <- read.csv("typeIerror_iseed423_200runs.csv")

results <- (first100 + second100 + third100 + 2*last200)[,2]/5
results <- data.frame(reject90 = results[1], reject95 = results[2])
results

#
# New approach
#

# This fails quite horribly. Based on the intermediate results, there are nearly
# no rejections on both levels.

n <- 1000
iseed <- 123
Zbin <- 1
Wbin <- 1
B <- 250
nruns <- 100

for (part in 1:3) {
  message("Busy doing part ", part, " out of 5.")
  message("")
  typeIerror_newapproach_results <- 
    GOF_SimulationTypeIerror_newapproach(parN, nruns, B, iseed + (part - 1)*nruns, Zbin, Wbin, parallel = TRUE)
  reject90 <- typeIerror_newapproach_results[["reject90"]]
  reject95 <- typeIerror_newapproach_results[["reject95"]]
  
  filename <- paste0("typeIerror_newapproach_iseed", iseed + (part - 1)*nruns, "_100runs.csv")
  
  write.csv(c(reject90, reject95), file = filename)
}



#
# omega-square distribution approach
#

# Get an idea of computing time of the full simulation without using bootstrap
# but comparing the statistic to the omega-square distribution instead.
GOF_noBootstrap_EstimateRunDuration(data, nruns, Zbin, Wbin)

# Simulate the type-I error of the Goodness-Of-Fit test which compares the
# obtained statistic to the omega-squared distribution.
typeIerror_results_noboot <- GOF_SimulationNoBootstrap_typeIerror(parN, nruns, iseed, Zbin, Wbin)
typeIerror_results_noboot[["reject90"]]
typeIerror_results_noboot[["reject95"]]
  # Using n = 1000, iseed = 123, B = 30 and nruns = 100, we obtain
  # reject90 = 0
  # reject95 = 0

################################################################################
#                   GOF test for misspecified models                           #
#                                                                              #
# Add this part to Section 4.2: Model under misspecification later and merge   #
# the code with Misspecification_functions.R                                   #
################################################################################

# Copied from 'Misspecification_main.R'... >>>..................................

################################################################################
#                              Set some parameters                             # 
################################################################################

n <- 1000
nsim <- 100
beta <- c(2.5, 2.6, 1.8, 2)
eta <- c(1.8,0.9,0.5,-2.2)
gamma <- c(-1, 0.6, 2.3)
sd <- c(1.1,1.4,0.75,1, 0.5)
iseed <- 12012000

parl <- length(beta)
totparl <- 2*parl
parlgamma <- (parl-1)

namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\alpha_T$","$\\lambda_T$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\alpha_C$","$\\lambda_C$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta_1$","$\\theta_2$")

# ...

################################################################################
#                       Misspecified error distributions                       # 
################################################################################

# 
# Skew-normal errors
#

# Set some more variables
psi.delta.theta <- c(-0.9, 0.97, 0.97, 1, 0.5)
par.skew_normal <- list(beta, eta, psi.delta.theta, gamma)

# ...

#
# Multivariate t distribution
#

par.df <- c(sd,3)
par.t <- list(beta, eta, par.df, gamma)

# ...

par.heteroscedastic <- list(beta, eta, sd, gamma)

# <<<...........................................................................

iseed <- 135468
type <- "heteroscedastic"
nruns <- 100
B <- 250
n <- 2000

for (part in 4:5) {
  message("Running part ", part, " out of 5.")
  message("")
  
  if (type == "t") {
    par.vector <- par.t 
  } else if (type == "skew"){
    par.vector <- par.skew_normal
  } else {
    par.vector <- par.heteroscedastic
  }
  misspec_results <- 
    GOF_SimulationMisspecification(type, par.vector, nruns, B, iseed + (part - 1)*100, 
                                   Zbin, Wbin, parallel = TRUE)
  reject90 <- misspec_results[["reject90"]]
  reject95 <- misspec_results[["reject95"]]
  
  filename <- paste0("misspec_", type, "_iseed", iseed + (part - 1)*100, "_",
                     nruns, "runs_B", B, ".csv")
  
  # put in proper folder
  if (type == "heteroscedastic") {
    folder <- paste0("Power/Misspec_hs_B", B, "_n", n, "/")
  } else {
    folder <- paste0("Power/Misspec_", type, "_B", B, "_n", n, "/")
  }
  
  write.csv(c(reject90, reject95), file = paste0(folder, filename))
  
}







