library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(nnet)

# For parallel computing
library(foreach)
library(doParallel)

source("Functions_multipleEV.R")

nEV <- 3-1
init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,-1,2,1.5),eta=c(1.8,0.9,0.5,1.5,-2.2,1.3),sd=c(1.1,1.4,0.75,1, 0.5),gamma1=c(-1,0.6,2.3,1.4),gamma2=c(1.3,2.6,-1.5,-1))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-nEV)
namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\alpha_{T,1}$","$\\alpha_{T,2}$","$\\lambda_{T,1}$","$\\lambda_{T,2}$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\alpha_{C,1}$","$\\alpha_{C,2}$","$\\lambda_{C,1}$","$\\lambda_{C,2}$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta_1$","$\\theta_2$")

samsize= c(500, 1000, 2000)

nsim <- 2500
myseed <- 472473
number.of.parts <- 250

# For now, we only evaluate 100 parts. We can extend this number later if
# possible.
parts.to.evaluate <- 241:250


MAX_NBR_CORES_TO_USE <- 10

for (part.to.evaluate in parts.to.evaluate) {
  
  # Prepare for parallel processing
  nbr.cores <- detectCores()
  cores.to.use <- min(MAX_NBR_CORES_TO_USE, nbr.cores - 1)
  
  clust <- makeCluster(cores.to.use)
  registerDoParallel(clust)
  
  for (n in samsize){
    
    # Start message
    message("Starting iteration ", part.to.evaluate - parts.to.evaluate[1] + 1,
            " out of ", length(parts.to.evaluate), " with sample size ",n)
    
    start.time <- Sys.time()
    
    # Simulation function here
    source("Functions_multipleEV.R")
    SimulationCI22_MEV(n, nsim, myseed, init.value.theta_1, init.value.theta_2,
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
  
  # Stop cluster
  stopCluster(clust)
}

# Summarize results into latex table
Summarize_results_full_MEV("CI22")
  # Add "{width=\linewidth}" after "\begin{adjustbox}", otherwise the resulting
  # LaTeX code won't compile.






