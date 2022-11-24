# Some notes:
#
# The code is now only (partially) done for the skew-normal dist., the other two
# distributions as well as misspecification of the control function is still to
# come.


rm(list = ls())

library(MASS)
library(squash)
library(LaplacesDemon)
library(nloptr)
library(numDeriv)
library(xtable)

library(expint) # For Gumbel distributed nu

source("Misspecification_functions.R")
source("Functions.R")


################################################################################
#                              Set some parameters                             # 
################################################################################

n <- 500
nsim <- 200
beta <- c(2.5, 2.6, 1.8, 2)
eta <- c(2.8, 1.9, 1.5, 1.2)
gamma <- c(-1, 0.6, 2.3)
sd <- c(1.1, 1.4, 0.75, 1)
iseed <- 12012000

parl <- length(beta)
totparl <- 2*parl
parlgamma <- (parl-1)

namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\alpha_T$","$\\lambda_T$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\alpha_C$","$\\lambda_C$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta$")

################################################################################
#                        Misspecified control function                         # 
################################################################################

# nu ~ N(0, 1)

par.probit = list(beta, eta, sd, gamma)

Misspecificaion.probit.sim(n = n, nsim = nsim, iseed = iseed, par = par.probit,
                           Wbin = 1)


# nu ~ Gumbel(0, 1)

par.cloglog = list(beta, eta, sd, gamma)

Misspecificaion.cloglog.sim(n = n, nsim = nsim, iseed = iseed, par = par.cloglog,
                            Wbin = 1)

################################################################################
#                       Misspecified error distributions                       # 
################################################################################

# 
# Skew-normal errors
#

# Set some more variables
psi.delta.theta <- c(-0.9, 0.8, 0.8, 0.5)
par.skew_normal <- list(beta, eta, psi.delta.theta, gamma)

# Generate data with skew-normal errors to see what it looks like
output <- misspecified.skew(n, par.skew_normal, iseed, 1, 1)
data <- output[[1]]
errors <- output[[2]]
omega <- output[[3]]

# Make a plot of the errors
hist2(errors[,1], errors[,2], main = bquote(omega == .(omega)))

# Simulation
output <- Misspecification.skew.sim(n, nsim = nsim, iseed = iseed, 
                                    par = par.skew_normal, Zbin = 1, Wbin = 1)


#
# Multivariate t distribution
#

par(mfrow = c(1, 2))
mu <- c(0, 0)
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

errs <- rmvt(n = n, mu, S = sigma, df = 10)
hist2(errs[,1], errs[,2])

# Compare with normal distribution
errs <- mvrnorm(n, mu=mu , Sigma=sigma)
hist2(errs[,1], errs[,2])

# By taking the degrees of freedom to be small, the t-distribution will have 
# more variance. Taking df = Inf, the variance-covariance matrix is equal to
# that of the normal distribution.


#
# Bimodal normal distribution
#

par(mfrow = c(1,1))
n <- 100000
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
errs1 <- mvrnorm(n/2, mu=c(-1, 1), Sigma=sigma)
errs2 <- mvrnorm(n/2, mu=c(1, -1), Sigma=sigma)
errs <- rbind(errs1, errs2)
hist2(errs[,1], errs[,2])


#
# Heteroskedastic errors?
#



