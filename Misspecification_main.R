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
library(VGAM)

library(expint) # For Gumbel distributed nu

source("Misspecification_functions.R")
source("Functions_ad.R")


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
psi.delta.theta <- c(-0.9, 0.97, 0.97, 1, 0.5)
par.skew_normal <- list(beta, eta, psi.delta.theta, gamma)

# Generate data with skew-normal errors to see what it looks like
output <- data.misspecified.skew(n, par.skew_normal, iseed, 1, 1)
data <- output[[1]]
errors <- output[[2]]
omega <- output[[3]]
lambda <- output[[5]]

# Make a plot of the errors
hist2(errors[,1], errors[,2], main = bquote(omega == .(omega)))

mean <- (2/pi)^(1/2)*psi.delta.theta[2]
curve(2*dnorm(x+mean,mean=0,sd=1)*pnorm(lambda[1]*(x+mean),mean=0,sd=1), xlim=c(-2,2),ylab="density", main="Skewed normal")



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

# Set degrees of freedom
par.df <- c(sd,3)
par.t <- list(beta, eta, par.df, gamma)


# Simulation
output_t <- Misspecification.t.sim(n, nsim = nsim, iseed = iseed, 
                                    par = par.t, Zbin = 1, Wbin = 1)


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

# Simulation
n <- 500
par.bimodal<-list(beta, eta, sd, gamma)
output_bimodal <- Misspecification.bimodal.sim(n, nsim = nsim, iseed = iseed, 
                                   par = par.bimodal, Zbin = 1, Wbin = 1)

#
# Heteroskedastic errors
#

par(mfrow = c(1,1))
errs=matrix(rep(1,2*n), nrow=n)
for (A in 1:n){
errs[A,] <- mvrnorm(1, mu=c(0, 0), Sigma=matrix(c(1.1*A,1.1*1.4*A*0.75,1.1*1.4*A*0.75,1.4*A),nrow=2,byrow=TRUE))
}
hist2(errs)
plot(errs[,1])

# Simulation
par.heteroscedastic<-list(beta, eta, sd, gamma)
output_heteroscedastic <- Misspecification.heteroscedastic.sim(n, nsim = nsim, iseed = iseed, 
                                                       par = par.heteroscedastic, Zbin = 1, Wbin = 1)

