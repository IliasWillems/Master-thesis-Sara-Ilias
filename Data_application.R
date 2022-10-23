
# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
fulldata<-read.csv("clean_dataset_JTPA.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
source("Functions.r")

# Select all males that have children (fathers)
dataset = fulldata[(fulldata$children == 1 & fulldata$male == 1), -c(1,5,8)]

Y = as.matrix(log(dataset$days))
Delta = as.matrix(dataset$delta)
dataset$intercept=rep(1,nrow(dataset))

# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, age, hsged, white, married]
X = as.matrix(subset(dataset, select = c(ncol(dataset),1,2,3,5)))

# Define some useful variables
parl=ncol(X)+2
totparl=2*parl
parlgamma=parl-1

# Create data matrix of confounded variable. Z is the indicator whether the
# participant actually participated the training program (0 for no participation,
# 1 otherwise).
Z = as.matrix(dataset$jtpa)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W = as.matrix(dataset$treatment)

# Create data matrix
XandW = as.matrix(cbind(X,W))
n=nrow(dataset)
data=as.matrix(cbind(Y,Delta,X,Z,W))
namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}","beta_{T,3}","beta_{T,5}",
               "alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "beta_{C,3}","beta_{C,5}","alpha_C","lambda_C","sigma_T",
               "sigma_C","rho","theta")


# We follow the approach as in Chapter 5 of Crommen, Van Keilegom (2022). Hence,
# we will be comparing the models (1) assuming independence (2) assuming
# dependence but no confounding and (3) our model. It would not make sense to
# compare our model with the model of Crommen, Van Keilegom (2022) directly.
# However, we will compare the estimated survival curves.


############################### IMPOARTANT NOTE: ###############################
# If in the estimation of the independence model, it should be                 #
# VargammaI = HgammaI[...] instead of VargammaI = Hgamma[...], this still needs#
# to be changed! Here we use VargammaI = Hgamma[...].                          #
################################################################################

init.value.theta <- 1
DataApplicationJPTA(data, init.value.theta) # Takes about 1 minute to run















