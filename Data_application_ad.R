
# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
fulldata<-read.csv("clean_dataset_JTPA.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
source("Functions_ad.r")

# Select all males that have children (fathers)
# dataset = fulldata[(fulldata$children == 1 & fulldata$married == 0 & fulldata$male==0 & fulldata$white==0), -c(1,4,5,7,8)]

dataset = fulldata[(fulldata$children == 1 & fulldata$married == 0 & fulldata$male==0 & fulldata$white==0), -c(1,4,5,7,8)]

Y = as.matrix(log(dataset$days))
Delta = as.matrix(dataset$delta)
dataset$intercept=rep(1,nrow(dataset))

# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, age, hsged, white, married]
X = as.matrix(subset(dataset, select = c(ncol(dataset),1,2)))

# Define some useful variables
parl=ncol(X)+2
totparl=2*parl
parlgamma=parl-1

# Create data matrix of confounded variable. Z is the indicator whether the
# participant actually participated in the study. (0 for no participation, 1 
# otherwise). 
Z = as.matrix(dataset$jtpa)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W = as.matrix(dataset$treatment)

# Create data matrix
XandW = as.matrix(cbind(X,W))
n=nrow(dataset)
data=as.matrix(cbind(Y,Delta,X,Z,W))
namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}",
               "alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "alpha_C","lambda_C","sigma_T",
               "sigma_C","rho","theta_1","theta_2")


# Fitting logistic regression model
fit_log<-glm(Delta ~X[,2]+X[,3], family="binomial")
summary(fit_log)

# We follow the approach as in Chapter 5 of Crommen, Van Keilegom (2022). Hence,
# we will be comparing the models (1) assuming independence (2) assuming
# dependence but no confounding and (3) our model. It would not make sense to
# compare our model with the model of Crommen, Van Keilegom (2022) directly.
# However, we will compare the estimated survival curves.


init.value.theta_1 <- 1
init.value.theta_2 <- 2
DataApplicationJPTA(data, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run










