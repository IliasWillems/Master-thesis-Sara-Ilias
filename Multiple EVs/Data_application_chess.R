# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
chess <- read.csv("duration_df.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
source("Functions.r")

Y <- as.matrix(log(chess$time))

# NOTE: It is possible to adapt the python code that scraped the data to provide
#       info about the type of censoring (administrative or not).
Delta <- as.matrix(chess$delta)

chess$intercept <- rep(1, nrow(chess))

# The variable indicating which opening was played is categorical and has a lot
# of levels, most of which only occur once. Therefore, we reclassify al
# uncommon openings (i.e. openings that are not played more than 2 times) under
# the same name
counts <- table(chess$opening)
often_played_openings <- names(counts[counts > 2])
chess$opening <- ifelse(chess$opening %in% often_played_openings, chess$opening,
                        "uncommon opening")
counts2 <- table(chess$opening)
counts2

# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, eval_depth, poi_is_white].

# NOTE: as soon as data set is fixed, we should also include ELO of white and
#       black player.
X = as.matrix(subset(chess, select = c(ncol(chess), 5, 8)))

# Define some useful variables
parl <- ncol(X) + 2
totparl <- 2*parl
parlgamma <- parl - 1

# Create data matrix of confounded variable. Z is the indicator of which opening
# was played in the game. Z is a multi-level categorical variable. This means
# that we need more than one dummy variable to represent it. Our model is not
# (yet) able to handle this.
Z = as.matrix(chess$opening)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W = as.matrix(chess$opponent_exp)

# Create data matrix
XandW = as.matrix(cbind(X,W))
n=nrow(dataset)
data=as.matrix(cbind(Y,Delta,X,Z,W))

namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}","beta_{T,3}","beta_{T,5}",
               "alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "beta_{C,3}","beta_{C,5}","alpha_C","lambda_C","sigma_T",
               "sigma_C","rho","theta")

namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\beta_{T,2}$",
               "$\\beta_{T,3}$","$\\beta_{T,5}$","$\\alpha_T$","$\\lambda_T$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\beta_{C,2}$",
               "$\\beta_{C,3}$","$\\beta_{C,5}$","$\\alpha_C$","$\\lambda_C$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta$")

# We follow the approach as in Chapter 5 of Crommen, Van Keilegom (2022). Hence,
# we will be comparing the models (1) assuming independence (2) assuming
# dependence but no confounding and (3) our model. It would not make sense to
# compare our model with the model of Crommen, Van Keilegom (2022) directly.
# However, we will compare the estimated survival curves.


init.value.theta <- 1
DataApplicationJPTA(data, init.value.theta) # Takes about 1-2 minute to run


