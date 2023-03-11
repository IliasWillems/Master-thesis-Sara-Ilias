# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
chess <- read.csv("duration_df_new.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)

chess <- chess[!(is.na(chess$time)),]
# time
Y <- as.matrix(log(chess$time))
# Administrative censoring
# Dependent censoring --> won by resignation (if someone was clearly winning)
Delta <- as.matrix(chess$delta)

chess$intercept <- rep(1, nrow(chess))

# The variable indicating which opening was played is categorical and has a lot
# of levels, most of which only occur once. Therefore, we reclassify al
# uncommon openings (i.e. openings that are not played more than 2 times) under
# the same name

chess$opening_cat <- ifelse(chess$opening=="1. e4 e5 2. Nf3 Nf6 3. Nxe5 Nc6 " |
                        grepl("1. e4 c6 2. d4 d5", chess$opening, fixed=TRUE)==TRUE,1,
                      ifelse(grepl("1. e4", chess$opening, fixed=TRUE)==TRUE,2,3))



#################### Multiple EV ###############################################

source("Functions_multipleEV.r")

namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}","alpha_{T,1}","alpha_{T,2}",
               "lambda_{T_1}","lambda_{T_2}","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "alpha_{C,1}","alpha_{C,2}","lambda_{C_1}","lambda_{C_2}","sigma_T",
               "sigma_C","rho","theta_1","theta_2")



# Define some useful variables
parl <- ncol(X) + 4
totparl <- 2*parl
parlgamma <- ncol(X)+1


# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, eval_depth, poi_is_white].

# NOTE: as soon as data set is fixed, we should also include ELO of white and
#       black player.
X = as.matrix(subset(chess, select = c(ncol(chess)-1, 5, 8)))


# Create data matrix of confounded variable. Z is the indicator of which opening
# was played in the game. Z is a multi-level categorical variable. This means
# that we need more than one dummy variable to represent it. 
Z = as.matrix(chess$opening_cat)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W = as.matrix(chess$opponent_exp)

# Create data matrix
XandW = as.matrix(cbind(X,W))
n=nrow(chess)
data=as.matrix(cbind(Y,Delta,X,Z,W))


init.value.theta_1 <- 1
init.value.theta_2 <- 2
DataApplicationChess.multipleEV(data, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run



#################### Stratification#############################################

source("Functions_ad.r")

namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}","alpha_{T}",
               "lambda_{T}","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "alpha_{C}","lambda_{C}","sigma_T",
               "sigma_C","rho","theta_1","theta_2")



# Define some useful variables
parl <- ncol(X) + 2
totparl <- 2*parl
parlgamma <- ncol(X)+1


chess$Z1 = ifelse(chess$opening_cat==2,1,0)
chess$Z2 = ifelse(chess$opening_cat==3,1,0)
chess.1 <- chess[(chess$Z1==1 | (chess$Z1==0 & chess$Z2==0)),]
chess.2 <- chess[(chess$Z2==1 | (chess$Z1==0 & chess$Z2==0)),]


# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, eval_depth, poi_is_white].

# NOTE: as soon as data set is fixed, we should also include ELO of white and
#       black player.
X.1 = as.matrix(subset(chess.1, select = c(ncol(chess)-3, 5, 8)))
X.2 = as.matrix(subset(chess.2, select = c(ncol(chess)-3, 5, 8)))

# Create data matrix of confounded variable. Z is the indicator of which opening
# was played in the game. Z is a multi-level categorical variable. This means
# that we need more than one dummy variable to represent it. 
Z.1 = as.matrix(chess.1$Z1)
Z.2 = as.matrix(chess.2$Z2)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W.1 = as.matrix(chess.1$opponent_exp)
W.2 = as.matrix(chess.2$opponent_exp)


# observed variables
Y.1 <- as.matrix(log(chess.1$time))
Y.2 <- as.matrix(log(chess.2$time))
Delta.1 <- as.matrix(chess.1$delta)
Delta.2 <- as.matrix(chess.2$delta)

# Create data matrix
XandW.1 = as.matrix(cbind(X.1,W.1))
n.1=nrow(chess.1)
data.1=as.matrix(cbind(Y.1,Delta.1,X.1,Z.1,W.1))

XandW.2 = as.matrix(cbind(X.2,W.2))
n.2=nrow(chess.2)
data.2=as.matrix(cbind(Y.2,Delta.2,X.2,Z.2,W.2))

init.value.theta_1 <- 1
init.value.theta_2 <- 2
DataApplicationChess(data.1, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run
DataApplicationChess(data.2, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run



