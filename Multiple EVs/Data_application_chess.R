# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
wd <- getwd()
parent_folder <- substr(wd, 1, nchar(wd) - 13)

chess <- read.csv(paste0(parent_folder, "/Data/duration_df_new.csv"))
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(dplyr)
library(nnet)

# Remove missing observations
chess <- chess[!(is.na(chess$time)),]
plot(chess[chess$delta==1,]$white_elo,log(chess[chess$delta == 1, ]$time))
plot(chess[chess$delta==1,]$black_elo,log(chess[chess$delta == 1, ]$time))

# Remove a couple of outliers
chess <- chess[chess$white_elo>=1000,]
chess <- chess[chess$black_elo>=1000,]

plot(chess[chess$delta==1,]$white_elo,log(chess[chess$delta == 1, ]$time))
plot(chess[chess$delta==1,]$black_elo,log(chess[chess$delta == 1, ]$time))

# Check for outliers
plot(log(chess[chess$delta == 1,]$time), rep(1.05, sum(chess$delta == 1)),
     main = "Plot of time until checkmate",
     xlab = "log time (in seconds)",
     ylab = "",
     yaxt = "n",
     xlim=c(0,12))
points(log(chess[chess$delta == 0,]$time), rep(1, sum(chess$delta == 0)), col = "red")
points(log(chess[chess$delta == 2,]$time), rep(0.95, sum(chess$delta == 2)), col = "green")
# chess <- chess[chess$time < 1500, ]

# Only work with games that took longer than 2 minutes to be sure that we
# exclude games that were played in the 1+0 format.
TWO_MINUTES <- 60*2
chess <- chess[which(chess$time > TWO_MINUTES), ]

n <- nrow(chess)

# time
Y <- as.matrix(log(chess$time))

# standardize
chess$opponent_exp <- (chess$opponent_exp-mean(chess$opponent_exp))/sd(chess$opponent_exp)
chess$white_elo <- (chess$white_elo-mean(chess$white_elo))/sd(chess$white_elo)
chess$black_elo <- (chess$black_elo-mean(chess$black_elo))/sd(chess$black_elo)

# Even though the data has a column with a censoring indicator, we update it to
# make it more accurate.

# Administrative censoring: "draw", "game abandoned", "unknown end", "won by resignation" (if neither player was clearly winning)
# Dependent censoring: "won by resignation" (if someone was clearly winning)
Delta <- as.numeric(chess$delta == 1)
Xi <- as.numeric((chess$delta == 0) & (chess$termination == "won by resignation"))

chess$delta <- Delta
chess$xi <- Xi

# Reorder so xi is the third column
chess <- chess[, c(1, 2, ncol(chess), 3:(ncol(chess) - 1))]

chess$intercept <- rep(1, nrow(chess))

# The variable indicating which opening was played is categorical and has a lot
# of levels, most of which only occur once. Therefore, we reclassify them in
# some way

# There are 320 unique openings in this data
length(unique(chess$opening))

# Create data frame of all openings with their respective count
opening_count_df <- data.frame("opening" = unique(chess$opening))
opening_count_df$count <- 0

opening_count_table <- table(chess$opening)
for (i in 1:nrow(opening_count_df)) {
  table_idx <- which(names(opening_count_table) == opening_count_df[i, "opening"])
  opening_count_df[i, "count"] <- opening_count_table[table_idx]
}

opening_count_df <- opening_count_df[order(opening_count_df[, "count"]), ]
row.names(opening_count_df) <- NULL # Reset indices
opening_count_df$cprop <- cumsum(opening_count_df[, "count"])/nrow(chess)

# index of the .33 quantile
q33_idx <- min(which(opening_count_df$cprop > 0.33))
# It can be noticed that this index is very close to the index of the last
# opening which only occurs once. Therefore we will take all openings which
# only occur once as a first category (option 2).

# index of the .66 quantile
q66_idx <- min(which(opening_count_df$cprop > 0.66))
# It can be noticed that this index is very close to the index of the first
# opening which occurs more than 10 times. Therefore we will take all openings
# with an occurence between 2 and 10 to be the second category (option 2).

# Pick a way to categorize openings

# Option 1
chess$opening_cat <- ifelse(chess$opening=="1. e4 e5 2. Nf3 Nf6 3. Nxe5 Nc6 " |
                              grepl("1. e4 c6 2. d4 d5", chess$opening, fixed=TRUE)==TRUE,1,
                            ifelse(grepl("1. e4", chess$opening, fixed=TRUE)==TRUE,2,3))

# Option 2
chess$opening_cat <- rep(0, nrow(chess))
for (i in 1:nrow(chess)) {
  idx_in_opening_count_df <- which(opening_count_df["opening"] == chess[i, "opening"])
  if (opening_count_df[idx_in_opening_count_df, "count"] == 1) {
    category <- 1
  } else if (opening_count_df[idx_in_opening_count_df, "count"] > 10) {
    category <- 3
  } else {
    category <- 2
  }
  chess[i, "opening_cat"] <- category
}

# As can be seen from the histogram, this leaves us with very similar counts
# for each category
hist(chess$opening_cat)


# Just to try out if the variance drops drastically, duplicate the data set 6
# times.

# chess <- rbind(chess, chess, chess, chess, chess, chess)
# Y <- as.matrix(log(chess$time))
# Delta <- as.matrix(chess$delta)

#################### Multiple EV: instrument = opponent_exp ####################

source("Functions_multipleEV.r")

namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}",
               "beta_{T,3}",
               "alpha_{T,1}","alpha_{T,2}",
               "lambda_{T_1}","lambda_{T_2}","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "beta_{C,3}",
               "alpha_{C,1}","alpha_{C,2}","lambda_{C_1}","lambda_{C_2}","sigma_T",
               "sigma_C","rho","theta_1","theta_2")

colnames(chess)
# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, eval_depth, poi_is_white].

# NOTE: as soon as data set is fixed, we should also include ELO of white and
#       black player.
X = as.matrix(subset(chess, select = c("intercept", "white_elo", "black_elo",
                                       "poi_is_white")))

# Define some useful variables
parl <- ncol(X) + 4 # The ones defined above + 2 dummy endogenous + 2 control function
totparl <- 2*parl
parlgamma <- ncol(X)+1 # The ones defined above + instrumental var

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
data=as.matrix(cbind(Y,Delta,Xi,X,Z,W))


init.value.theta_1 <- 1
init.value.theta_2 <- 2

#
# NOTE: check code when estimating gamma: why not just use glm?
#

source("Functions_multipleEV.r")

# Takes about 1 minute to run
parhat_full <- DataApplicationChess.multipleEV(data, init.value.theta_1, init.value.theta_2,
                                               only_2step_model = TRUE)

#################### Multiple EV: instrument = opponent_exp + poi_is_white #####

namescoef.multipleIV <- namescoef[-c(4, 12)]

X.multipleIV = as.matrix(subset(chess, select = c("intercept", "white_elo", "black_elo")))

# Define some useful variables
parl.multipleIV <- ncol(X.multipleIV) + 4 # The ones defined above + 2 dummy endogenous + 2 control function
totparl.multipleIV <- 2*parl.multipleIV
parlgamma.multipleIV <- ncol(X.multipleIV) + 2 # The ones defined above + 2 instrumental vars

# Create data matrix of confounded variable. Z is the indicator of which opening
# was played in the game. Z is a multi-level categorical variable. This means
# that we need more than one dummy variable to represent it. 
Z.multipleIV = as.matrix(chess$opening_cat)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control or treatment group (0 and 1, respectively).
W.multipleIV = as.matrix(cbind(chess$opponent_exp, chess$poi_is_white))

# Create data matrix
XandW.multipleIV = as.matrix(cbind(X.multipleIV, W.multipleIV))
data.multipleIV = as.matrix(cbind(Y, Delta, Xi, X.multipleIV, Z.multipleIV,
                                  W.multipleIV))


init.value.theta_1 <- 1
init.value.theta_2 <- 2

#
# NOTE: check code when estimating gamma: why not just use glm?
#

source("Functions_multipleEV.r")

# Takes about 1 minute to run
parhat_full2 <- DataApplicationChess.multipleEV.multipleIV(data.multipleIV,
                                                           init.value.theta_1,
                                                           init.value.theta_2,
                                                           parl.multipleIV,
                                                           totparl.multipleIV,
                                                           parlgamma.multipleIV,
                                                           namescoef.multipleIV)

#################### Stratification#############################################

source(paste0(parent_folder, "/Functions_ad.r"))

namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}",
               "beta_{T,3}",
               "alpha_{T,1}",
               "lambda_{T_1}","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "beta_{C,3}",
               "alpha_{C,1}","lambda_{C_1}","sigma_T",
               "sigma_C","rho","theta_1","theta_2")



# Define some useful variables
parl <- ncol(X) + 2 # The variables defined above + 1 endogenous var + 1 control function
totparl <- 2*parl
parlgamma <- ncol(X)+1 # The variables defined above + 1 instrumental var


chess$Z1 = ifelse(chess$opening_cat==2,1,0)
chess$Z2 = ifelse(chess$opening_cat==3,1,0)
chess.1 <- chess[(chess$Z1==1 | (chess$Z1==0 & chess$Z2==0)),]
chess.2 <- chess[(chess$Z2==1 | (chess$Z1==0 & chess$Z2==0)),]


# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, eval_depth, poi_is_white].

X.1 = as.matrix(subset(chess.1, select = c("intercept", "white_elo", "black_elo",
                                           "poi_is_white")))
X.2 = as.matrix(subset(chess.2, select = c("intercept", "white_elo", "black_elo",
                                           "poi_is_white")))

# Create data matrix of confounded variable. Z is the indicator of which opening
# was played in the game. Z is a multi-level categorical variable. This means
# that we need more than one dummy variable to represent it. 
Z.1 = as.matrix(chess.1$Z1)
Z.2 = as.matrix(chess.2$Z2)

W.1 = as.matrix(chess.1$opponent_exp)
W.2 = as.matrix(chess.2$opponent_exp)


# observed variables
Y.1 <- as.matrix(log(chess.1$time))
Y.2 <- as.matrix(log(chess.2$time))
Delta.1 <- as.matrix(chess.1$delta)
Delta.2 <- as.matrix(chess.2$delta)
Xi.1 <- as.matrix(chess.1$xi)
Xi.2 <- as.matrix(chess.2$xi)

# Create data matrix
XandW.1 = as.matrix(cbind(X.1,W.1))
n.1=nrow(chess.1)
data.1=as.matrix(cbind(Y.1,Delta.1,Xi.1,X.1,Z.1,W.1))

XandW.2 = as.matrix(cbind(X.2,W.2))
n.2=nrow(chess.2)
data.2=as.matrix(cbind(Y.2,Delta.2,Xi.2,X.2,Z.2,W.2))

init.value.theta_1 <- 1
init.value.theta_2 <- 2

#
# NOTE: check code when estimating gamma: why not just use glm?
#

parhat1=DataApplicationChess(data.1, init.value.theta_1, init.value.theta_2)
parhat2=DataApplicationChess(data.2, init.value.theta_1, init.value.theta_2)


############# Plotting the curves ##############################################

S1_full <- NULL
S2_full <- NULL
par <- parhat_full[[1]]
gamma_1 <-parhat_full[[2]]
gamma_2 <- parhat_full[[3]]
theta <- par[length(par)-1]
s <- par[length(par)-4]

# parameter vector:
# c(intercept, white_elo,black_elo,poi_is_white)
# Suppose Z=2 and Z=3
# W = median(W)
Zobs.2 <- c(1,0)
Zobs.3 <- c(0,1)
dd <- c(1, median(chess$white_elo), median(chess$black_elo),1,median(chess$opponent_exp))

# Calculating control functions
a = dd%*%(gamma_1-gamma_2)
b = dd%*%gamma_1
c = dd%*%gamma_2

E1.nu1 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))/(exp(dd%*%gamma_1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
E2.nu1 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))/(exp(dd%*%gamma_2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
E3.nu1 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))

E1.nu2 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))/(exp(dd%*%gamma_1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
E2.nu2 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))/(exp(dd%*%gamma_2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
E3.nu2 = (1+exp(dd%*%gamma_1)+exp(dd%*%gamma_2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))


V11 = Zobs.2[1]*E1.nu1+Zobs.2[2]*E2.nu1+(1-Zobs.2[1]-Zobs.2[2])*E3.nu1
V12 = Zobs.3[1]*E1.nu1+Zobs.3[2]*E2.nu1+(1-Zobs.3[1]-Zobs.3[2])*E3.nu1
V21 = Zobs.2[1]*E1.nu2+Zobs.2[2]*E2.nu2+(1-Zobs.2[1]-Zobs.2[2])*E3.nu2
V22 = Zobs.3[1]*E1.nu2+Zobs.3[2]*E2.nu2+(1-Zobs.3[1]-Zobs.3[2])*E3.nu2

# Final vectors of covariates
dd.1  <- c(dd[-length(dd)],Zobs.2,V11,V21) #X,Z,V
dd.2  <- c(dd[-length(dd)],Zobs.3,V12,V22) #X,Z,V

Time <- 1:10000


for (i in 1:length(Time)) {
  sd1 = (YJtrans(log(Time[i]), theta) - t(par[1:length(dd.1)]) %*% dd.1)/s
  S1_full[i] = 1 - pnorm(sd1)
  sd2 = (YJtrans(log(Time[i]), theta) - t(par[1:length(dd.2)]) %*% dd.2)/s
  S2_full[i] = 1 - pnorm(sd2)
}

print(paste("Multiple EV: ",exp(IYJtrans(qnorm(0.5)*s+t(par[1:length(dd.1)]) %*% dd.1,theta))))
print(paste("Multiple EV: ",exp(IYJtrans(qnorm(0.5)*s+t(par[1:length(dd.1)]) %*% dd.2,theta))))



S1 <- NULL
S2 <- NULL
theta1 <- parhat1[2*ncol(X) + 8]
theta2 <- parhat2[2*ncol(X) + 8]
s1 <- parhat1[2*ncol(X) + 5]
s2 <- parhat2[2*ncol(X) + 5]
gamma1 <- parhat1[(2*ncol(X)+10):length(parhat1)]
gamma2 <- parhat2[(2*ncol(X)+10):length(parhat2)]


# parameter vector:
# c(intercept, white_elo,black_elo,poi_is_white)
# Suppose Z=2 and Z=3
# W = median(W)
Zobs <- 1
dd <- c(1, median(chess$white_elo), median(chess$black_elo),1,median(chess$opponent_exp))
V1 <- (1-Zobs)*((1+exp(dd%*%gamma1))*log(1+exp(dd%*%gamma1))-(dd%*%gamma1)*exp(dd%*%gamma1))-Zobs*((1+exp(-(dd%*%gamma1)))*log(1+exp(-(dd%*%gamma1)))+(dd%*%gamma1)*exp(-(dd%*%gamma1)))
V2 <- (1-Zobs)*((1+exp(dd%*%gamma2))*log(1+exp(dd%*%gamma2))-(dd%*%gamma2)*exp(dd%*%gamma2))-Zobs*((1+exp(-(dd%*%gamma2)))*log(1+exp(-(dd%*%gamma2)))+(dd%*%gamma2)*exp(-(dd%*%gamma2)))
dd.1  <- c(dd[-length(dd)],Zobs,V1) #X,Z,V
dd.2  <- c(dd[-length(dd)],Zobs,V2) #X,Z,V

for (i in 1:length(Time)) {
  sd1 = (YJtrans(log(Time[i]), theta1) - t(parhat1[1:length(dd.1)]) %*% dd.1)/s1
  S1[i] = 1 - pnorm(sd1)
  sd2 = (YJtrans(log(Time[i]), theta2) - t(parhat2[1:length(dd.2)]) %*% dd.2)/s2
  S2[i] = 1 - pnorm(sd2)
}

print(paste("Stratified: ",exp(IYJtrans(qnorm(0.5)*s1+t(parhat1[1:length(dd.1)]) %*% dd.1,theta1))))
print(paste("Stratified: ",exp(IYJtrans(qnorm(0.5)*s2+t(parhat2[1:length(dd.1)]) %*% dd.2,theta2))))


plot(Time,S1_full, type = 's', col = 1, ylab="Probability", main="Time until checkmate")
lines(Time,S1, type = 's', col = 1, lty=3)
legend(x = 7000, y = 0.85, c("Multiple EV", "Stratified"),
       col = c(1, 1), lty = c(1,3))


plot(Time,S2_full, type = 's', col = 1, ylab="Probability", main="Time until checkmate")
lines(Time,S2, type = 's', col = 1, lty=3)
legend(x = 7000, y = 0.85, c("Multiple EV", "Stratified"),
       col = c(1, 1), lty = c(1,3))
