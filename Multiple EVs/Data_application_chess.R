# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
wd <- getwd()
parent_folder <- substr(wd, 1, nchar(wd) - 13)

chess <- read.csv(paste0(parent_folder, "/duration_df_new.csv"))
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(dplyr)
library(nnet)

chess <- chess[!(is.na(chess$time)),]

# time
Y <- as.matrix(log(chess$time))

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
# It can be notived that this index is very close to the index of the first
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

#################### Multiple EV ###############################################

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

DataApplicationChess.multipleEV(data, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run



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

DataApplicationChess(data.1, init.value.theta_1, init.value.theta_2)
DataApplicationChess(data.2, init.value.theta_1, init.value.theta_2)



