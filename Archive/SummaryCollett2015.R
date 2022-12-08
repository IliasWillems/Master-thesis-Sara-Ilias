rm(list = ls())

library(survival)

data <- read.csv("clean_dataset_JTPA.csv")
data <- data[(data$children == 1 & data$male == 1), ]
data <- data[, -c(5, 8)]

# Plot of observed survival times against some explanatory variables. 
plot(data[data$delta == 0, 9], data[data$delta == 0, 7] - 0.03,
     xlab = "Days", ylab = "Treatment", ylim = c(-0.1, 1.1), xlim = c(0, 1500))
points(data[data$delta == 1, 9], data[data$delta == 1, 7] + 0.03, col = 2)
legend(x = 1200, y = 0.6, c("Event", "Censored"), col = c(2, 1), pch = 1)

plot(data[data$delta == 0, 9], data[data$delta == 0, 3] - 0.03,
     xlab = "Days", ylab = "hsged", ylim = c(-0.1, 1.1), xlim = c(0, 1500))
points(data[data$delta == 1, 9], data[data$delta == 1, 3] + 0.03, col = 2)
legend(x = 1200, y = 0.6, c("Event", "Censored"), col = c(2, 1), pch = 1)

plot(data[data$delta == 0, 9], data[data$delta == 0, 2],
     xlab = "Days", ylab = "age", xlim = c(0, 1500))
points(data[data$delta == 1, 9], data[data$delta == 1, 2], col = 2)
legend(x = 1200, y = 0.6, c("Event", "Censored"), col = c(2, 1), pch = 1)

# Logistic regression model of censoring indicator dependent on the explanatory
# variables.
fit1 <- glm(delta ~ age + hsged + white + married + treatment, data = data,
            family = binomial(link = "logit"))
summary(fit1)
# We find that the probability of censoring is higher for older people. This may
# indicate that survival time and censoring time are dependent.

transplant
