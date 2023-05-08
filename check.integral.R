#
# This script is used to numerically check the solution of the integrals in the
# section "Multiple endogenous variables".
#

library(expint)

rm(list = ls())

# Set parameters to determine the grid over which the integral will be evaluated
grid.lower.bound <- -20
grid.upper.bound <- 20
grid.length <- grid.upper.bound - grid.lower.bound

number.of.grid.cells <- 1000
cell.width <- grid.length/number.of.grid.cells


# Define the points in which the integrand will be evaluated
nu1star.eval.points <- seq(grid.lower.bound + cell.width/2,
                           grid.upper.bound - cell.width/2,
                           length.out = number.of.grid.cells)
nu2star.eval.points <- nu1star.eval.points
nu3star.eval.points <- nu1star.eval.points

# Probability density function of standard Gumbel distribution
Gumbel_pdf <- function (x) {
  exp(-x-exp(-x))
}

# Cumulative density function of standard Gumbel distribution
Gumbel_cdf <- function(x) {
  exp(-exp(-x))
}

################################################################################

# Check for integrals in E[nu1 | Z1 = 1, Z2 = 0], hence note that we leave out
# the reciprocal of the probability P(Z1 = 1, Z2 = 0) in the following.

# Set parameters
a <- 5
b <- -2

# theoretical value (according to calculations in paper)
theoretical.value.1 <- a/(1+exp(-a)+exp(-b)) - (1/(1+exp(-b)))*log(1+exp(a)*(1+exp(-b)))

integrand.1 <- function (nu1star, nu2star, b) {
  (nu2star - nu1star)*Gumbel_pdf(nu1star)*Gumbel_pdf(nu2star)*Gumbel_cdf(b + nu1star)
}

int.approx.1 <- 0

# This takes some time
for (nu1star in nu1star.eval.points) {
  for (nu2star in nu2star.eval.points) {
    if (nu2star < a + nu1star) {
      int.approx.1 <- int.approx.1 + integrand.1(nu1star, nu2star, b)*(cell.width^2)
    }
  }
}

data.frame("theoretical value" = c(theoretical.value.1), "numeric value" = c(int.approx.1))


################################################################################

# Check for integrals in E[nu1 | Z1 = 0, Z2 = 1]

# Set parameters
a <- -2
c <- -5

# theoretical value (according to calculations in paper)
theoretical.value.2 <- (1/(1+exp(-c)))*(-(a/(exp(-a)+exp(-a-c)+1))+log(exp(a)+exp(-c)+1))


integrand.2 <- function (nu1star, nu2star, c) {
  (nu2star - nu1star)*Gumbel_pdf(nu1star)*Gumbel_pdf(nu2star)*Gumbel_cdf(c + nu2star)
}

int.approx.2 <- 0

# This takes some time
for (nu1star in nu1star.eval.points) {
  for (nu2star in nu2star.eval.points) {
    if (nu2star > a + nu1star) {
      int.approx.2 <- int.approx.2 + integrand.2(nu1star, nu2star, c)*(cell.width^2)
    }
  }
}

data.frame("theoretical value" = c(theoretical.value.2), "numeric value" = c(int.approx.2))

################################################################################

# Check for integrals in E[nu1 | Z1 = 0, Z2 = 0]

# Set parameters
b <- -2
c <- -5

# theoretical value (according to calculations in paper)
theoretical.value.3 <- 1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b))

integrand.3 <- function (nu1star, nu2star,nu3star) {
  (nu2star - nu1star)*Gumbel_pdf(nu1star)*Gumbel_pdf(nu2star)*Gumbel_pdf(nu3star)
}

int.approx.3 <- 0

# This takes a long time
for (nu1star in nu1star.eval.points) {
  if (which(nu1star.eval.points == nu1star) %% 100 == 0) {
    message(100*which(nu1star.eval.points == nu1star)/length(nu1star.eval.points), "% completion")
  }
  for (nu2star in nu2star.eval.points) {
    for (nu3star in nu3star.eval.points) {
      if (nu3star > b + nu1star) {
        if(nu2star < nu3star-c)
          int.approx.3 <- int.approx.3 + integrand.3(nu1star, nu2star, nu3star)*(cell.width^3)
      }
    }
  }
}

data.frame("theoretical value" = c(theoretical.value.3), "numeric value" = c(int.approx.3))

# This gives
#
#   theoretical.value numeric.value
# 1         0.2243771     0.2263654
