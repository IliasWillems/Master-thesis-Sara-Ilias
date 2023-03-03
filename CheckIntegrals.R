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


# Check if inner integral is computed correctly. It is
nu1star <- -5

subintegrand <- function (nu2star) {
  (nu2star - nu1star)*Gumbel_pdf(nu2star)*Gumbel_cdf(c + nu2star)
}

gamma_EM <- -digamma(1)
theoretical.subvalue <- (1/(1+exp(-c)))*(gamma_EM + log(1+exp(-c)) - a*exp(-exp(-a-nu1star)*(1+exp(-c))) + expint_E1((1+exp(-c))*exp(-a-nu1star)) - nu1star)

subint.approx <- integrate(subintegrand, a+nu1star, Inf)$value
data.frame("theoretical value" = c(theoretical.subvalue), "numeric value" = c(subint.approx))


# Check if each of the terms in the outer integral are integrated correctly

# First term seems correct
integrand.first.term <- function (nu1star) {
  Gumbel_pdf(nu1star)*a*exp(-exp(-a-nu1star)*(1+exp(-c)))
}

integrate(integrand.first.term, -Inf, Inf)
a/(exp(-a)+exp(-c-a)+1)

# Second term seems correct, EXCEPT THERE IS A SIGN DIFFERENCE!
integrand.second.term <- function (nu1star) {
  Gumbel_pdf(nu1star)*expint_E1((1+exp(-c))*exp(-a-nu1star))
}

integrate(integrand.second.term, -100, +100)
-log(1+exp(a)/(1+exp(-c)))

# Second equation
integrand.second.equation <- function (nu1star, u) {
  (exp(u)/u)*Gumbel_pdf(nu1star)
}
second.equation.approx <- 0
for (nu1star in nu1star.eval.points) {
  for (u in nu2star.eval.points) {
    if (u < -(1+exp(-c))*exp(-a-nu1star)) {
      second.equation.approx <- second.equation.approx + integrand.second.equation(nu1star, u)*(cell.width^2)
    }
  }
}
second.equation.approx


# Fifth equation already has the mistake
fifth.equation <- function (u) {
  exp(u)/u - exp((exp(a)/(1+exp(-c)) + 1)*u)/u
}
integrate(fifth.equation, -Inf, 0)


# Third term
integrand.third.term <- function (nu1star) {
  Gumbel_pdf(nu1star)*nu1star
}

integrate(integrand.third.term, -Inf, Inf)
gamma_EM












