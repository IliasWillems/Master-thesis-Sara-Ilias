
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

# This now gives
#
#   theoretical.value numeric.value
# 1         0.2243771     0.2263654


# Check first integral
nu3star<-5
nu1star <- 3

# theoretical value (according to calculations in paper)
theoretical.value.3.1 <- (nu3star-nu1star-c)*Gumbel_cdf(nu3star-c)-expint_E1(exp(c-nu3star))

integrand.3.1 <- function (nu2star) {
  (nu2star - nu1star)*Gumbel_pdf(nu2star)
}

int.approx.3.1 <- integrate(integrand.3.1,lower=-Inf, upper=nu3star-c)$value
data.frame("theoretical value" = c(theoretical.value.3.1), "numeric value" = c(int.approx.3.1))
# Okay

# Check second integral


integrand.3.2 <- function(nu3star){
  ((nu3star-nu1star-c)*Gumbel_cdf(nu3star-c)-expint_E1(exp(c-nu3star)))*Gumbel_pdf(nu3star)
}

sub_function <- function(w){
  (exp(-w)-exp(-w*(1+exp(-c))))/w
}

theoretical.value.3.2.1 <- 1/(1+exp(c))*(0.5772156649+log(1+exp(c))-(b-c)*exp(-exp(-b-nu1star)*(1+exp(c)))+expint_E1((1+exp(c))*exp(-b-nu1star))-(nu1star+c))
theoretical.value.3.2.2 <- integrate(sub_function,0,exp(c-b-nu1star))$value+(1-Gumbel_cdf(b+nu1star))*expint_E1(exp(c-b-nu1star))
theoretical.value.3.2 <- theoretical.value.3.2.1-theoretical.value.3.2.2

int.approx.3.2 <- integrate(integrand.3.2,lower=b+nu1star, upper=700)$value
data.frame("theoretical value" = c(theoretical.value.3.2), "numeric value" = c(int.approx.3.2))
# Okay




# Check third integral

b <- 5
c <- 2

integrand.3.3 <- function(nu1star){
  integrand.3.3.1 <-1/(1+exp(c))*(0.5772156649+log(1+exp(c))-(b-c)*exp(-exp(-b-nu1star)*(1+exp(c)))+expint_E1((1+exp(c))*exp(-b-nu1star))-(nu1star+c))*Gumbel_pdf(nu1star)
  integrand.3.3.2 <- (integrate(sub_function,0,exp(c-b-nu1star))$value+(1-Gumbel_cdf(b+nu1star))*expint_E1(exp(c-b-nu1star)))*Gumbel_pdf(nu1star)
  return(integrand.3.3.1-integrand.3.3.2)
  }

theoretical.value.3.3 <- 1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b))

int.approx.3.3 <- 0

# Shouldn't this just be multiplying with cell.width instead of cell.width^2
# Trying with cell.width instead gives much more similar, though maybe still
# not satisfactory results.
for (nu1star in nu1star.eval.points) {
  int.approx.3.3 <- int.approx.3.3 + integrand.3.3(nu1star)*(cell.width)
}

data.frame("theoretical value" = c(theoretical.value.3.3), "numeric value" = c(int.approx.3.3))

# This now gives:
# 
#   theoretical.value numeric.value
# 1         0.2243771     0.2243771


# Check third integral in pieces
integrand.3.3.1 <- function(nu1star){
  1/(1+exp(c))*(0.5772156649+log(1+exp(c))-c)*Gumbel_pdf(nu1star)
}

theoretical.value.3.3.1 <- (0.5772156649+log(1+exp(c))-c)/(1+exp(c))
int.approx.3.3.1 <- integrate(integrand.3.3.1,lower=-Inf, upper=Inf)$value
data.frame("theoretical value" = c(theoretical.value.3.3.1), "numeric value" = c(int.approx.3.3.1))
# Okay


integrand.3.3.2 <- function(nu1star){
  1/(1+exp(c))*((b-c)*exp(-exp(-b-nu1star)*(1+exp(c))))*Gumbel_pdf(nu1star)
}

theoretical.value.3.3.2 <- (b-c)/(1+exp(c))/(exp(-b)+exp(c-b)+1)
int.approx.3.3.2 <- integrate(integrand.3.3.2,lower=-Inf, upper=Inf)$value
data.frame("theoretical value" = c(theoretical.value.3.3.2), "numeric value" = c(int.approx.3.3.2))
# Okay

integrand.3.3.3 <- function(nu1star){
  1/(1+exp(c))*expint_E1((1+exp(c))*exp(-b-nu1star))*Gumbel_pdf(nu1star)
}

theoretical.value.3.3.3 <- log(exp(b)/(1+exp(c))+1)/(1+exp(c))
int.approx.3.3.3 <- integrate(integrand.3.3.3,lower=-700, upper=700)$value
data.frame("theoretical value" = c(theoretical.value.3.3.3), "numeric value" = c(int.approx.3.3.3))
# Okay

integrand.3.3.4 <- function(nu1star){
  nu1star/(1+exp(c))*Gumbel_pdf(nu1star)
}

theoretical.value.3.3.4 <- 0.5772156649/(1+exp(c))
int.approx.3.3.4 <- integrate(integrand.3.3.4,lower=-Inf, upper=Inf)$value
data.frame("theoretical value" = c(theoretical.value.3.3.4), "numeric value" = c(int.approx.3.3.4))
# Okay

InnerFunc <- function(w){(exp(-w)-exp(-w*(1+exp(-c))))/w}
InnerIntegral = function(nu1star) { sapply(nu1star, 
                                           function(nu1star) { integrate(InnerFunc, 0, exp(c-b-nu1star))$value })*Gumbel_pdf(nu1star) }
int.approx.3.3.5 <- integrate(InnerIntegral,-100,100)$value

theoretical.value.3.3.5 <- log((1+exp(-c)+exp(b-c))/(1+exp(b-c)))

data.frame("theoretical value" = c(theoretical.value.3.3.5), "numeric value" = c(int.approx.3.3.5))
# Okay


InnerFunc2 <- function(w){exp(-w)/w}
InnerIntegral2 = function(nu1star) { sapply(nu1star, 
                                           function(nu1star) { integrate(InnerFunc2, max(exp(c-b-nu1star), 0.000000001), 50, subdivisions = 20000)$value})*Gumbel_pdf(nu1star)*(1-Gumbel_cdf(b+nu1star)) }
int.approx.3.3.6 <- integrate(InnerIntegral2,0,500)$value


theoretical.value.3.3.6 <- log(1+exp(b-c))-log(1+exp(b-c)+exp(-c))/(1+exp(-b))

data.frame("theoretical value" = c(theoretical.value.3.3.6), "numeric value" = c(int.approx.3.3.6))
# Still different values

# Second step page 65
InnerFunc3 <- function(nu1star){exp(-nu1star-exp(-nu1star))-exp(-nu1star-exp(-nu1star)*(1+exp(-b)))}
InnerIntegral3 = function(w) { sapply(w, 
                                            function(w) { integrate(InnerFunc3, c-b-log(w), 10000000)$value })*exp(-w)/w }
integrate(InnerIntegral3,0,Inf)$value
# Okay (stil zero)

# Third step page 65
func1 <- function(w){(1-exp(-exp(b-c)*w))*exp(-w)/w}
InnerFunc4 <- function(nu1star){1/(1+exp(-b))*exp(-(nu1star-log(1+exp(-b))-exp(-(nu1star-log(1+exp(-b))))))}
InnerIntegral4 = function(w) { sapply(w, 
                                      function(w) { integrate(InnerFunc4, c-b-log(w), 10000000)$value })*exp(-w)/w }

integrate(func1,0,Inf)$value-integrate(InnerIntegral4,0,Inf)$value
# Not zero anymore --> Mistake in this step?

