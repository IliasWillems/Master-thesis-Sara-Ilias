
###################### data generation functions ###############################


data.misspecified.probit = function(n, par, iseed, Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bivariate normal distribution of error terms
  mu = c(0, 0)
  sigma = matrix(c(sd[1]^2, sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2), ncol=2)
  err = mvrnorm(n, mu=mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  # This is the only difference with dat.sim.reg: v is now sampled from a
  # standard normal distribution and the conditional expectations are now based
  # on a normally distributed random variable.
  v <- rnorm(n)
  Z <- as.matrix(as.numeric(XandW %*% gamma - v > 0))
  realV <- (1-Z)*(dnorm(XandW %*% gamma)/pnorm(-(XandW %*% gamma))) - Z*(dnorm(XandW %*% gamma)/pnorm(XandW %*% gamma))
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4])# model YJ transformed time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model YJ transformed censoring with real covariates
  A = runif(n,0,8)
  
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed YJ transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}


data.misspecified.cloglog = function(n, par, iseed, Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bivariate normal distribution of error terms
  mu = c(0, 0)
  sigma = matrix(c(sd[1]^2, sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2), ncol=2)
  err = mvrnorm(n, mu=mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  # This is the only difference with dat.sim.reg: v is now sampled from a
  # standard Gumbel distribution and the conditional expectations are now based
  # on a Gumbel distributed random variable.
  
  # Generate random sample from Gumbel distribution using probability integral
  # transformation.
  v <- -log(-log(runif(n)))
  Z <- as.matrix(as.numeric(XandW %*% gamma - v > 0))
  
  gamma.EM <- -digamma(1) # Euler-Mascheroni constant
  
  # Compute E[v|v < W^T \gamma]
  E.v_lt_Wg <- (XandW %*% gamma) - exp(exp(-XandW %*% gamma))*expint_E1(x = (exp(- XandW %*% gamma)))
  
  # Compute E[v|v > W^T \gamma]
  E.v_gt_Wg <- (gamma.EM - (XandW %*% gamma)*exp(-exp(-XandW %*% gamma)) + expint_E1(x = (exp(-XandW %*% gamma)))) / (1-exp(-exp(-XandW %*% gamma)))
  
  # Compute real V
  realV <- (1-Z)*(E.v_gt_Wg) - Z*(E.v_lt_Wg)
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model YJ transformed time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5])  # model YJ transformed censoring with real covariates
  A = runif(n,0,8)
  
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed YJ transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}


data.misspecified.skew <- function(n, par, iseed, Zbin, Wbin) {
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  psi.delta.theta = par[[3]]
  gamma = par[[4]]
  
  offdiag.correlation <- psi.delta.theta[1]
  delta_1 <- psi.delta.theta[2]
  delta_2 <- psi.delta.theta[3]
  theta_1 <- psi.delta.theta[4]
  theta_2 <- psi.delta.theta[5]
  
  parl = length(par[[1]])
  totparl = 2*parl
  parlgamma = (parl-1)
  
  # Some input validation
  if ((offdiag.correlation <= -1) || (offdiag.correlation >= 1)) {
    stop("Invalid value of the offdiagonal element in the correlation matrix")
  }
  
  if ((delta_1 <= -1) || (delta_1 >= 1)) {
    stop("All deltas must be between -1 and 1")
  }
  
  if ((delta_2 <= -1) || (delta_2 >= 1)) {
    stop("All deltas must be between -1 and 1")
  }
  
  #
  # Creation of errors and computation of some relevant parameters
  #
  
  # Draw (Y_1, Y_2, Y_3) from a normal distribution with mean zero and specified
  # variance-covariance matrix.
  mu <- c(0, 0, 0)
  Psi <- matrix(c(1, offdiag.correlation, offdiag.correlation, 1), ncol = 2) 
  sigma <- matrix(c(1, 0, 0, 0, Psi[1], Psi[2], 0, Psi[3], Psi[4]), nrow = 3)
  Y <- mvrnorm(n, mu=mu , Sigma=sigma)
  
  # Compute the value for omega.
  lambda <- c(delta_1/sqrt(1-delta_1^2), delta_2/sqrt(1-delta_2^2))
  Delta <- diag(c(sqrt(1-delta_1^2), sqrt(1-delta_2^2)), nrow = 2, ncol = 2)
  Omega <- Delta %*% (Psi + (lambda %*% t(lambda))) %*% Delta
  offdiag.Omega <- Omega[2]
  
  # Create the skew-normal distributed errors (we use the variable Z for these
  # errors to be consistent with the notation used in the paper of Azzalini).
  Z_1 <- delta_1*abs(Y[,1]) + sqrt(1 - delta_1^2)*Y[,2]-(2/pi)^(1/2)*delta_1
  Z_2 <- delta_2*abs(Y[,1]) + sqrt(1 - delta_2^2)*Y[,3]-(2/pi)^(1/2)*delta_2

  # Theoretical correlation
  numerator <- offdiag.correlation*sqrt((1-delta_1^2)*(1-delta_2^2)) + delta_1*delta_2*(1-2*(1/pi))
  denominator <- sqrt((1-2*(1/pi)*delta_1^2)*(1-2*(1/pi)*delta_2^2))
  theoretical.correlation <- numerator/denominator
  
  #
  # Creation of event and censoring times
  #
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T =IYJtrans(Mgen%*%beta+Z_1,theta_1) # model YJ transformed time with real covariates
  C = IYJtrans(Mgen%*%eta+Z_2,theta_2) # model YJ transformed censoring with real covariates
  A = runif(n,0,8)
  
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed YJ transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  errors <- cbind(Z_1, Z_2)
  
  return(list(data, errors, offdiag.Omega, theoretical.correlation, lambda))
}


data.misspecified.t = function(n,par,iseed,Zbin,Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bivariate t distribution of error terms
  # scale matrix S such that S*df/(df-2) is variance-covariance matrix if df>2
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  err = rmvt(n, mu =mu , S=sigma*(sd[6]-2)/sd[6], sd[6])
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model censoring with real covariates
  A = runif(n,0,8)
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed non-transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}


data.misspecified.bimodal = function(n,par,iseed,Zbin,Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bimodal normal distribution
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  errs1 <- mvrnorm(n/2, mu=c(-1, 1), Sigma=sigma)
  errs2 <- mvrnorm(n/2, mu=c(1, -1), Sigma=sigma)
  err <- rbind(errs1, errs2)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model censoring with real covariates
  A = runif(n,0,8)
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed non-transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}


data.misspecified.heteroscedastic = function(n,par,iseed,Zbin,Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
 # simulate data
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  
  # heteroscedastic normal distribution
  # variance increases with T (and C)
  # higher value of Mgen%*%beta (or Mgen%*%eta) results in higher variance
  mu = c(0,0)

  err=matrix(rep(1,2*n), nrow=n)
  for (A in 1:n){
    # take exponential such that product always positive
    # divide by 10 to fade out differences
    prod_t <- exp((Mgen%*%beta)[A]/5)
    prod_c <- exp((Mgen%*%eta)[A]/5)
    err[A,] <- mvrnorm(1, mu=mu, Sigma=matrix(c(sd[1]^2*(1+prod_t),sd[1]*sd[2]*sd[3]*sqrt((1+prod_t)*(1+prod_c)),sd[1]*sd[2]*sd[3]*sqrt((1+prod_t)*(1+prod_c)),sd[2]^2*(1+prod_c)),nrow=2,byrow=TRUE))
  }
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model censoring with real covariates
  A = runif(n,0,8)
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed non-transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}

############################# Simulation functions #############################

Misspecificaion.probit.sim = function(n, nsim, iseed, par, Wbin) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = data.misspecified.probit(n, par, iseed+i, Wbin)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2 = 1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta = parhat[12]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    psii = -WMI%*%t(mi)
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(par[[1]],par[[2]],par[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  sum = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  sum <- sum[,-3]
  
  # Make nice LaTeX table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_probit",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
}


Misspecificaion.cloglog.sim = function(n, nsim, iseed, par, Wbin) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = data.misspecified.cloglog(n, par, iseed+i, Wbin)
    
    Y = data[,1]
    Delta = data[,2]
    Xi=data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2=1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta = parhat[12]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    psii = -WMI%*%t(mi)
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))  
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(par[[1]],par[[2]],par[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  sum = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  sum <- sum[,-3]
  
  # Make nice LaTeX table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_cloglog",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
}


Misspecification.skew.sim = function(n, nsim, iseed, par, Zbin, Wbin) {
  
  # Some input validation
  if ((Zbin != 1) || (Wbin != 1)) {
    stop("Not implemented yet")
  }
  
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {

    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    # Generate data with skew-normal error distributions
    output <- data.misspecified.skew(n, par, iseed + i, Zbin, Wbin)
    data <- output[[1]]
    theoretical.correlation <- output[[4]]
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2 = 1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0

    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta_1 = parhat[12]
    # - theta_2 = parhat[13]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector. Note that we did not
  # specify variances, nor did we specify directly the correlation between the 
  # transformed log-event and log-censoring times. Therefore, we temporarily put
  # a "-1" for the variances and use the theoretical correlation.
  
  par0 = c(par[[1]], par[[2]], c(-1, -1, theoretical.correlation, par[[3]][4], par[[3]][5]))
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(summary) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(summary) = namescoef
  
  # Make nice Latex table
  xtab = xtable(summary[,-3])
  print(xtab, sanitize.text.function = function(x){x})
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_skew",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
  
  return(summary)
}



Misspecification.t.sim = function(n, nsim, iseed, par, Zbin, Wbin) {
  
  # Some input validation
  if ((Zbin != 1) || (Wbin != 1)) {
    stop("Not implemented yet")
  }
  
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    # Generate data with skew-normal error distributions
    data <- data.misspecified.t(n, par, iseed + i, Zbin, Wbin)

    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2 = 1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta_1 = parhat[12]
    # - theta_2 = parhat[13]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector. Note that we did not 
  # specify the degrees of freedom. We did specify variances but it will not
  # make sense to compare them (since distribution error term misspecified)
  
  par0 = c(par[[1]], par[[2]], par[[3]][1:length(par[[3]])-1])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(summary) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(summary) = namescoef
  
  # Make nice Latex table
  xtab = xtable(summary[,-3])
  print(xtab, sanitize.text.function = function(x){x})
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_t",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
  
  return(summary)
}

Misspecification.bimodal.sim = function(n, nsim, iseed, par, Zbin, Wbin) {
  
  # Some input validation
  if ((Zbin != 1) || (Wbin != 1)) {
    stop("Not implemented yet")
  }
  
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    # Generate data with skew-normal error distributions
    data <- data.misspecified.bimodal(n, par, iseed + i, Zbin, Wbin)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2 = 1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta_1 = parhat[12]
    # - theta_2 = parhat[13]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector. We did specify
  # variances but it will not make sense to compare them 
  # (since distribution error term misspecified)
  
  par0 = c(par[[1]], par[[2]], par[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(summary) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(summary) = namescoef
  
  # Make nice Latex table
  xtab = xtable(summary[,-3])
  print(xtab, sanitize.text.function = function(x){x})
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_bimodal",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
  
  return(summary)
}

Misspecification.heteroscedastic.sim = function(n, nsim, iseed, par, Zbin, Wbin) {
  
  # Some input validation
  if ((Zbin != 1) || (Wbin != 1)) {
    stop("Not implemented yet")
  }
  
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    # Generate data with skew-normal error distributions
    data <- data.misspecified.heteroscedastic(n, par, iseed + i, Zbin, Wbin)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = 1
    # - theta_2 = 1
    init = c(rep(0,totparl), 1, 1, 1, 1)
    
    # Independent model for starting values sigma and theta.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta_1 = parhat[12]
    # - theta_2 = parhat[13]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector. We did specify
  # variances but it will not make sense to compare them 
  # (since error term is heteroscedastic)
  
  par0 = c(par[[1]], par[[2]], par[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V
  
  colnames(summary) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(summary) = namescoef
  
  # Make nice Latex table
  xtab = xtable(summary[,-3])
  print(xtab, sanitize.text.function = function(x){x})
  
  # set to 3 significant digits
  digits(xtab) = rep(3,5)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, file = paste0("YJ_misspec_heteroscedastic",n,".txt"),
               add.to.row = addtorow, append = FALSE, table.placement="!",
               sanitize.text.function = function(x){x})
  print(xtab, include.colnames = TRUE, add.to.row = addtorow, 
        table.placement="!", sanitize.text.function = function(x){x})
  
  return(summary)
}