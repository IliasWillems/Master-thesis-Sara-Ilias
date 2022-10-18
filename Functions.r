

######################### functions Yeo-Johnson ################################

log_transform = function(y)   
{ 
  transform_y = (y>0)*y+(y<=0)*1 # sy>0 --> condition is a check
  return(log(transform_y)) 
} 

power_transform = function(y,pw) 
{ 
  transform_y = (y>0)*y+(y<=0)*1  # sy>0
  return(transform_y^pw) 
} 

YJtrans = function(y,theta) # Yeo-Johnson transformation 
{ 
  sg = y>=0 
  if (theta==0) {temp = log_transform(y+1)*sg+(1-sg)*(0.5-0.5*(y-1)^2)} 
  if (theta==2) {temp = sg*(-0.5+0.5*(y+1)^2)-log_transform(-y+1)*(1-sg)} 
  if ((theta!=0) & (theta!=2)) {temp = 
    sg*(power_transform(y+1,theta)-1)/theta+(1-sg)*
    (1-power_transform(-y+1,2-theta))/(2-theta)} 
  return(temp) 
} 

IYJtrans = function(y,theta) # Inverse of Yeo-Johnson transformation 
{ 
  sg = y>=0 
  if (theta==0) {temp =(exp(y)-1)*sg+(1-sg)*(1-power_transform(-2*y+1,0.5))} 
  if (theta==2) {temp = sg*(-1+power_transform(2*y+1,0.5))+(1-exp(-y))*(1-sg)} 
  if ((theta!=0) & (theta!=2)) {temp = 
    sg*(power_transform(abs(theta)*y+1,1/theta)-1)
  +(1-sg)*(1-power_transform(1-(2-theta)*y,1/(2-theta)))} 
  return(temp) 
} 

DYJtrans = function(y,theta) # Derivative of Yeo-Johnson transformation 
{ 
  sg = y>=0 
  temp = power_transform(y+1,theta-1)*sg+power_transform(-y+1,1-theta)*(1-sg) 
  return(temp) 
} 


###################### data simulating function ################################

dat.sim.reg = function(n,par,iseed,Zbin,Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bivariate normal distribution of error terms
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  err = mvrnorm(n, mu =mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) # Bernoulli with p =0.5
  {
    W = sample(c(0,1), n, replace = TRUE)} # sample 0 and 1 with equal probability
  else if (Wbin==1)
  {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2)  # nu is standard logistic
  {
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))
                -(XandW%*%gamma)*exp(XandW%*%gamma))
                -Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))
                +(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  }
  else if (Zbin==1) # nu is standard normal
  {
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = Mgen%*%beta+err1 # model YJ transformed time with real covariates
  C = Mgen%*%eta+err2 # model YJ transformed censoring with real covariates
  
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C) # observed YJ transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  Y = IYJtrans(Y,sd[4]) # observed time
  data = cbind(Y,d1,M,realV) # data consisting of observed time,
                            # censoring indicator, all data and the control function
  
  return(data)
}


######################## likelihood function ###################################


# maximum likelihood for Gamma (Z continuous)
LikGamma1 = function(par,Y,M){ 
  W=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (Y-W%*%gamma)^2
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(Logn)
}

# maximum likelihood for Gamma (Z discrete)
LikGamma2 = function(par,Y,M){ 
  W=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (plogis(W%*%gamma)^Y)*((1-plogis(W%*%gamma))^(1-Y))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

###########################################################
# joint model with dependent censoring and transformation #
###########################################################

LikF = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  rho = par[l+3]
  theta = par[l+4]
  
  transY=YJtrans(Y,theta)
  DtransY=DYJtrans(Y,theta)
  
  z1 = (transY-(M%*%beta))/sigma1 # b_T
  z2 = ((1-(rho*sigma2/sigma1))*transY-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5)) #  term within Phi for T
  z3 = (Y-(M%*%eta))/sigma2 # b_C
  z4 = ((1-(rho*sigma1/sigma2))*transY-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5) #  term within Phi for C
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)*DtransY # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

# Needed for Hessian matrix when Z is continuous
LikFG1 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2 # remove Z and W
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1] 
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  theta = par[l+6]
  gamma = as.matrix(par[(l+7):(l+6+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-XandW%*%gamma
  
  transY=YJtrans(Y,theta)
  DtransY=DYJtrans(Y,theta)
  
  # likelihood for estimated parameters theta(beta,eta,alpha, lambda, sigma,rho, theta) and gamma
  z1 = (transY-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*transY-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (transY-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*transY-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)*DtransY 
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn) # returns likelihood (the same as likF but including gamma as parameter)
}

# needed for Hessian matrix when Z is binary
LikFG2 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  theta = par[l+6]
  gamma = as.matrix(par[(l+7):(l+6+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
  transY=YJtrans(Y,theta)
  DtransY=DYJtrans(Y,theta)
  
  z1 = (transY-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*transY-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (transY-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*transY-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)*DtransY
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}


################################
# Independent model assumption #
################################

LikI = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M) 
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k]) # parameters related to T
  eta = as.matrix(par[v:l]) # parameters related to C
  sigma1 = par[l+1] # sigma of T
  sigma2 = par[l+2] # sigma of C
  theta = par[l+3]
  
  transY = YJtrans(Y,theta)
  DtransY = DYJtrans(Y,theta)
  
  z1 = (transY-(M%*%beta))/sigma1 # term within Phi for T (rho = 0)
  z2 = (transY-(M%*%eta))/sigma2 # term withing Phi for C (rho = 0)
  
  # likelihood (when rho = 0)
  # tot gives sample size number of rows
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))*DtransY
  p1 = pmax(tot,1e-100) 
  Logn = sum(log(p1)); # calculate loglikelihood
  return(-Logn)
}

# needed for Hessian matrix

# no dependent censoring (Z continuous)
LikIGamma1 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  theta = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-(XandW%*%gamma) #estimating the V 
  
  transY = YJtrans(Y,theta)
  DtransY = DYJtrans(Y,theta)
  
  z1 = (transY-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = (transY-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))*DtransY
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

# no dependent censoring (Z discrete)
LikIGamma2 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  theta = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
  transY = YJtrans(Y,theta)
  DtransY = DYJtrans(Y,theta)
  
  z1 = (transY-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = (transY-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta)*DtransY
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}


################################################
# No transformation assumption (model Gilles)  #
################################################

LikFNT = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  rho = par[l+3]

  z1 = (Y-(M%*%beta))/sigma1 # b_T
  z2 = ((1-(rho*sigma2/sigma1))*Y-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5)) #  term within Phi for T
  z3 = (Y-(M%*%eta))/sigma2 # b_C
  z4 = ((1-(rho*sigma1/sigma2))*Y-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5) #  term within Phi for C
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta) # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

# Needed for Hessian matrix when Z is continuous
LikFNTG1 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2 # remove Z and W
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1] 
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-XandW%*%gamma

  
  # likelihood for estimated parameters theta(beta,eta,alpha, lambda, sigma,rho, theta) and gamma
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*Y-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*Y-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn) # returns likelihood (the same as likF but including gamma as parameter)
}

# needed for Hessian matrix when Z is binary
LikFNTG2 = function(par,Y,Delta,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*Y-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*Y-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}



