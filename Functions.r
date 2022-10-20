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
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(y+1,theta)-1)/theta+(1-sg)*(1-power_transform(-y+1,2-theta))/(2-theta)} 
  return(temp) 
} 

IYJtrans = function(y,theta) # Inverse of Yeo-Johnson transformation 
{ 
  sg = y>=0 
  if (theta==0) {temp =(exp(y)-1)*sg+(1-sg)*(1-power_transform(-2*y+1,0.5))} 
  if (theta==2) {temp = sg*(-1+power_transform(2*y+1,0.5))+(1-exp(-y))*(1-sg)} 
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(abs(theta)*y+1,1/theta)-1)+(1-sg)*(1-power_transform(1-(2-theta)*y,1/(2-theta)))} 
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
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))
                 -(XandW%*%gamma)*exp(XandW%*%gamma))
    -Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))
        +(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
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

# - Doesn't assume endogeneity
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

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
  z3 = (transY-(M%*%eta))/sigma2 # b_C
  z4 = ((1-(rho*sigma1/sigma2))*transY-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5) #  term within Phi for C
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)*DtransY # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

# Needed for Hessian matrix when Z is continuous

# - Does assume endogeneity
#   - Z is continuous
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

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

# - Does assume endogeneity
#   - Z is binary
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

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

# - Doesn't assume endogeneity
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

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

# - Does assume endogeneity
#   - Z is continuous
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

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

# - Does assume endogeneity
#   - Z is binary
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

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

# - Doesn't assume endogeneity
# - Does assume dependent censoring
# - Does not use Yeo-Johnson transformation

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

# - Does assume endogeneity
#   - Z is continuous
# - Does assume dependent censoring
# - Does not use Yeo-Johnson transformation

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

# - Does assume endogeneity
#   - Z is binary
# - Does assume dependent censoring
# - Does not use Yeo-Johnson transformation

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

######################## simulation function ###################################

# Some notes:
#   - How to set initial value for theta?
#   - In Deresa, the CI for theta is just constructed as
#         theta-hat +/- 1.96 se(theta-hat).
#     Shouldn't we impose CI(theta) \subset [0,2] (using, f.e. fisher-z transf)?
#   - I didn't change anything to the code to compute the variances on lines
#     ~ 671 - 730 (kan zijn dat de lijnnummers niet meer exact overeen komen
#     omdat, zoals nu, code/uitleg aan het bijschrijven ben). Not sure if 
#     something should change there when also including theta.
#   - If the code is ran with n = 100 matrices become singular. Needing so many
#     observations for an identified model seems like a downside.
#   - We kunnen de naam van deze functie altijd nog veranderen. Voor 't moment 
#     dacht ik dat het wel handig was om een andere naam te gebruiken dan in de 
#     code van Gilles of Deresa.
################################################################################

SimulationCI11_SaraIlias = function(n, nsim, iseed, init.value.theta) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,1)
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,3:(1+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,3:(2+parl)]
    
    # True value for V
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta = init.value.theta
    init = c(rep(0,totparl), 1, 1, init.value.theta)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0),ub=c(rep(Inf,totparl),Inf,Inf, 2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta to initE and replace the original theta (now second-to-last
    # element) with the initial value for rho.
    initE = c(initE,initE[length(initE)])
    initE[length(initE) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl-1]
    t_s2 = 1/parhatE[totparl]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    # Confidence interval for theta
    
    use.fisher.z <- FALSE
    
    if (use.fisher.z) {
      z1theta <- 0.5*log((1+parhatE[length(parhatE)])/(1-parhatE[length(parhatE)]))
      se1theta_z <- (1/(1-parhatE[length(parhatE)]^2))*se1[length(parhatE)]
      z1theta_l <- z1theta-1.96*(se1theta_z)
      z1theta_u <- z1theta+1.96*(se1theta_z)
      
      # Back transform
      
      r1theta_l <- (exp(2*z1theta_l)-1)/(exp(2*z1theta_l)+1)      
      r1theta_u <- (exp(2*z1theta_u)-1)/(exp(2*z1theta_u)+1)
      
    } else {
      r1theta_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
      r1theta_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
    }
    
    # Matrix of all the confidence intervals
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l, r1theta_l),ncol=1),
                matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u, r1theta_u),ncol=1)) 
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta = parhat1[11]
    
    initd <- c(parhat1, parhat1[length(parhat1)])
    initd[length(initd) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1,0),ub=c(rep(Inf,totparl),Inf,Inf,1,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta = parhat[12]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
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
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
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
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      ztheta <- 0.5*log((1+parhat[length(parhat)])/(1-parhat[length(parhat)]))
      setheta_z <- (1/(1-parhat[length(parhat)]^2))*se[length(parhat)]
      ztheta_l <- ztheta-1.96*(setheta_z)
      ztheta_u <- ztheta+1.96*(setheta_z)
      
      # Back transform
      
      rtheta_l <- (exp(2*ztheta_l)-1)/(exp(2*ztheta_l)+1)      
      rtheta_u <- (exp(2*ztheta_u)-1)/(exp(2*ztheta_u)+1)
      
    } else {
      rtheta_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
      rtheta_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    }
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta = parhat1[11]
    
    initd <- c(parhat1, parhat1[length(parhat1)])
    initd[length(initd) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1,0),ub=c(rep(Inf,totparl),Inf,Inf,1,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      zretheta <- 0.5*log((1+parhatre[length(parhatre)])/(1-parhatre[length(parhatre)]))
      seretheta_z <- (1/(1-parhatre[length(parhatre)]^2))*sere[length(parhatre)]
      zretheta_l <- zretheta-1.96*(seretheta_z)
      zretheta_u <- zretheta+1.96*(seretheta_z)
      
      # Back transform
      
      rretheta_l <- (exp(2*zretheta_l)-1)/(exp(2*zretheta_l)+1)      
      rretheta_u <- (exp(2*zretheta_u)-1)/(exp(2*zretheta_u)+1)
      
    } else {
      rretheta_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
      rretheta_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    }
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta
    # [12:14] : params for (gamma_0, gamma_X, gamma_W)
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      zItheta <- 0.5*log((1+parhat1[length(parhat1)])/(1-parhat1[length(parhat1)]))
      seItheta_z <- (1/(1-parhat1[length(parhat1)]^2))*seI[length(parhat1)]
      zItheta_l <- zItheta-1.96*(seItheta_z)
      zItheta_u <- zItheta+1.96*(seItheta_z)
      
      # Back transform
      
      rItheta_l <- (exp(2*zItheta_l)-1)/(exp(2*zItheta_l)+1)      
      rItheta_u <- (exp(2*zItheta_u)-1)/(exp(2*zItheta_u)+1)
      
    } else {
      rItheta_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
      rItheta_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
    }
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta_l),ncol=1),
                matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  
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
  Bias = apply(results[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+4)
  datacp = results[,(2*(totparl+4)+1):(4*(totparl+4))]
  for(i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+2)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+2)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+2) + 1):(2*(totparl+2))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+2))
  datacp = results1[,(2*(totparl+2)+1):(4*(totparl+2))]
  for(i in 1:(totparl+2)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
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
  Bias = apply(results2[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+4)
  datacp = results[,(2*(totparl+4)+1):(4*(totparl+4))]
  for(i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+3)+1):(2*(totparl+3))],2, mean)
  
  CP = rep(0,totparl+3)
  datacp = results3[,(2*(totparl+3) + 1):(4*(totparl+3))]
  for (i in 1:(totparl+3)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-2)],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

# Some notes:
#   - Same notes as for SimulationCI11_SaraIlias
#   - The only difference between this code and Simulation CI11_SaraIlias is 
#     on line 1170 where data is generated (of course the names of the txt-files
#     storing the results is also different.)
################################################################################

SimulationCI12_SaraIlias = function(n, nsim, iseed, init.value.theta) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,2)
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,3:(1+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,3:(2+parl)]
    
    # True value for V
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta = init.value.theta
    init = c(rep(0,totparl), 1, 1, init.value.theta)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0),ub=c(rep(Inf,totparl),Inf,Inf, 2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta to initE and replace the original theta (now second-to-last
    # element) with the initial value for rho.
    initE = c(initE,initE[length(initE)])
    initE[length(initE) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl-1]
    t_s2 = 1/parhatE[totparl]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    # Confidence interval for theta
    
    use.fisher.z <- TRUE
    
    if (use.fisher.z) {
      z1theta <- 0.5*log((1+parhatE[length(parhatE)])/(1-parhatE[length(parhatE)]))
      se1theta_z <- (1/(1-parhatE[length(parhatE)]^2))*se1[length(parhatE)]
      z1theta_l <- z1theta-1.96*(se1theta_z)
      z1theta_u <- z1theta+1.96*(se1theta_z)
      
      # Back transform
      
      r1theta_l <- (exp(2*z1theta_l)-1)/(exp(2*z1theta_l)+1)      
      r1theta_u <- (exp(2*z1theta_u)-1)/(exp(2*z1theta_u)+1)
      
    } else {
      r1theta_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
      r1theta_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
    }
    
    # Matrix of all the confidence intervals
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l, r1theta_l),ncol=1),
                matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u, r1theta_u),ncol=1)) 
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta = parhat1[11]
    
    initd <- c(parhat1, parhat1[length(parhat1)])
    initd[length(initd) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1,0),ub=c(rep(Inf,totparl),Inf,Inf,1,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    # - beta (4 params) = First 4 params of parhat
    # - eta (4 params) = Next 4 params of parhat
    # - sigma1 = parhat[9]
    # - sigma2 = parhat[10]
    # - rho = parhat[11]
    # - theta = parhat[12]
    # - gamma = (intercept, gamma_X, gamma_W)
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
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
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
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
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      ztheta <- 0.5*log((1+parhat[length(parhat)])/(1-parhat[length(parhat)]))
      setheta_z <- (1/(1-parhat[length(parhat)]^2))*se[length(parhat)]
      ztheta_l <- ztheta-1.96*(setheta_z)
      ztheta_u <- ztheta+1.96*(setheta_z)
      
      # Back transform
      
      rtheta_l <- (exp(2*ztheta_l)-1)/(exp(2*ztheta_l)+1)      
      rtheta_u <- (exp(2*ztheta_u)-1)/(exp(2*ztheta_u)+1)
      
    } else {
      rtheta_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
      rtheta_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    }
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta = parhat1[11]
    
    initd <- c(parhat1, parhat1[length(parhat1)])
    initd[length(initd) - 1] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1,0),ub=c(rep(Inf,totparl),Inf,Inf,1,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      zretheta <- 0.5*log((1+parhatre[length(parhatre)])/(1-parhatre[length(parhatre)]))
      seretheta_z <- (1/(1-parhatre[length(parhatre)]^2))*sere[length(parhatre)]
      zretheta_l <- zretheta-1.96*(seretheta_z)
      zretheta_u <- zretheta+1.96*(seretheta_z)
      
      # Back transform
      
      rretheta_l <- (exp(2*zretheta_l)-1)/(exp(2*zretheta_l)+1)      
      rretheta_u <- (exp(2*zretheta_u)-1)/(exp(2*zretheta_u)+1)
      
    } else {
      rretheta_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
      rretheta_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    }
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta
    # [12:14] : params for (gamma_0, gamma_X, gamma_W)
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI)
    
    # Confidence interval for theta
    
    if (use.fisher.z) {
      zItheta <- 0.5*log((1+parhat1[length(parhat1)])/(1-parhat1[length(parhat1)]))
      seItheta_z <- (1/(1-parhat1[length(parhat1)]^2))*seI[length(parhat1)]
      zItheta_l <- zItheta-1.96*(seItheta_z)
      zItheta_u <- zItheta+1.96*(seItheta_z)
      
      # Back transform
      
      rItheta_l <- (exp(2*zItheta_l)-1)/(exp(2*zItheta_l)+1)      
      rItheta_u <- (exp(2*zItheta_u)-1)/(exp(2*zItheta_u)+1)
      
    } else {
      rItheta_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
      rItheta_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
    }
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta_l),ncol=1),
                matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta_u), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  
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
  Bias = apply(results[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+4)
  datacp = results[,(2*(totparl+4)+1):(4*(totparl+4))]
  for(i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+2)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+2)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+2) + 1):(2*(totparl+2))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+2))
  datacp = results1[,(2*(totparl+2)+1):(4*(totparl+2))]
  for(i in 1:(totparl+2)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
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
  Bias = apply(results2[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+4)
  datacp = results[,(2*(totparl+4)+1):(4*(totparl+4))]
  for(i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+3)+1):(2*(totparl+3))],2, mean)
  
  CP = rep(0,totparl+3)
  datacp = results3[,(2*(totparl+3) + 1):(4*(totparl+3))]
  for (i in 1:(totparl+3)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-2)],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}