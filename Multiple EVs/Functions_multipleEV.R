# The data simulation function is changed so that it can (only) generate data
# where Z is a multi-level categorical variate.



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

dat.sim.reg.Zbin = function(n,par,iseed){
  
  parl = length(par[[1]])
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma1 = par[[4]]
  gamma2 = par[[5]]
  
  # bivariate normal distribution of error terms
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  err = mvrnorm(n, mu =mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  

 # Bernoulli with p =1/3
  W = sample(c(0,1,2), n, replace = TRUE) # sample 0 and 1,2 with equal probability
  W1 = ifelse(W==1,1,0)
  W2 = ifelse(W==2,1,0)
  XandW=as.matrix(cbind(x0,x1,W1,W2))

  
  nu1star <- rgumbel(n)
  nu2star <- rgumbel(n)
  nu3star <- rgumbel(n)
  
  nu1 <- nu2star-nu1star
  nu2 <- nu3star-nu1star
  
  #gamma1 = gamma1star-gamma3star
  #gamma2 = gamma2star-gamma3star
  
  a = XandW%*%(gamma1-gamma2)
  b = XandW%*%gamma1
  c = XandW%*%gamma2

  Z1 = as.matrix(as.numeric(a>nu1)*as.numeric(b>nu2))
  Z2 = as.matrix(as.numeric(a<nu1)*as.numeric(c>(nu2-nu1)))
  E1.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
  E2.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
  E3.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))
  
  E1.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
  E2.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
  E3.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))
  
  
  realV1 = Z1*E1.nu1+Z2*E2.nu1+(1-Z1-Z2)*E3.nu1
  realV2 = Z1*E1.nu2+Z2*E2.nu2+(1-Z1-Z2)*E3.nu2

  
  Mgen = matrix(c(x0,x1,Z1,Z2,realV1,realV2),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model censoring with real covariates
  A = runif(n,0,15)
  M = matrix(c(x0,x1,Z1,Z2,W1,W2),ncol=parl,nrow=n)    # data matrix


  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed non-transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV1,realV2) # data consisting of observed time,
                            # censoring indicator, all data and the control function
  
  return(data)
}


######################## likelihood function ###################################


# maximum likelihood for Gamma (Z discrete)
LikGamma2 = function(par,Y,M){ 
  W=as.matrix(M)
  gamma1= as.matrix(par[1:parlgamma])
  gamma2= as.matrix(par[(parlgamma+1):(2*parlgamma)])
  Y1 = Y[,1]
  Y2 = Y[,2]
  
  tot = (exp(W%*%gamma1)/(1+exp(W%*%gamma1)+exp(W%*%gamma2)))^Y1*(exp(W%*%gamma2)/(1+exp(W%*%gamma1)+exp(W%*%gamma2)))^Y2*(1/(1+exp(W%*%gamma1)+exp(W%*%gamma2)))^(1-Y1-Y2)
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

###########################################################
# joint model with dependent censoring and transformation #
###########################################################

# - Used in second step of 2-step estimation
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

LikF = function(par,Y,Delta,Xi,M){ 
  M = as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  rho = par[l+3]
  theta_1 = par[l+4]
  theta_2 = par[l+5]
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(M%*%beta))/sigma1 # b_T
  z2 = ((transY.C-(rho*sigma2/sigma1)*transY.T)-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5)) #  term within Phi for T
  z3 = (transY.C-(M%*%eta))/sigma2 # b_C
  z4 = ((transY.T-(rho*sigma1/sigma2)*transY.C)-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5) #  term within Phi for C
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4))*DtransY.C)^Xi*(pbinorm(q1=-z1,q2=-z3,cov12=rho))^(1-(Delta+Xi)) # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}


# needed for Hessian matrix when Z is binary

# - Does assume endogeneity
#   - Z is binary
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

LikFG2 = function(par,Y,Delta,Xi,M,discrete){ 
  M = as.matrix(M)
  if (discrete==1){
    k = ncol(M)-4
    W=as.matrix(M[,(k+3):(k+4)])
  }else{
    k=ncol(M)-3
    W=as.matrix(M[,(k+3)])
  }
  l = 2*(k+2)
  v = k+5
  
  parlgamma <- (length(par) - (l+9))/2

  beta = as.matrix(par[1:k])
  alphaT.1 = par[k+1]
  alphaT.2 = par[k+2]
  lambdaT.1 = par[k+3]
  lambdaT.2 = par[k+4]
  eta = as.matrix(par[v:l])
  alphaC.1 = par[l+1]
  alphaC.2 = par[l+2]
  lambdaC.1 = par[l+3]
  lambdaC.2 = par[l+4]
  sigma1 = par[l+5]
  sigma2 = par[l+6]
  rho = par[l+7]
  theta_1 = par[l+8]
  theta_2 = par[l+9]
  gamma1 = as.matrix(par[(l+10):(l+9+parlgamma)])
  gamma2 = as.matrix(par[(l+10+parlgamma):(l+9+2*parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z1=as.matrix(M[,(k+1)])
  Z2=as.matrix(M[,(k+2)])
  XandW=as.matrix(cbind(X,W))
  
  a = XandW%*%(gamma1-gamma2)
  b = XandW%*%gamma1
  c = XandW%*%gamma2
  
  E1.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
  E2.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
  E3.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))
  
  E1.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
  E2.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
  E3.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))
  
  
  Vest1 = Z1*E1.nu1+Z2*E2.nu1+(1-Z1-Z2)*E3.nu1
  Vest2 = Z1*E1.nu2+Z2*E2.nu2+(1-Z1-Z2)*E3.nu2
  
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(X%*%beta+Z1*alphaT.1+Z2*alphaT.2+Vest1*lambdaT.1+Vest2*lambdaT.2))/sigma1
  z2 = ((transY.C-rho*sigma2/sigma1*transY.T)-((X%*%eta+Z1*alphaC.1+Z2*alphaC.2+Vest1*lambdaC.1+Vest2*lambdaC.2)-rho*(sigma2/sigma1)*(X%*%beta+Z1*alphaT.1+Z2*alphaT.2+Vest1*lambdaT.1+Vest2*lambdaT.2)))/(sigma2*(1-rho^2)^0.5)
  z3 = (transY.C-(X%*%eta+Z1*alphaC.1+Z2*alphaC.2+Vest1*lambdaC.1+Vest2*lambdaC.2))/sigma2
  z4 = ((transY.T-rho*sigma1/sigma2*transY.C)-((X%*%beta+Z1*alphaT.1+Z2*alphaT.2+Vest1*lambdaT.1+Vest2*lambdaT.2)-rho*(sigma1/sigma2)*(X%*%eta+Z1*alphaC.1+Z2*alphaC.2+Vest1*lambdaC.1+Vest2*lambdaC.2)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z3)*(1-pnorm(z4))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z3,cov12=rho))^(1-(Delta+Xi))) # likelihood
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

LikI = function(par,Y,Delta,Xi,M){ 
  M = as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k]) # parameters related to T
  eta = as.matrix(par[v:l]) # parameters related to C
  sigma1 = par[l+1] # sigma of T
  sigma2 = par[l+2] # sigma of C
  theta_1 = par[l+3]
  theta_2 = par[l+4]
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(M%*%beta))/sigma1 # term within Phi for T (rho = 0)
  z2 = (transY.C-(M%*%eta))/sigma2 # term withing Phi for C (rho = 0)
  
  # likelihood (when rho = 0)
  # tot gives sample size number of rows
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z2,cov12=0))^(1-(Delta+Xi)))
  p1 = pmax(tot,1e-100) 
  Logn = sum(log(p1)); # calculate loglikelihood
  return(-Logn)
}

# needed for Hessian matrix

# no dependent censoring (Z discrete)

# - Does assume endogeneity
#   - Z is binary
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

LikIGamma2 = function(par,Y,Delta,Xi,M,discrete){ 
  M = as.matrix(M)
  if (discrete==1){
    k = ncol(M)-4
    W=as.matrix(M[,(k+3):(k+4)])
  }else{
    k=ncol(M)-3
    W=as.matrix(M[,(k+3)])
  }
  l = 2*(k+2)
  v = k+5
  
  parlgamma <- (length(par) - (l+8))/2
  
  beta = as.matrix(par[1:k])
  alphaT.1 = par[k+1]
  alphaT.2 = par[k+2]
  lambdaT.1 = par[k+3]
  lambdaT.2 = par[k+4]
  eta = as.matrix(par[v:l])
  alphaC.1 = par[l+1]
  alphaC.2 = par[l+2]
  lambdaC.1 = par[l+3]
  lambdaC.2 = par[l+4]
  sigma1 = par[l+5]
  sigma2 = par[l+6]
  theta_1 = par[l+7]
  theta_2 = par[l+8]
  gamma1 = as.matrix(par[(l+9):(l+8+parlgamma)])
  gamma2 = as.matrix(par[(l+9+parlgamma):(l+8+2*parlgamma)])
  
  
  X=as.matrix(M[,1:k])
  Z1=as.matrix(M[,(k+1)])
  Z2=as.matrix(M[,(k+2)])
  XandW=as.matrix(cbind(X,W))
  
  a = XandW%*%(gamma1-gamma2)
  b = XandW%*%gamma1
  c = XandW%*%gamma2
  
  E1.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
  E2.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
  E3.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))
  
  E1.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
  E2.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
  E3.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))
  
  
  Vest1 = Z1*E1.nu1+Z2*E2.nu1+(1-Z1-Z2)*E3.nu1
  Vest2 = Z1*E1.nu2+Z2*E2.nu2+(1-Z1-Z2)*E3.nu2
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(X%*%beta+Z1*alphaT.1+Z2*alphaT.2+Vest1*lambdaT.1+Vest2*lambdaT.2))/sigma1
  z2 = (transY.C-(X%*%eta+Z1*alphaC.1+Z2*alphaC.2+Vest1*lambdaC.1+Vest2*lambdaC.2))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z2,cov12=0))^(1-(Delta+Xi)))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}


######################## simulation function ###################################


SimulationCI22_MEV <- function(n, nsim, iseed, init.value.theta_1, init.value.theta_2,
                               part.to.evaluate, number.of.parts) {
  
  if (nsim %% number.of.parts != 0) {
    stop("nsim needs to be a multiple of number.of.parts.")
  }
  if ((part.to.evaluate > number.of.parts) || (part.to.evaluate <= 0)) {
    stop("part.to.evaluate is not valid.")
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
  
  # Create the appropriate set of i's to check
  part.size <- nsim / number.of.parts
  i.to.check <- 1:nsim
  i.to.check <- i.to.check[(part.size*(part.to.evaluate-1) + 1):(part.size*part.to.evaluate)]
  
  package.vector <- c("MASS", "nloptr", "numDeriv", "VGAM", "nnet")
  export.vector <- c("parN", "nEV", "parlgamma")
  
  result_list <- foreach(i = i.to.check,
                         .packages = package.vector,
                         .export = export.vector) %dopar% 
    {
                           
    source("Functions_multipleEV.R")
    
    # cat((i - i.to.check[1] + 1)*100/length(i.to.check), "% Completion \n")
    
    parl = length(parN[[1]])
    totparl = 2*parl
    # parlgamma = (parl-nEV)
    
    data = dat.sim.reg.Zbin(n,parN,iseed+i)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl-1))]
    Z1 = data[,parl]
    Z2 = data[,(parl+1)]
    Z = data[,(parl:(parl+1))]
    W = data[,(parl+2):(parl+3)]
    XandW = cbind(data[,4],X,W)
    
    
    # gammaest <- nloptr(x0=rep(0,2*parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma*2)),ub=c(rep(Inf,parlgamma*2)),
    #                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution

    end_var_class <- relevel(as.factor(Z1 + 2*Z2), ref = 1)
    gammaest <- summary(multinom(end_var_class ~ -1 + XandW, ref = 1, trace = FALSE))$coefficients
    gammaest <- c(gammaest[1,], gammaest[2,])
    
    gamma1 <- gammaest[1:parlgamma]  
    gamma2 <- gammaest[(parlgamma+1):(2*parlgamma)]
    
    a = XandW%*%(gamma1-gamma2)
    b = XandW%*%gamma1
    c = XandW%*%gamma2
    
    E1.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
    E2.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
    E3.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))
    
    E1.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
    E2.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
    E3.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))
    
    
    V1 = Z1*E1.nu1+Z2*E2.nu1+(1-Z1-Z2)*E3.nu1
    V2 = Z1*E1.nu2+Z2*E2.nu2+(1-Z1-Z2)*E3.nu2
    
    
    # Estimated V
    M = cbind(data[,4:(1+parl)],V1,V2)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(1+parl)],data[,(ncol(data)-1)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-((ncol(M)-1):ncol(M))]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:4] : beta
    # [5:8] : eta
    # [9]   : sigma1
    # [10]  : sigma2
    # [11]  : rho
    # [12]  : theta_1
    # [13]  : theta_2
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(parl-1)]
    initE = initE[-(2*parl-2)]
    initE = initE[-(2*parl-3)]
    
    # Append theta's to initE and replace the original theta_1 (now third-to-last
    # element) with the initial value for rho.
    initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-4)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-4)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)))
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
    t_s1 = 1/parhatE[totparl-3]*se1[totparl-3]
    t_s2 = 1/parhatE[totparl-2]*se1[totparl-2]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-3])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-3])+1.96*t_s1 
    ms2_l = log(parhatE[totparl-2])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl-2])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl-1])/(1-parhatE[totparl-1])))     # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl-1]^2))*se1[totparl-1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    # Confidence interval for theta
    
    r1theta1_l <- parhatE[length(parhatE)-1] - 1.96 * se1[length(parhatE)-1]
    r1theta1_u <- parhatE[length(parhatE)-1] + 1.96 * se1[length(parhatE)-1]
    r1theta2_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
    r1theta2_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
    
    # Matrix of all the confidence intervals
    EC2 = cbind(matrix(c(parhatE[1:(totparl-4)]-1.96*(se1)[1:(totparl-4)],S1_l,S2_l,r1_l, r1theta1_l, r1theta2_l),ncol=1),
                matrix(c(parhatE[1:(totparl-4)]+1.96*(se1)[1:(totparl-4)],S1_u,S2_u,r1_u, r1theta1_u, r1theta2_u),ncol=1)) 
    
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (6 params) = First 6 params of parhat1
    # - eta (6 params) = Next 6 params of parhat1
    # - sigma1 = parhat1[13]
    # - sigma2 = parhat1[14]
    # - rho = 0
    # - theta_1 = parhat1[15]
    # - theta_2 = parhat1[16]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Stack the vectors of estimates.
    parhatG = c(parhat,as.vector(gamma1),as.vector(gamma2))
    # - beta (6 params) = First 6 params of parhat
    # - eta (6 params) = Next 6 params of parhat
    # - sigma1 = parhat[13]
    # - sigma2 = parhat[14]
    # - rho = parhat[15]
    # - theta_1 = parhat[16]
    # - theta_2 = parhat[17]
    # - gamma1 = (intercept, gamma_X, gamma_W1, gamma_W2)
    # - gamma2 = (intercept, gamma_X, gamma_W1, gamma_W2)
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,discrete=1,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 

    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
    # (i.e. H_delta).
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+2*parlgamma)]
    
    #-M matrix 
    WM <- matrix(rep(0,(parlgamma*2)^2),nrow=parlgamma*2)
    for (i in 1:n){
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      matrix <- matrix(c(p1*(1-p1),-p1*p2, -p1*p2, p2*(1-p2)),ncol=2)
      xx <- XandW[i,]%*%t(XandW[i,])
      WM <- WM+kronecker(matrix,xx)
    }

    MI = ginv(WM)
    
    Epartvar2 = 0
    
    for (i in 1:n){
      # h_m matrix
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      Zp <- c(Z1[i]-p1, Z2[i]-p2)
      h_mi <- kronecker(Zp,XandW[i,])
      
      # psi matrix
      psii = MI%*%h_mi
      
      # h_l matrix
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))

      partvar = t(J1) + Vargamma%*%psii
      Epartvar2 = Epartvar2+(partvar%*%t(partvar))
      
    }
    
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
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (6 params) = First 6 params of parhat1
    # - eta (6 params) = Next 6 params of parhat1
    # - sigma1 = parhat1[13]
    # - sigma2 = parhat1[14]
    # - rho = 0
    # - theta_1 = parhat1[15]
    # - theta_2 = parhat1[16]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
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
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta_1
    # [12]    : param for theta_2
    # [13:15] : params for (gamma_0, gamma_X, gamma_W)
    parhatGI = c(parhat1,as.vector(gamma1),as.vector(gamma2))
    
    HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,discrete=1,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+(2*parlgamma)-1)]
   
    
    
    Epartvar2I = 0
    
    for (i in 1:n){
      # h_m matrix
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      Zp <- c(Z1[i]-p1, Z2[i]-p2)
      h_mi <- kronecker(Zp,XandW[i,])
      
      # psi matrix
      psii = MI%*%h_mi
      
      # h_l matrix
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      
      partvarI = t(J1I) + VargammaI%*%psii
      Epartvar2I = Epartvar2I+(partvarI%*%t(partvarI))
      
    }

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
    
    rItheta1_l <- parhat1[length(parhat1)-1] - 1.96 * seI[length(parhat1)-1]
    rItheta1_u <- parhat1[length(parhat1)-1] + 1.96 * seI[length(parhat1)-1]
    rItheta2_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
    rItheta2_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta1_l, rItheta2_l),ncol=1),
                matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta1_u, rItheta2_u), ncol=1))
    
    return(list(c(parhat,se,c(t(EC1))), c(parhatE,se1,c(t(EC2))),
                c(parhatre,sere,c(t(EC3))), c(parhat1,seI,c(t(EC4)))))
  }
  
  # Combine all separate results
  for (i in 1:length(result_list)) {
    results = rbind(results, result_list[[i]][[1]])
    results1 = rbind(results1, result_list[[i]][[2]])
    results2 = rbind(results2, result_list[[i]][[3]])
    results3 = rbind(results3, result_list[[i]][[4]])
  }
  
  percentage1 <- per/(n*length(i.to.check))     #percentage of censoring
  percentage2 <- per2/(n*length(i.to.check))
  
  df.estV <- data.frame(results, row.names = i.to.check)
  df.naive <- data.frame(results1, row.names = i.to.check)
  df.realV <- data.frame(results2, row.names = i.to.check)
  df.indep <- data.frame(results3, row.names = i.to.check)
  df.percentage <- data.frame(per1 = percentage1, per2 = percentage2,
                              row.names = part.to.evaluate)
  
  write.csv(df.estV, file = paste0("Simulations_multipleEV/Size ",n,"/df_MEV_estV_",
                                   part.to.evaluate, "_out_of_",
                                   number.of.parts, ".csv"),
            row.names = TRUE)
  
  write.csv(df.naive, file = paste0("Simulations_multipleEV/Size ",n,"/df_MEV_naive_",
                                    part.to.evaluate, "_out_of_",
                                    number.of.parts, ".csv"),
            row.names = TRUE)
  
  write.csv(df.realV, file = paste0("Simulations_multipleEV/Size ",n,"/df_MEV_realV_",
                                    part.to.evaluate, "_out_of_",
                                    number.of.parts, ".csv"),
            row.names = TRUE)
  
  write.csv(df.indep, file = paste0("Simulations_multipleEV/Size ",n,"/df_MEV_indep_",
                                    part.to.evaluate, "_out_of_",
                                    number.of.parts, ".csv"),
            row.names = TRUE)
  
  write.csv(df.percentage,
            file = paste0("Simulations_multipleEV/Size ",n,"/df_MEV_percentage_",
                          part.to.evaluate, "_out_of_", number.of.parts, ".csv"),
            row.names = TRUE)
}


Summarize_results_full_MEV <- function(CI) {
  
  # Some input validation
  if (CI != "CI22") {
    stop("Only implemented for CI22")
  }
  
  # For each sample size, create a data set of the results. (The code in this 
  # for-loop is just a copy-paste from the summarize_results function above)
  sum_list = list()
  sum1_list = list()
  sum2_list = list()
  sum3_list = list()
  
  for (samsize_index in 1:length(samsize)) {

    n <- samsize[samsize_index]
    
    # Create correct folder name
    folder.name <- paste0("Simulations_multipleEV/Size ",n)
    
    #
    # Create full data sets
    #
    
    # Get all file names in 'Data 2500 simulations ...' folder
    files <- list.files(folder.name)
    
    # Read all files starting with "df_estV_". Store them in a list object. Do the
    # same for files starting with "df_naive_", "df_realV_", "df_indep" and
    # "df_percentage_"
    data_estV <- list()
    data_naive <- list()
    data_realV <- list()
    data_indep <- list()
    data_percentage <- list()
    
    for (file_name in files) {
      if (grepl("df_MEV_estV_", file_name)) {
        data_estV[[length(data_estV) + 1]] <- read.csv(paste0(folder.name, "/", file_name))
      }
      
      if (grepl("df_MEV_naive_", file_name)) {
        data_naive[[length(data_naive) + 1]] <- read.csv(paste0(folder.name, "/", file_name))
      }
      
      if (grepl("df_MEV_realV_", file_name)) {
        data_realV[[length(data_realV) + 1]] <- read.csv(paste0(folder.name, "/", file_name))
      }
      
      if (grepl("df_MEV_indep_", file_name)) {
        data_indep[[length(data_indep) + 1]] <- read.csv(paste0(folder.name, "/", file_name))
      }
      
      if (grepl("df_MEV_percentage_", file_name)) {
        data_percentage[[length(data_percentage) + 1]] <- read.csv(paste0(folder.name, "/", file_name))
      }
    }
    
    # Create empty data frames
    results.estV <- data_estV[[1]]
    results.naive <- data_naive[[1]]
    results.realV <- data_realV[[1]]
    results.indep <- data_indep[[1]]
    results.percentage <- data_percentage[[1]]
    
    # Append all separate data frames
    if (length(data_estV) > 1) {
      for (i in 2:length(data_estV)) {
        results.estV <- rbind(results.estV, data_estV[[i]])
        results.naive <- rbind(results.naive, data_naive[[i]])
        results.realV <- rbind(results.realV, data_realV[[i]])
        results.indep <- rbind(results.indep, data_indep[[i]])
        results.percentage <- rbind(results.percentage, data_percentage[[i]])
      }
    }
    
    # Sanity check for data: check if each index only appears once. If sanity
    # check is okay, remove the extra index column.
    if (length(results.estV$X) != length(unique(results.estV$X))) {
      stop("Duplicates in results.estV detected")
    } else {
      results.estV <- results.estV[,-1]
    }
    if (length(results.naive$X) != length(unique(results.naive$X))) {
      stop("Duplicates in results.naive detected")
    } else {
      results.naive <- results.naive[,-1]
    }
    if (length(results.realV$X) != length(unique(results.realV$X))) {
      stop("Duplicates in results.realV detected")
    } else {
      results.realV <- results.realV[,-1]
    }
    if (length(results.indep$X) != length(unique(results.indep$X))) {
      stop("Duplicates in results.indep detected")
    } else {
      results.indep <- results.indep[,-1]
    }
    if (length(results.percentage$X) != length(unique(results.percentage$X))) {
      stop("Duplicates in results.percentage detected")
    } else {
      results.percentage <- results.percentage[,-1]
    }
    
    #
    # Results of model with estimated V
    #
    
    # Put all parameters (except gamma) into a vector
    par0 = c(parN[[1]],parN[[2]],parN[[3]])
    par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
    
    # par0:
    # - [1:6] : beta
    # - [7:12]: eta
    # - [13]  : sigma1
    # - [14]  : sigma2
    # - [15]  : rho
    # - [16]  : theta_1
    # - [17]  : theta_2
    #
    # - totparl = 12
    
    # Statistics on the parameter estimates
    Bias = apply(results.estV[,1:(totparl+5)]-par0m,2,mean)
    ESE = apply(results.estV[,1:(totparl+5)],2,sd)
    RMSE = sqrt(apply((results.estV[,1:(totparl+5)]-par0m)^2,2,mean))
    
    # Statistics on the parameter standard deviations
    MSD  = apply(results.estV[,((totparl+5)+1):(2*(totparl+5))],2, mean)
    
    # Statistics on the parameter CI's: for each parameter, check how many times the
    # true value is contained in the estimated confidence interval. We divide by
    # nsim to obtain a percentage.
    CP = rep(0,totparl+5)
    datacp = results.estV[,(2*(totparl+5)+1):(4*(totparl+5))]
    for(i in 1:(totparl+5)) {
      index=c(2*i-1,2*i)
      CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nrow(results.estV)
    } 
    
    sum_list[[samsize_index]] <- cbind(Bias,ESE,MSD,RMSE,CP)
    
    #
    # Model with no V
    #
    
    # Put all parameters (except gamma) into a vector
    par0 = c(parN[[1]],parN[[2]],parN[[3]])
    
    # Remove parameters pertaining to V
    par0 = par0[-(parl)]
    par0 = par0[-(parl-1)]
    par0 = par0[-(2*parl-2)]
    par0 = par0[-(2*parl-3)]
    par0m = matrix(par0,nsim,(totparl+1),byrow=TRUE)
    
    # par0:
    # - [1:4] : beta
    # - [5:8] : eta
    # - [9]   : sigma1
    # - [10]  : sigma2
    # - [11]  : rho
    # - [12]  : theta_1
    # - [13]  : theta_2
    #
    # - totparl = 12
    
    # Statistics on the parameter estimates
    Bias = apply(results.naive[,1:(totparl+1)]-par0m,2,mean)
    ESE = apply(results.naive[,1:(totparl+1)],2,sd)
    RMSE = sqrt(apply((results.naive[,1:(totparl+1)]-par0m)^2,2,mean))
    
    # Statistics on the parameter standard deviations
    MSD  = apply(results.naive[,((totparl+1) + 1):(2*(totparl+1))],2,mean)
    
    # Statistics on the parameter CI's: for each parameter, check how many times the
    # true value is contained in the estimated confidence interval. We divide by
    # nsim to obtain a percentage.
    CP = rep(0,(totparl+1))
    datacp = results.naive[,(2*(totparl+1)+1):(4*(totparl+1))]
    for(i in 1:(totparl+1)){
      index = c(2*i-1,2*i)
      CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nrow(results.naive)
    } 
    
    sum1_list[[samsize_index]] <- cbind(Bias,ESE,MSD,RMSE,CP) 
    
    #
    # Model with real V
    #
    
    par0 = c(parN[[1]],parN[[2]],parN[[3]])
    par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
    # par0:
    # - [1:6] : beta
    # - [7:12]: eta
    # - [13]  : sigma1
    # - [14]  : sigma2
    # - [15]  : rho
    # - [16]  : theta_1
    # - [17]  : theta_2
    #
    # - totparl = 12
    
    # Statistics on the parameter estimates
    Bias = apply(results.realV[,1:(totparl+5)]-par0m,2,mean)
    ESE = apply(results.realV[,1:(totparl+5)],2,sd)
    RMSE = sqrt(apply((results.realV[,1:(totparl+5)]-par0m)^2,2,mean))
    
    # Statistics on the standard deviation estimates
    MSD  = apply(results.realV[,((totparl+5)+1):(2*(totparl+5))],2, mean)
    
    # Statistics on the parameter CI's: for each parameter, check how many times the
    # true value is contained in the estimated confidence interval. We divide by
    # nsim to obtain a percentage.
    CP = rep(0,totparl+5)
    datacp = results.realV[,(2*(totparl+5)+1):(4*(totparl+5))]
    for(i in 1:(totparl+5)) {
      index=c(2*i-1,2*i)
      CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nrow(results.realV)
    } 
    
    sum2_list[[samsize_index]] <- cbind(Bias,ESE,MSD,RMSE,CP) 
    
    #
    # Results of model with estimated V but independence
    #
    
    par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4], parN[[3]][5])
    par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
    # par0:
    # - [1:6] : beta
    # - [7:12]: eta
    # - [13]   : sigma1
    # - [14]  : sigma2
    # - [15]  : theta_1
    # - [16]  : theta_2
    #
    # - totparl = 12
    
    # Statistics on the parameter estimates.
    Bias = apply(results.indep[,1:(totparl+4)]-par0m,2,mean)
    ESE = apply(results.indep[,1:(totparl+4)],2,sd)
    RMSE = sqrt(apply((results.indep[,1:(totparl+4)]-par0m)^2,2,mean))
    
    # Statistics on the standard deviation estimates
    MSD  = apply(results.indep[,((totparl+4)+1):(2*(totparl+4))],2, mean)
    
    CP = rep(0,totparl+4)
    datacp = results.indep[,(2*(totparl+4) + 1):(4*(totparl+4))]
    for (i in 1:(totparl+4)) {
      index=c(2*i-1,2*i)
      CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nrow(results.indep)
    } 
    
    sum3_list[[samsize_index]] <- cbind(Bias,ESE,MSD,RMSE,CP) 
  }
  
  # Remove the MSD columns
  for (i in 1:3) {
    sum_list[[i]] <- subset(sum_list[[i]], select = -MSD)
    sum1_list[[i]] <- subset(sum1_list[[i]], select = -MSD)
    sum2_list[[i]] <- subset(sum2_list[[i]], select = -MSD)
    sum3_list[[i]] <- subset(sum3_list[[i]], select = -MSD)
  }
  
  # Two-step estimator
  sum <- cbind(sum_list[[1]], sum_list[[2]], sum_list[[3]])
  
  # naive estimator
  sum1 <- cbind(sum1_list[[1]], sum1_list[[2]], sum1_list[[3]])
  
  # Oracle estimator
  sum2 <- cbind(sum2_list[[1]], sum2_list[[2]], sum2_list[[3]])
  
  # Independence estimator
  sum3 <- cbind(sum3_list[[1]], sum3_list[[2]], sum3_list[[3]])

  
  ## Results of model with estimated V
  
  colnames(sum) = rep(c("Bias", "ESD", "RMSE", "CR"), 3)
  rownames(sum) = namescoef
  sum_all_names <- namescoef
  
  ## Results of naive model
  colnames(sum1) <- rep(c("Bias", "ESD", "RMSE", "CR"), 3)
  namescoefr <- namescoef[-((parl-1):parl)]
  namescoefr <- namescoefr[-((2*parl-3):(2*parl-2))]
  rownames(sum1) <- namescoefr
  sum_all_names <- c(sum_all_names, namescoefr)
  
  ## Results of independence model
  colnames(sum3) <- rep(c("Bias", "ESD", "RMSE", "CR"), 3)
  namescoefi <- namescoef[-(length(namescoef) - 2)]
  rownames(sum3) <- namescoefi
  sum_all_names <- c(sum_all_names, namescoefi)
  
  
  ## Results of oracle estimator
  colnames(sum2) <- rep(c("Bias", "ESD", "RMSE", "CR"), 3)
  rownames(sum2) <- namescoef
  sum_all_names <- c(sum_all_names, namescoef)
  
  # Stack all results in a single data frame
  sum_all <- as.data.frame(rbind(sum, sum1, sum3, sum2))
  
  # Since the R does not allow duplicate row names, we also store
  # the row names in a column
  sum_all$'' <- sum_all_names
  
  # Reorder columns so to put names column in front
  sum_all <- sum_all[c(13, 1:12)]
  
  # Make nice Latex table
  xtab <- xtable(sum_all)
  align(xtab) <- c("c|c|cccc|cccc|cccc|")
  
  # set to 3 significant digits
  digits(xtab) = rep(3,14)
  
  header <- paste0('\\hline \\multicolumn{5}{|c}{$n = 500$} & \\multicolumn{4}{|c}{$n = 1000$}',
                   '& \\multicolumn{4}{|c|}{$n = 2000$} \\\\')
  title_2step <- "\\hline \\multicolumn{13}{|c|}{two-step estimator} \\\\ "
  title_naive <- "\\hline \\multicolumn{13}{|c|}{naive estimator} \\\\ \\hline "
  title_indep <- "\\hline \\multicolumn{13}{|c|}{independence estimator} \\\\ \\hline "
  title_oracle <- "\\hline \\multicolumn{13}{|c|}{oracle estimator} \\\\ \\hline "
  fit_to_page <- "\begin{adjustbox}{width=1\textwidth}"
  addtorow <- list()
  addtorow$pos <- list(-1, -1, 13, 24, 36)
  addtorow$command <- c(header, title_2step, title_naive, title_indep, title_oracle)
  
  # For correct naming of the file
  numbers <- substr(CI, 3, 4)
  
  
  
  ##############################################################################
  # Note that the generated table will result in the file not compiling. This  #
  # is because the 'adjustbox'-environment needs the argument                  #
  #                                                                            #
  # \begin{adjustbox}{width=\linewidth}                                        #
  #                                                                            #
  # but unfortunately, print.xtable does not allow for this argument to be     #
  # specified.                                                                 #
  ##############################################################################
  
  
  print.xtable(xtab, include.rownames = F,
               add.to.row = addtorow, append=TRUE, table.placement = "H",
               sanitize.text.function = function(x){x},
               sanitize.colnames.function = function(x){gsub("[[:punct:], [:digit:]]", "", x)},
               latex.environments = c("adjustbox", "center"),
               comment = F)
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab, include.rownames = F,
               add.to.row = addtorow, append=TRUE, table.placement = "H",
               sanitize.text.function = function(x){x},
               sanitize.colnames.function = function(x){gsub("[[:punct:], [:digit:]]", "", x)},
               latex.environments = c("adjustbox", "center"),
               file = paste0("Simulation results_multipleEV/YJ_ad_all_CI22.txt"))
}

######################## Application ###########################################

DataApplicationChess.multipleEV = function(data,init.value.theta_1, init.value.theta_2) {

  
  n = nrow(data)
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  intercept <- data[,4]
  X = data[,(5:(parl-1))]
  Z1 = ifelse(data[,parl]==2,1,0) # Z = 1 is the reference level.
  Z2 = ifelse(data[,parl]==3,1,0)
  Z = cbind(Z1,Z2)
  W = data[,parl + 1]
  XandW = cbind(intercept,X,W)
  
  # gammaest <- nloptr(x0=rep(0,2*parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma*2)),ub=c(rep(Inf,parlgamma*2)),
  #                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  ##############################################################################
  # Maybe better to estimate with a known function since the above does not seem
  # to work very well.
  
  end_var_class <- relevel(as.factor(data[,parl]), ref = 1)
  gammaest <- summary(multinom(end_var_class ~ -1 + XandW, ref = 1, trace = FALSE))$coefficients
  gammaest <- c(gammaest[1,], gammaest[2,])
  
  ##############################################################################
  
  gamma1 <- gammaest[1:parlgamma]  
  gamma2 <- gammaest[(parlgamma+1):(2*parlgamma)]
    
  a = XandW%*%(gamma1-gamma2)
  b = XandW%*%gamma1
  c = XandW%*%gamma2
    
  E1.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(a/(1+exp(-a)+exp(-b))-1/(1+exp(-b))*log(1+exp(a)*(1+exp(-b))))
  E2.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2)*(1+exp(-c)))*(-a/(exp(-a)+exp(-a-c)+1)+log(exp(a)+exp(-c)+1))
  E3.nu1 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))*(1/(1+exp(c))*(-c-(b-c)/(exp(-b)+exp(c-b)+1)+log(1+exp(c)+exp(b)))-log(1+exp(-c)+exp(b-c))+log(1+exp(b-c)+exp(-c))/(1+exp(-b)))
    
  E1.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma1))*(b/(1+exp(-a)+exp(-b))-1/(1+exp(-a))*log(1+exp(b)*(1+exp(-a))))
  E2.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(exp(XandW%*%gamma2))*(1/(1+exp(-c))*(c-(a+c)/(exp(-a)+exp(-c-a)+1)+log(1+exp(-c)+exp(a)))-log(1+exp(c)+exp(a+c))+log(1+exp(a+c)+exp(c))/(1+exp(-a)))
  E3.nu2 = (1+exp(XandW%*%gamma1)+exp(XandW%*%gamma2))/(1+exp(c))*(-b/(exp(-b)+exp(-b+c)+1)+log(exp(b)+exp(c)+1))
    
    
  V1 = Z1*E1.nu1+Z2*E2.nu1+(1-Z1-Z2)*E3.nu1
  V2 = Z1*E1.nu2+Z2*E2.nu2+(1-Z1-Z2)*E3.nu2
    
    
  # Estimated V
  M = cbind(intercept,X,Z,V1,V2)
    
  # No V (using W instead)
  MnoV = cbind(intercept,X,Z,W)
    
    
  # Assign starting values:
  # - beta = zero-vector
  # - eta = zero-vector
  # - sigma1 = 1
  # - sigma2 = 1 
  # - theta_1 = init.value.theta_1
  # - theta_2 = init.value.theta_2
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
  # Independent model for starting values sigma and theta.
  #
  # Note the difference with the version of Gilles: the likelihood function now
  # takes an extra argument (= theta), so the vector of initial values needs
  # to take this into account. Also the vectors for the lower -and upper bound
  # of the parameters ('lb' and 'ub') should take this into account. Note that
  # theta is a value between 0 and 2.
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  # Model with no V --> remove data for v from data matrix
  ME = M[,-((ncol(M)-1):ncol(M))]
    
  # Remove coefficients for v
  initE = parhat1[-parl]
  initE = initE[-(parl-1)]
  initE = initE[-(2*parl-2)]
  initE = initE[-(2*parl-3)]
    
  # Append theta's to initE and replace the original theta_1 (now third-to-last
  # element) with the initial value for rho.
  initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
  initE[length(initE) - 2] <- 0
    
  # Again we make sure to properly adapt the upper -and lower bound values of
  # theta.
  parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-4)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-4)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  H1I = ginv(H1)
  se1 = sqrt(abs(diag(H1I)));
    
  # Delta method variance (makes sure no negative values in CI for variance)
  # --> Take the log of the estimates, construct CI for the log-estimates and
  #     backtransform to the original estimates by exponentiating.
    
  t_s1 = 1/parhatE[totparl-3]*se1[totparl-3]
  t_s2 = 1/parhatE[totparl-2]*se1[totparl-2]
    
  # Conf. interval for transf. sigma's
    
  ms1_l = log(parhatE[totparl-3])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-3])+1.96*t_s1 
  ms2_l = log(parhatE[totparl-2])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl-2])+1.96*t_s2 
    
  # Back transform
    
  S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
  # Confidence interval for rho
    
  z1t = 0.5*(log((1+parhatE[totparl-1])/(1-parhatE[totparl-1])))     # Fisher's z transform
  se1_z = (1/(1-parhatE[totparl-1]^2))*se1[totparl-1]
  z1t_l = z1t-1.96*(se1_z)
  z1t_u = z1t+1.96*(se1_z)
    
    
  # Back transform
    
  r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
  r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
  # Confidence interval for theta
    
  r1theta1_l <- parhatE[length(parhatE)-1] - 1.96 * se1[length(parhatE)-1]
  r1theta1_u <- parhatE[length(parhatE)-1] + 1.96 * se1[length(parhatE)-1]
  r1theta2_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
  r1theta2_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
    
  # Matrix of all the confidence intervals
  EC2 = cbind(matrix(c(parhatE[1:(totparl-4)]-1.96*(se1)[1:(totparl-4)],S1_l,S2_l,r1_l, r1theta1_l, r1theta2_l),ncol=1),
                matrix(c(parhatE[1:(totparl-4)]+1.96*(se1)[1:(totparl-4)],S1_u,S2_u,r1_u, r1theta1_u, r1theta2_u),ncol=1)) 
    
    
  # Model with estimated V

    
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
    
  # Again we make sure to properly adapt the upper -and lower bound values of
  # theta.
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  parhatG = c(parhat,as.vector(gamma1),as.vector(gamma2))

  Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,discrete=0,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    
  # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
  # (i.e. H_delta).
  H = Hgamma[1:length(initd),1:length(initd)]
  HI = ginv(H)
    
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+2*parlgamma)]
    
  #-M matrix 
  WM <- matrix(rep(0,(parlgamma*2)^2),nrow=parlgamma*2)
  for (i in 1:n){
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      matrix <- matrix(c(p1*(1-p1),-p1*p2, -p1*p2, p2*(1-p2)),ncol=2)
      xx <- XandW[i,]%*%t(XandW[i,])
      WM <- WM+kronecker(matrix,xx)
  }
    
  MI = ginv(WM)
    
  Epartvar2 = 0
    
  for (i in 1:n){
      # h_m matrix
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      Zp <- c(Z1[i]-p1, Z2[i]-p2)
      h_mi <- kronecker(Zp,XandW[i,])
      
      # psi matrix
      psii = MI%*%h_mi
      
      # h_l matrix
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      
      partvar = t(J1) + Vargamma%*%psii
      Epartvar2 = Epartvar2+(partvar%*%t(partvar))
      
  }
    
    
    
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
    
  rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
  rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
  rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
  rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
  # Matrix with all confidence intervals
  EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
  # Model with estimated V but assuming independence

  parhatGI = c(parhat1,as.vector(gamma1),as.vector(gamma2))
    
  HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,discrete=0,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
  HIInd = ginv(HInd)
    
  VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+(2*parlgamma)-1)]
    
    
    
  Epartvar2I = 0
    
  for (i in 1:n){
      # h_m matrix
      p1 <- exp(XandW[i,]%*%gamma1)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      p2 <- exp(XandW[i,]%*%gamma2)/(1+exp(XandW[i,]%*%gamma1)+exp(XandW[i,]%*%gamma2))
      Zp <- c(Z1[i]-p1, Z2[i]-p2)
      h_mi <- kronecker(Zp,XandW[i,])
      
      # psi matrix
      psii = MI%*%h_mi
      
      # h_l matrix
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      
      partvarI = t(J1I) + VargammaI%*%psii
      Epartvar2I = Epartvar2I+(partvarI%*%t(partvarI))
      
  }
    
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
    
  rItheta1_l <- parhat1[length(parhat1)-1] - 1.96 * seI[length(parhat1)-1]
  rItheta1_u <- parhat1[length(parhat1)-1] + 1.96 * seI[length(parhat1)-1]
  rItheta2_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
  rItheta2_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
    
  EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta1_l, rItheta2_l),ncol=1),
                matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta1_u, rItheta2_u), ncol=1))


  
  
  
  # Results of model assuming confounding and dependence between T and C.
  pvalue <- 2*pmin((1-pnorm(parhat/se)),pnorm(parhat/se))
  significant <- ifelse(pvalue < 0.10,
                        ifelse(pvalue < 0.05,
                               ifelse(pvalue < 0.01, "**", "*"),"."), "")
  results.confound_dep <- cbind(parhat, se, pvalue, EC1)
  colnames(results.confound_dep) <- c("Estimate", "St.Dev.", "p", "CI.lb", "CI.ub")
  rownames(results.confound_dep) <- namescoef
  
  summary <- data.frame(round(results.confound_dep, 4))
  summary$sign <- significant
  summary <- summary[,c(1:3, 6, 4:5)]
  summary
  
  # Results of naive model
  pvalue.naive <- 2*pmin((1-pnorm(parhatE/se1)),pnorm(parhatE/se1))
  significant.naive <- ifelse(pvalue.naive < 0.10,
                              ifelse(pvalue.naive < 0.05,
                                     ifelse(pvalue.naive < 0.01, "**", "*"),"."), "")
  results.naive <- cbind(parhatE, se1, pvalue.naive, EC2)
  colnames(results.naive) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(parl-1)]
  namescoefr=namescoefr[-(2*parl-2)]
  namescoefr=namescoefr[-(2*parl-3)]
  
  rownames(results.naive) <- namescoefr
  
  summary1 <- data.frame(round(results.naive, 4))
  summary1$sign <- significant.naive
  summary1 <- summary1[,c(1:3, 6, 4:5)]
  summary1
  
  # Results of independence model
  pvalue.indep <- 2*pmin((1-pnorm(parhat1/seI)),pnorm(parhat1/seI))
  significant.indep <- ifelse(pvalue.indep < 0.10,
                              ifelse(pvalue.indep < 0.05,
                                     ifelse(pvalue.indep < 0.01, "**", "*"),"."), "")
  results.indep <- cbind(parhat1, seI, pvalue.indep, EC4)
  colnames(results.indep) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  rownames(results.indep) <- namescoef[-(length(namescoef) - 2)]
  
  summary2 <- data.frame(round(results.indep, 4))
  summary2$sign <- significant.indep
  summary2 <- summary2[,c(1:3, 6, 4:5)]
  summary2
  
  ## Create LaTeX tables of results
  xtab = xtable(summary, digits = 3)
  header= c("sample size",n,"Results 2-step_Estimation with YT-transformation multiple EV")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_2-step_Estimation_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary1, digits = 3)
  header= c("sample size",n,"Results naive model with YT-transformation multiple EV")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_naive_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary2, digits = 3)
  header= c("sample size",n,"Results independence model with YT-transformation multiple EV")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  return(list(parhat,gamma1,gamma2))
  
  
}