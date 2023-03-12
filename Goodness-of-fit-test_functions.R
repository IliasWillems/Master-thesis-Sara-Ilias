
GOF_test <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # structure of 'initd':
  #
  # [1:4] beta (4 params) = First 4 params of parhat1
  # [5:8] eta (4 params) = Next 4 params of parhat1
  # [9]   sigma1 = parhat1[9]
  # [10]  sigma2 = parhat1[10]
  # [11]  rho = 0
  # [12]  theta_1 = parhat1[11]
  # [13]  theta_2 = parhat1[12]
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1 or less.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      } else {
        break
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    surv_km[["time"]] <- c(surv_km[["time"]][1:matched_until], surv_km[["time"]][matched_until:length(surv_km[["time"]])])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    for (i in 1:n) {
      A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1)
    cdf_km.star <- 1 - surv_km.star[["surv"]]
    
    # Same issue as before
    while (length(cdf_km.star) < length(FK.star)) {
      
      Y.star.sorted <- sort(Y.star)
      matched_until <- 0
      for (i in 1:length(surv_km.star[["time"]])) {
        if (Y.star.sorted[i] == surv_km.star[["time"]][i]) {
          matched_until <- i
        } else {
          break
        }
      }
      
      surv_km.star[["time"]] <- c(surv_km.star[["time"]][1:matched_until], surv_km.star[["time"]][matched_until:length(surv_km.star[["time"]])])
      cdf_km.star <- c(cdf_km.star[1:matched_until], cdf_km.star[matched_until:length(cdf_km.star)])
    }
    
    TCMb_vector[b] <- n*sum((cdf_km.star - FK.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}


GOF_test_parallel <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[totparl+4])
  Y_trans_theta2 <- YJtrans(Y, parhat[totparl+5])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:parl])/parhat[totparl+1]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[(parl+1):totparl])/parhat[totparl+2]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[totparl+3]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1 or less.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      } else {
        break
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    surv_km[["time"]] <- c(surv_km[["time"]][1:matched_until], surv_km[["time"]][matched_until:length(surv_km[["time"]])])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  package.vector <- c("MASS", "nloptr", "VGAM", "pbivnorm", "survival")
  export.vector <- c("parl", "totparl")
  TCMb_vector <- foreach(b = 1:B,
                         .packages = package.vector,
                         .export = export.vector,
                         .combine = c) %dopar% {

    source("Functions_ad.R")
    source("Goodness-of-fit-test_functions.R")
    
    MAX_SEED_SIZE <- 2147483647
    
    # Create bootstrap data sample
    bootstrapseed <- (B*iseed + b) %% MAX_SEED_SIZE
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    
    # If there is not administrative censoring in the data, then the Kaplan - 
    # Meier estimator 'surv_A' will always be 1. In this case, we set all
    # simulated administrative censoring times equal to infinity. Otherwise,
    # proceed normally.
    if (all(surv_A$surv == 1)) {
      A.star <- rep(Inf, n)
    } else {
      for (i in 1:n) {
        A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
      }
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[totparl+4])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[totparl+5])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:parl])/parhat.star[totparl+1]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[(parl + 1):totparl])/parhat.star[totparl+2]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[totparl+3]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1)
    cdf_km.star <- 1 - surv_km.star[["surv"]]
    
    # Same issue as before
    while (length(cdf_km.star) < length(FK.star)) {
      
      Y.star.sorted <- sort(Y.star)
      matched_until <- 0
      for (i in 1:length(surv_km.star[["time"]])) {
        if (Y.star.sorted[i] == surv_km.star[["time"]][i]) {
          matched_until <- i
        } else {
          break
        }
      }
      
      surv_km.star[["time"]] <- c(surv_km.star[["time"]][1:matched_until], surv_km.star[["time"]][matched_until:length(surv_km.star[["time"]])])
      cdf_km.star <- c(cdf_km.star[1:matched_until], cdf_km.star[matched_until:length(cdf_km.star)])
    }
    
    # Output of loop iteration
    n*sum((cdf_km.star - FK.star)^2 * w.star)
    
  } 
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}


GOF_test_noBootstrap <- function(data, Zbin, Wbin, plot.comparison) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # structure of 'initd':
  #
  # [1:4] beta (4 params) = First 4 params of parhat1
  # [5:8] eta (4 params) = Next 4 params of parhat1
  # [9]   sigma1 = parhat1[9]
  # [10]  sigma2 = parhat1[10]
  # [11]  rho = 0
  # [12]  theta_1 = parhat1[11]
  # [13]  theta_2 = parhat1[12]
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on the KM-estimator
  
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # Plot both estimates if needed
  if (plot.comparison) {
    plot(surv_km[["time"]], surv_km[["surv"]], type = 's', xlab = "t (time)",
         ylab = "S(t)", main = "Comparison of the CDF's used in CM statistic")
    lines(surv_km[["time"]], 1 - FK, type = 'l', col = 'red')
    legend("topright", c("nonparametric", "model"), col = c("black", "red"), lty = 1)
  }

  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # 0.90-quantile
  cutoff90 <- 0.34730
  
  # 0.95-quantile
  cutoff95 <- 0.46136
  
  signif90 <- 0
  signif95 <- 0
  if (TCM > cutoff90) {
    signif90 <- 1
  }
  if (TCM > cutoff95) {
    signif95 <- 1
  }
  
  list(signif90, signif95)
}


GOF_EstimateRunDuration <- function(data, B, nruns, Zbin, Wbin, parallel) {
  testseed <- 123321
  start.time <- Sys.time()
  if (!parallel) {
    out <- GOF_test(data, B, testseed, Zbin, Wbin, FALSE)
  } else {
    clust <- makeCluster(10)
    registerDoParallel(clust)
    
    out <- GOF_test_parallel(data, B, testseed, Zbin, Wbin, FALSE)
    
    stopCluster(clust)
  }
  
  diff <- difftime(Sys.time(), start.time, units = "secs")
  diff <- as.integer(round(diff))
  duration <- round((diff*nruns)/60)
  if (duration > 60) {
    message("The simulation will take approximately ", floor(duration / 60),
            " hours and ", duration %% 60, " minutes")
  } else {
    message("The simulation will take approximately ", duration, " minutes")
  } 
}


GOF_noBootstrap_EstimateRunDuration <- function(data, nruns, Zbin, Wbin) {
  start.time <- Sys.time()
  out <- GOF_test_noBootstrap(data, Zbin, Wbin, plot.comparison = F)
  diff <- difftime(Sys.time(), start.time, units = "secs")
  diff <- as.integer(round(diff))
  duration <- round((diff*nruns)/60)
  if (duration > 60) {
    message("The simulation will take approximately ", floor(duration / 60),
            " hours and ", duration %% 60, " minutes")
  } else {
    message("The simulation will take approximately ", duration, " minutes")
  } 
}


GOF_SimulationTypeIerror <- function(parN, nruns, B, iseed, Zbin, Wbin, parallel) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    data <- dat.sim.reg(n, parN, iseed + run, Zbin, Wbin)
    
    if (parallel) {
      out <- GOF_test_parallel(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    } else {
      out <- GOF_test(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    }
    
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


GOF_SimulationTypeIerror_newapproach <- function(parN, nruns, B, iseed, Zbin, Wbin, parallel) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    data <- dat.sim.reg(n, parN, iseed + run, Zbin, Wbin)
    
    if (parallel) {
      out <- GOF_test_newapproach_parallel(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    } else {
      out <- GOF_test_newapproach(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    }
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


GOF_SimulationMisspecification <- function(type, par, nruns, B, iseed, Zbin, Wbin, parallel) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    if (parallel) {
      if (type == "skew") {
        data.skew <- data.misspecified.skew(n, par, iseed + run, Zbin, Wbin)[[1]]
        out <- GOF_test_parallel(data.skew, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "t") {
        data.t <- data.misspecified.t(n, par, iseed + run, Zbin, Wbin)
        out <- GOF_test_parallel(data.t, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "heteroscedastic") {
        data.heteroscedastic <- data.misspecified.heteroscedastic(n, par, iseed + run, Zbin, Wbin)
        out <- GOF_test_parallel(data.heteroscedastic, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'probit') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as probit control function requires ", 
                  "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.probit(n, par, iseed + run, Wbin)
        out <- GOF_test_parallel(data.probit, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'cloglog') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as cloglog control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.cloglog(n, par, iseed + run, Wbin)
        out <- GOF_test_parallel(data.probit, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else {
        stop(paste0("'type' parameter must be either 'skew', 't', 'heteroscedastic'",
                    ", 'probit' or 'cloglog'."))
      }
    } else {
      if (type == "skew") {
        data.skew <- data.misspecified.skew(n, par, iseed + run, Zbin, Wbin)[[1]]
        out <- GOF_test(data.skew, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "t") {
        data.t <- data.misspecified.t(n, par, iseed + run, Zbin, Wbin)
        out <- GOF_test(data.t, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "heteroscedastic") {
        data.heteroscedastic <- data.misspecified.heteroscedastic(n, par, iseed + run, Zbin, Wbin)
        out <- GOF_test(data.heteroscedastic, B, iseed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'probit') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as probit control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.probit(n, par, iseed + run, Wbin)
        out <- GOF_test(data.probit, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'cloglog') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as cloglog control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.cloglog <- data.misspecified.cloglog(n, par, iseed + run, Wbin)
        out <- GOF_test(data.cloglog, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
        
      } else {
        stop(paste0("'type' parameter must be either 'skew', 't', 'heteroscedastic'",
                    ", 'probit' or 'cloglog'."))
        
      }
    }
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
    }
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


GOF_SimulationNoBootstrap_typeIerror <- function(parN, nruns, iseed, Zbin, Wbin) {
  nbr_reject90_noboot <- 0
  nbr_reject95_noboot <- 0
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    data <- dat.sim.reg(n, parN, iseed + run, Zbin, Wbin)
    out <- GOF_test_noBootstrap(data, Zbin, Wbin, plot.comparison = F)
    
    nbr_reject90_noboot <- nbr_reject90_noboot + out[[1]]
    nbr_reject95_noboot <- nbr_reject95_noboot + out[[2]]
  }
  
  list(reject90 = nbr_reject90_noboot/nruns, reject95 = nbr_reject95_noboot/nruns)
}



# These functions use the first approach considered to obtain the Cramer -
# von Mises statistic. This approach will model F_Y instead of F_K.
# --> It was shown that they do not work and therefore shouldn't be used. These
#     functions will also no longer be updated.
GOF_test_newapproach <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # structure of 'initd':
  #
  # [1:4] beta (4 params) = First 4 params of parhat1
  # [5:8] eta (4 params) = Next 4 params of parhat1
  # [9]   sigma1 = parhat1[9]
  # [10]  sigma2 = parhat1[10]
  # [11]  rho = 0
  # [12]  theta_1 = parhat1[11]
  # [13]  theta_2 = parhat1[12]
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_Y(y; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FY <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FY[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
    
    # Take into account the administrative censoring. Note that, in theory, the
    # value for Y[i] should match exactly to some value in surv_A[["time"]], but
    # probably due to rounding, this is not always the case. Hence the necessity
    # for slightly complicated code below. The maximum with 1 is taken to keep
    # the code from throwing an error should, for some reason (maybe again due
    # to rounding), there would not be a value in surv_A[["time"]] which is
    # smaller than or equal to Y[i].
    FY[i] <- FY[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]) + 
      cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]
  }
  
  # Order the values
  FY <- FY[order(FY)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FY[1]
  w[2:n] = FY[2:n]-FY[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  ord <- 1:n
  EFY <- (ord-1)/n
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((EFY - FY)^2 * w)
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    for (i in 1:n) {
      A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FY.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FY.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
      
      # Take into account the administrative censoring
      FY.star[i] <- FY.star[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]) + 
        cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]
    }
    
    FY.star <- FY.star[order(FY.star)]
    
    w.star = rep(0,n)
    w.star[1] = FY.star[1]
    w.star[2:n] = FY.star[2:n]-FY.star[1:(n-1)]
    
    ord <- 1:n
    EFY.star <- (ord-1)/n
    
    TCMb_vector[b] <- n*sum((EFY.star - FY.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}
GOF_test_newapproach_parallel <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # structure of 'initd':
  #
  # [1:4] beta (4 params) = First 4 params of parhat1
  # [5:8] eta (4 params) = Next 4 params of parhat1
  # [9]   sigma1 = parhat1[9]
  # [10]  sigma2 = parhat1[10]
  # [11]  rho = 0
  # [12]  theta_1 = parhat1[11]
  # [13]  theta_2 = parhat1[12]
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_Y(y; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FY <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FY[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
    
    # Take into account the administrative censoring. Note that, in theory, the
    # value for Y[i] should match exactly to some value in surv_A[["time"]], but
    # probably due to rounding, this is not always the case. Hence the necessity
    # for slightly complicated code below. The maximum with 1 is taken to keep
    # the code from throwing an error should, for some reason (maybe again due
    # to rounding), there would not be a value in surv_A[["time"]] which is
    # smaller than or equal to Y[i].
    FY[i] <- FY[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]) + 
      cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]
  }
  
  # Order the values
  FY <- FY[order(FY)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FY[1]
  w[2:n] = FY[2:n]-FY[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  ord <- 1:n
  EFY <- (ord-1)/n
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((EFY - FY)^2 * w)
  
  package.vector <- c("MASS", "nloptr", "VGAM", "pbivnorm", "survival")
  export.vector <- c("parl", "totparl")
  TCMb_vector <- foreach(b = 1:B,
                         .packages = package.vector,
                         .export = export.vector,
                         .combine = c) %dopar% {
                           
                           source("Functions_ad.R")
                           source("Goodness-of-fit-test_functions.R")
                           
                           # Create bootstrap data sample
                           bootstrapseed <- B*iseed + b
                           set.seed(bootstrapseed)
                           
                           beta <- parhat[1:parl]
                           eta <- parhat[(parl+1):totparl]
                           sd <- parhat[(totparl+1):(totparl+5)]
                           gamma <- gammaest
                           
                           mu <- c(0,0)
                           sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
                           err <- mvrnorm(n, mu=mu , Sigma=sigma)
                           err1 = err[,1]
                           err2 = err[,2]
                           
                           T.star = IYJtrans(M %*% beta + err1, sd[4])
                           C.star = IYJtrans(M %*% eta + err2, sd[5])
                           A.star <- rep(0, n)
                           for (i in 1:n) {
                             A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
                           }
                           
                           Y.star = pmin(T.star,C.star,A.star)
                           Delta.star = as.numeric(Y.star==T.star)
                           Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
                           data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
                           
                           initd1 <- parhat
                           parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                                                 ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
                           
                           Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
                           Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
                           
                           # For each observed time y, compute F_K(k; gamma, delta)
                           FY.star <- rep(0, n)
                           
                           for (i in 1:n) {
                             
                             # Vector of arguments of first term and second term
                             arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
                             arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
                             
                             # Evaluate them in standard normal CDF, and compute the average
                             FY.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                                                 pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
                             
                             # Take into account the administrative censoring
                             FY.star[i] <- FY.star[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]) + 
                               cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]
                           }
                           
                           FY.star <- FY.star[order(FY.star)]
                           
                           w.star = rep(0,n)
                           w.star[1] = FY.star[1]
                           w.star[2:n] = FY.star[2:n]-FY.star[1:(n-1)]
                           
                           ord <- 1:n
                           EFY.star <- (ord-1)/n
                           
                           n*sum((EFY.star - FY.star)^2 * w.star)
                         }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}

# This function is exactly like GOF_test, but instead of bootstrap sampling A
# from its estimated distribution, we sample from its true distribution.
# --> This function is no longer used and will no longer be updated.
GOF_test_oracle <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # structure of 'initd':
  #
  # [1:4] beta (4 params) = First 4 params of parhat1
  # [5:8] eta (4 params) = Next 4 params of parhat1
  # [9]   sigma1 = parhat1[9]
  # [10]  sigma2 = parhat1[10]
  # [11]  rho = 0
  # [12]  theta_1 = parhat1[11]
  # [13]  theta_2 = parhat1[12]
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1 or less.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    surv_km[["time"]] <- c(surv_km[["time"]][1:matched_until], surv_km[["time"]][matched_until:length(surv_km[["time"]])])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star = runif(n, 0, 8)
    # A.star <- rep(0, n)
    # for (i in 1:n) {
    #   A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    # }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1)
    cdf_km.star <- 1 - surv_km.star[["surv"]]
    
    # Same issue as before
    while (length(cdf_km.star) < length(FK.star)) {
      
      Y.star.sorted <- sort(Y.star)
      matched_until <- 0
      for (i in 1:length(surv_km.star[["time"]])) {
        if (Y.star.sorted[i] == surv_km.star[["time"]][i]) {
          matched_until <- i
        }
      }
      
      surv_km.star[["time"]] <- c(surv_km.star[["time"]][1:matched_until], surv_km.star[["time"]][matched_until:length(surv_km.star[["time"]])])
      cdf_km.star <- c(cdf_km.star[1:matched_until], cdf_km.star[matched_until:length(cdf_km.star)])
    }
    
    TCMb_vector[b] <- n*sum((cdf_km.star - FK.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}






