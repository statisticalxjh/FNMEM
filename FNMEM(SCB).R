rm(list=ls(all=TRUE))

library(MASS)
library(nlme)
library(locpol)

#setwd("C:/Users/Administrator/Desktop/functional nonlinear mixed effects models/Rcode")
source("sigma.R")
source("Kernelest.R")
source("CV.R")
source("Kh.R")
source("SCB.R")

#################### Speed up by mix programming with C++ ####################
#library(Rcpp)
#library(RcppArmadillo)
#sourceCpp("SCB.cpp")


set.seed(20221207)
#################### set parameters ####################
p <- 2                   	    		#dimension of fixed effects
q <- 2						#dimension of random effects
t.start <- 0					#beginning of grid point
t.end <- 1						#end of grid point
M <- 25						#number of grid points
n <- 50               	   		    	#number of subjects
n.rep <- 1                	    		#number of repeated trials
G <- 500						#number of bootstrap
mu.X <- rep(0, p)              		#mean of design matrix
rho <- 0.3					
Sigma.X <- sigma2(rho, p)			#variance of design matrix
alpha <- c(0.05, 0.01)					#confidence level
ni <- 5						#number of time points

f <- function(x1, x2, phi1, phi2){		#nonlinear function
  y <-  exp(x1*phi1 + x2*phi2)
  return(y)
}

fdbeta <- function(x1, x2, phi1, phi2){	#partial derivative of fixed effects
  y1 <- x1*exp(x1*phi1)
  y2 <- x2*exp(x2*phi2)
  y <- rbind(y1, y2)
  return(y)
}

fdb <- function(x1, x2, phi1, phi2){	#partial derivative of random effects
  y1 <- x1*exp(x1*phi1)
  y2 <- x2*exp(x2*phi2)
  y <- rbind(y1, y2)
  return(y)
}

#################### data generation ####################
t <- seq(t.start, t.end, length.out=M)		#gride points
CP <- matrix(0, p, length(alpha))  		#cover probability

options(warn=-1)
time.begin<-proc.time()
n.rep.obs <- 0
for (k in 1:n.rep){					#beginning of replications
  n.rep.obs <- n.rep.obs+1
  print(n.rep.obs)
  
  beta <- rbind(t^2, (1-t)^2)/2			#fixed effects
  Ni <- c() 								#time points
  Seed <- Y <- X <- NULL
  der.beta.hat <- beta.hat.smooth <- beta.hat <- matrix(0, p, M)
  b.hat.smooth <- b.hat <- b <- matrix(0, n*q, M)
  Phi.hat <- matrix(0, q*M, q)
  sigma.hat <- matrix(0, M, 1)
  
  for (i in 1:n){
    #ni <- 4+rpois(1, 1)					#number of time points
    Ni[i] <- ni
    Seed <- c(Seed, rep(i, ni))   #group
    Xi <- mvrnorm(ni, mu.X, Sigma.X)
    X <- rbind(X, Xi)						#covariates
    Yi <- NULL
    for (m in 1:M){
      bm <- mvrnorm(1, rep(0, q), 0.1*sigma2(rho, q))*sin(2*pi*t[m]) + mvrnorm(1, rep(0, q), 0.2*sigma2(rho, q))*cos(2*pi*t[m])		#random effects
      b[(q*i - 1):(q*i), m] <- bm		
      phi1 <- beta[1, m] + bm[1]
      phi2 <- beta[2, m] + bm[2]	
      Yim<- f(Xi[, 1], Xi[, 2], phi1, phi2) + rnorm(1, 0, 0.1)
      Yi <- cbind(Yi, t(t(Yim)))
    }
    Y <- rbind(Y, Yi)						#reponses
  }
  
  for (m in 1:M){
    Ym <- Y[, m]
    mydata <- cbind(X, Ym, Seed)
    mydata <- groupedData(Ym ~ V1 + V2 | Seed, as.data.frame(mydata))
    
    

    
    #################### estimator of fixed effects ####################
    repeat{
      init<-runif(1,0,beta[1,m])
      fm <- try(nlme(Ym ~ f(V1, V2, phi1, phi2),  	#nonlinear mixed model
                     data = mydata,
                     fixed = phi1 + phi2 ~ 1,
                     random =  phi1 + phi2  ~ 1,
                     start = c(phi1 = init, phi2 = beta[2, m]),
                     method = "REML",
                     control=nlmeControl(returnObject = TRUE)
      ), silent=TRUE)
      if(!('try-error'%in%class(fm))) break
    }
    beta.hat[, m] <- fixed.effects(fm)			#estimators of fixed effects
    temp <- random.effects(fm)
    temp1 <- as.numeric(rownames(temp)) 
    temp1 <- sort(temp1, index.return=T)$ix 
    temp <- temp[temp1,]
    b.hat[, m] <- array(t(temp))  #estimators of random effects
    temp <- matrix(as.numeric(VarCorr(fm)[8]), q, q)*as.numeric(VarCorr(fm)[4])*as.numeric(VarCorr(fm)[5])
    diag(temp) <- as.numeric(VarCorr(fm)[1:q])
    Phi.hat[(q*(m - 1) + 1):(q*m), ] <- temp   #estimators of variances of random effects
    sigma.hat[m, 1] <- as.numeric(VarCorr(fm)[3]) #estimators of variances of errors
  }
  
  
  ########## kernel smooth of fixed effects ##########
  h1 <- CV(beta.hat, t)  #bandwidth
  h1 <- max(h1/2, (t[2]-t[1])*1.05)
  for (i in 1:p){
    beta.hat.smooth[i,] <- kernelest(beta.hat[i,], t, t, h1)	#smoothed estimators of fixed effects
    d <- data.frame(t)
    d$y <- beta.hat.smooth[i,]
    kernel.beta <- locpol(y ~ t, d, deg=2, xeval=t)
    der.beta.hat[i, ]<- kernel.beta$lpFit[, 4]  #second order derivatives of fixed effects
  }

  
  #################### simultaneous confidence bands ####################
  
  ########## unbias ##########
  beta.hat.unbias <- beta.hat.smooth
  mu2K <- 0.2
  for (m in 1:M){
    for (i in 1:p){
      bias <- 0.5*der.beta.hat[i,m]*mu2K*h1^2
      beta.hat.unbias[i, m] <- beta.hat.smooth[i, m] - bias #unbias estimators of fixed effects
    }
  }
  
  ########## C ##########
  allC <- SCB(n, M, p, q, ni, G, h1, t, X, Y, Seed, beta.hat, b.hat, Phi.hat, sigma.hat)
  C <- matrix(0, p, length(alpha))
  for (i in 1:p){
    temp <- allC[i, ]
    temp <- sort(temp)
    for (j in 1:length(alpha)){
      C[i, j] <- temp[G*(1-alpha[j])]
    }
  }

  ########## cover probability ##########
  for (i in 1:p){
    for (j in 1:length(alpha)){
      temp <- 0
      for (m in 1:M){		
        L.beta.hat <- beta.hat.unbias[i, m] - C[i, j] #lower bound
        U.beta.hat <- beta.hat.unbias[i, m] + C[i, j] #upper bound 
        if ((L.beta.hat<=beta[i, m]) && (U.beta.hat>=beta[i, m])) {temp <- temp + 1}
      }
      CP[i, j] <- CP[i, j] + (temp==M)/n.rep
    }
  }
  
}				# end of replications

colnames(CP) <- c('alpha=0.05','alpha=0.01')
rownames(CP) <- c('beta1','beta2')
CP

time.end<-proc.time()
time.end-time.begin

