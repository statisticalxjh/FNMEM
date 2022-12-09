#################### calculate Sigmahat ####################
calcu.sm2 <- function(X, Y, Ni, Seed, beta.hat, b.hat, Phi.hat, sigma.hat, m){
  Sxx <- 0
  score <- matrix(0, p, n)
  for (i in 1:n){
    XX <- X[Seed==i, ]
    YY <- t(t(Y[Seed==i, m]))
    phi1.hat <- beta.hat[1, m] + b.hat[(2*i - 1), m]
    phi2.hat <- beta.hat[2, m] + b.hat[(2*i), m]
    XXi <- fdbeta(XX[, 1], XX[, 2], phi1.hat, phi2.hat)
    ZZi <- fdb(XX[, 1], XX[, 2], phi1.hat, phi2.hat)
    wi <- YY - f(XX[, 1], XX[, 2], phi1.hat, phi2.hat) + t(XXi)%*%beta.hat[, m] + t(ZZi)%*%b.hat[(2*i - 1):(2*i), m]
    ddelta.inv <- Phi.hat[(q*(m - 1) + 1):(q*m), ]/sigma.hat[m, 1]
    Sigmai <- diag(1, Ni[i]) +  t(ZZi)%*%ddelta.inv%*%ZZi  
    Sigmai <- solve(Sigmai)
    Sxx <- Sxx + XXi%*%Sigmai%*%t(XXi)          
    score[, i] <- XXi%*%Sigmai%*%(wi - t(XXi)%*%beta.hat[, m])
  }
  Hn.inv <- solve(Sxx)
  Sigma.hat <- Hn.inv%*%score%*%t(score)%*%Hn.inv
  return(Sigma.hat)
}

#################### calculate Sigmahat (bootstrap) ####################
calcu.sm.boot2 <- function(X, Y, Ni, Seed, beta.hat, beta.hat.star, b.hat, Phi.hat, sigma.hat, m, tau){
  Sxx <- 0
  score.boot <- score <- matrix(0, p, n)
  for (i in 1:n){
    XX <- X[Seed==i, ]
    YY <- t(t(Y[Seed==i, m]))
    phi1.hat <- beta.hat[1, m] + b.hat[(2*i - 1), m]
    phi2.hat <- beta.hat[2, m] + b.hat[(2*i), m]
    XXi <- fdbeta(XX[, 1], XX[, 2], phi1.hat, phi2.hat)
    ZZi <- fdb(XX[, 1], XX[, 2], phi1.hat, phi2.hat)
    wi <- YY - f(XX[, 1], XX[, 2], phi1.hat, phi2.hat) + t(XXi)%*%beta.hat[, m] + t(ZZi)%*%b.hat[(2*i - 1):(2*i), m]
    ddelta.inv <- Phi.hat[(q*(m - 1) + 1):(q*m), ]/sigma.hat[m, 1]
    Sigmai <- diag(1, Ni[i]) +  t(ZZi)%*%ddelta.inv%*%ZZi  
    Sigmai <- solve(Sigmai)
    Sxx <- Sxx + XXi%*%Sigmai%*%t(XXi)
    score[, i] <- XXi%*%Sigmai%*%(wi - t(XXi)%*%beta.hat[, m])*tau[i] 
    score.boot[, i] <- XXi%*%Sigmai%*%(wi - t(XXi)%*%beta.hat.star[, m])*tau[i]
  }
  Hn.inv <- solve(Sxx)
  Sigma.hat <- Hn.inv%*%score%*%t(score)%*%Hn.inv
  score.boot <- apply(score.boot, 1, sum)
  return(list(Hn.inv, score.boot, Sigma.hat))
}

##################################################
# hyp test
##################################################
test.p <- function(X, Y, R, b0, h1, Ni, Seed, beta.hat, beta.hat.unbias, b.hat, Phi.hat, sigma.hat ){
  
  ########## calculate Sn ##########
  sn <- 0
  for (m in 1:M){
    Sigma.hat <- calcu.sm2(X, Y, Ni, Seed, beta.hat, b.hat, Phi.hat, sigma.hat, m)
    d <- R%*%beta.hat[, m] - b0[m] #NMEM
    sn <- sn + d^2*solve(R%*%Sigma.hat%*%t(R))/M
  }
  
  ########## fit null ##########
  beta.hat.null <- beta.hat.star <- matrix(0, p, M)  
  for (m in 1:M){
    Ym <- Y[, m]
    mydata <- cbind(X, Ym, Seed)
    mydata <- groupedData(Ym ~ V1 + V2 | Seed, as.data.frame(mydata))
    repeat{
      init <- runif(1,0,0.5)
      fm <- try(nlme(Ym ~ f(V1, V2, phi1, phi2),  	#nonlinear mixed model
                     data = mydata,
                     fixed = phi2  ~ 1,
                     random =  phi1 + phi2  ~ 1,
                     start = c( phi2 = init),
                     method = "REML",
                     control=nlmeControl( returnObject = TRUE)
      ), silent=TRUE)
      if(!('try-error'%in%class(fm))) break
    }
    beta.hat.star[2, m] <- fixed.effects(fm)			#estimators of fixed effects   
  }
  
#  for (i in 1:p){
#    beta.hat.null[i,] <- kernelest(beta.hat.star[i,], t, t, h1)  #smoothed estimators of fixed effects
#  }
  
  p.value <- 0
  for (g in 1:G){
    #print(g)	
    ########## calculate Sng ##########
    sng <- 0
    tau <- rnorm(n, 0, 1)
    #tau <- rep(1,n)
    for (m in 1:M){      
      temp <- calcu.sm.boot2(X, Y, Ni, Seed, beta.hat, beta.hat.star, b.hat, Phi.hat, sigma.hat, m, tau) #NMEM
      Hn.inv <- temp[[1]]
      score.boot <- temp[[2]]
      Sigma.hat <- temp[[3]]
      temp <- R%*%Hn.inv%*%score.boot
      sng <- sng + temp^2*solve(R%*%Sigma.hat%*%t(R))/M
    }
    ########## calculate p ##########
    if (sng>=sn){p.value <- p.value+1/G}
    
  }
  return(p.value)
}
