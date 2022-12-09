#################### stochastic process G ####################
spG <- function(X, Y, Ni, Seed, beta.hat, b.hat, Phi.hat, sigma.hat, s, tau){

sp <- 0
K <- rep(0, M)
for(m in 1:M){	
	if (abs(t[m]-t[s])<h1) {
		Sxx <- Sxy <- 0
		Kh1 <- Kh(t[m]-t[s], h1)
		K[m] <- Kh1
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
			Sxy <- Sxy + XXi%*%Sigmai%*%(wi - t(XXi)%*%beta.hat[, m])*tau[i]
		}
	sp <- sp + Kh1*solve(Sxx)%*%Sxy
	}	
}

sp <- sp/sum(K)

return(sp)

}

#################### simultaneous confidence bands ####################
SCB <- function(n, M, p, q, ni, G, h1, t, X, Y, Seed, beta.hat, b.hat, Phi.hat, sigma.hat){

allC <- matrix(0, p, G)
for (g in 1:G){
	#print(g)
	tau <- rnorm(n, 0, 1)
	sp <- matrix(0, p, M)
	for (s in 1:M){
		sp[, s] <- spG(X, Y, Ni, Seed, beta.hat, b.hat, Phi.hat, sigma.hat, s, tau)
	}
	
	for (i in 1:p){
		allC[i, g] <- max(abs(sp[i, ]))
	}
}

return(allC)

}