#################### weighted nonlinear least sequare ####################
wnls <- function(Xi, Yi, ni, h2, beta.hat, s, eps, max.step){

bi <- matrix(0, q, 1)
delta <- 1000
ite <- 0

while ((delta>eps)&&(ite<max.step)){
	ite <- ite + 1
	Ssh2 <- 0
	SS <- 0
	for (j in 1:ni){
		for(m in 1:M){
		  if (abs(t[m]-t[s])<h2) {
			phi10 <- beta.hat[1, m] + bi[1, 1]
			phi20 <- beta.hat[2, m] + bi[2, 1]			
			A <- fdb(Xi[j, 1], Xi[j, 2], phi10, phi20)
			W <- Yi[j, m] - f(Xi[j, 1], Xi[j, 2], phi10, phi20) + t(A)%*%bi
			Z <- rbind(1, t[m]-t[s])
			Kh2 <- Kh(t[m]-t[s], h2)
			Ssh2 <- Ssh2 + Kh2*kronecker(A%*%t(A), Z%*%t(Z))
			SS <- SS + Kh2*kronecker(A, Z)%*%W
		  }
		}
	}
	bi.hat <- kronecker(diag(1, q),t(c(1, 0)))%*%ginv(Ssh2)%*%SS
	delta <- sum(abs(bi - bi.hat))
	bi <- bi.hat
}

return(bi)
}