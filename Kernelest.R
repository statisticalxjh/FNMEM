#################### Kernel estimator ####################
kernelest <- function(input, t.in, t.out, h){
	n.in <- length(t.in)
  n.out <- length(t.out)
	output <- c()
	for (i in 1:n.out){
		K <- c()
		for (j in 1:n.in){
			delta <- (t.in[j] - t.out[i])
			K[j] <-  Kh(delta, h)
		}
		output[i] <- sum(K*input)/sum(K)
	}
	return(output)
}