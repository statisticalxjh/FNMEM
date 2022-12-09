#################### Kernel function ####################
Kh <- function(mu, h){
	output <- (1 - (mu/h)^2)*0.75*(abs(mu/h)<=1)/h
	return(output)
}

