  #################### leave-one-out CV for choosing h ####################
CV <- function(inputall, t){

n <- nrow(inputall)
n.input <- length(t)
temp <- c()
for (i in 1:(n.input-1)) {temp[i] <- t[i + 1] - t[i]}
h.min <- max(temp)*1.05
h.max <- (t[n.input] - t[1])/5
n.cho <- 50
h.cho <- seq(h.min, h.max, length.out = n.cho)
bias <- matrix(0, n, n.cho)

for (nn in 1:n){
	input <- inputall[nn,]
	for (i in 1:n.cho){
		bias2 <- 0
		for (j in 1:n.input){
			bias2 <- bias2 + (kernelest(input[-j], t[-j], t[j], h.cho[i]) - input[j])^2
		}
		bias[nn, i] <- bias2
	}
}
biasall <- apply(bias, 2, sum)
index <- which.min(biasall)
return(h.cho[index])
}

#################### * folds CV for choosing h ####################
CV2 <- function(inputall, t, n.fold){
  
  n <- nrow(inputall)
  n.input <- length(t)
  temp <- c()
  for (i in 1:(n.input-1)) {temp[i] <- t[i + 1] - t[i]}
  h.min <- max(temp)*1.05
  h.max <- (t[n.input] - t[1])/5
  n.cho <- 50
  h.cho <- seq(h.min, h.max, length.out = n.cho)
  bias <- matrix(0, n, n.cho)
  index <- sample(n.input)
  index <- 1:n.input
  n.each <- floor(n.input/n.fold)
  for (nn in 1:n){
    input <- inputall[nn,]
    for (i in 1:n.cho){
      bias2 <- 0
      for (j in 1:n.fold){      
        index.tes <- index[(n.each*(j-1)+1):(n.each*j)]
        bias2 <- bias2 + sum((kernelest(inputall[nn,-index.tes], t[-index.tes], t[index.tes], h.cho[i]) - inputall[nn,index.tes])^2)
      }
      bias[nn, i] <- bias2
    }
  }
  biasall <- apply(bias, 2, sum)
  index <- which.min(biasall)
  return(h.cho[index])
}