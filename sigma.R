######################################## Sigma of model 1 ########################################
sigma1<-function(rho,p){
sigma<-matrix(rho,p,p)
for (i in 1:p){
	sigma[i,i]<-1
}
return(sigma)
}

######################################## Sigma of model 2 ########################################
sigma2<-function(rho,p){
  sigma <- rho^abs(outer(1:p,1:p,'-')) 
return(sigma)
}
