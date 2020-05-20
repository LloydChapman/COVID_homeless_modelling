calc_distance <- function(D_S,D_S_star,D_T,D_T_star){
	sre1 <- sqrt(sum((D_S - D_S_star)^2))
	sre2 <- sqrt(sum((D_T - D_T_star)^2))
	return(c(sre1,sre2))
}

# Perturbation kernel 
rK <- function(mean, sigma, lm.low, lm.upp){   
	return(rtmvnorm(1,mean=mean, sigma=sigma, lower=lm.low, upper=lm.upp)) 
}

#  Heaviside function: H(x)=1 if x>0
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero <- function(par, lm.low, lm.upp){
	prod(sapply(1:length(par), function(a) H(par[a]-lm.low[a])* H(lm.upp[a]-par[a])))
}
		    
Norm.Eucl.dist <- function(p1, p2, lm.low, lm.upp){
	sqrt(sum(((p1-p2)/(lm.upp-lm.low))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours <- function(N, M, theta, Theta, lm.low, lm.upp){
	dist <- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,]), lm.low, lm.upp))
	temp <- data.frame(no=seq(1,N), dist)
	temp <- temp[order(temp$dist),]
	if (length(theta)==1){
	  sigma <- var(Theta[temp$no[1:(M+1)],])
	}
	else{
	  sigma <- cov(Theta[temp$no[1:(M+1)],])
	}
	return(sigma)
}

# Effective sample size
ESS <- function(w){1/sum(w^2)}
  
  
