# install.packages(c("tmvtnorm","actuar"))
library(tmvtnorm)
library(actuar)
library(splines)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(gsubfn)

setwd("~/Dropbox/Homeless/Code") 
source("ABC_SMC_MNN_functions.R")
source("COVID_homeless_calibration2.R")
source("COVID_homeless_functions.R")

# Load CCMS data from MSC South outbreak
CCMS_data <- read.csv("../Data/CCMS_data.csv",stringsAsFactors = F)
names(CCMS_data)[1] <- "Date"
CCMS_data$Date <- as.Date(CCMS_data$Date,format = "%d-%b")
# Remove empty rows from after Apr 10
CCMS_data <- CCMS_data[CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(CCMS_data))

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
SF_case_data <- aggregate(Case.Count ~ Date,SF_data,sum)

T_sim <- 13
N_pop <- 415

# sim_data <- read.csv("sim_data.csv")
# 
# # Length of observed outbreak (weeks)
# n_obs <- nrow(sim_data)

# True (simulated) data
# D_T <- sim_data$new_infections
D_S <- c(2,1,2,5)
D_T <- c(29,33) #70 - sum(D_S) # number of positives out of 134 tested during mass testing on 4/8-4/9
N_tested <- rep((144-10)/2,2)

# Number of particles
N <- 1000

# Number of neighbours for covariance matrix calculations
M <- 50

# Epsilon tolerance values for number PCR positive in early testing of symptomatics
epsilon_S <- round(seq(sqrt(sum((2*D_S)^2)),sqrt(sum((0.3*D_S)^2)),length.out = 10))

# Epsilon tolerance values for number PCR positive in mass testing
CI1 <- binom.test(D_T[1],N_tested[1])$conf.int
CI2 <- binom.test(D_T[2],N_tested[2])$conf.int
tol1 <- diff(CI1*N_tested[1])/2
tol2 <- diff(CI2*N_tested[2])/2
epsilon_T <- round(seq(sqrt(sum((2*D_T)^2)),sqrt(tol1^2+tol2^2),length.out = 10))

# Number of generations
G <- length(epsilon_T)

# Number of simulations for each parameter set
n <- 1

#  Lower and upper boundaries for priors
lm.low <- 0
lm.upp <- 0.1
npar <- length(lm.low)

# Empty matrices to store results (population x number of model parameters)
res.old <- matrix(nrow=N,ncol=npar)
res.new <- matrix(nrow=N,ncol=npar)

# Empty vectors to store weights
w.old <- matrix(nrow=N,ncol=npar)
w.new <- matrix(nrow=N,ncol=npar)

# # Empty list to store output of each accepted simulation
# res_list <- vector("list",N)

# Empty vectors and matrices to store outputs for each simulation
new_infections_mat <- matrix(nrow=N,ncol=T_sim) # daily numbers of new infections (by symptom onset)
PCRpos_sx_testing_mat <- matrix(nrow=N,ncol=length(D_S)) # daily number PCR positive during symptomatic testing 4/4-4/7
PCRpos_mat <- matrix(nrow=N,ncol=length(D_T)) # number PCR positive in mass testing
state_arr <- array(dim=c(N_pop,T_sim,N)) # states of individuals over time
sim_pop_list <- vector("list",N) # individual-level information from end of simulation

# Empty vector for storing acceptance rates
acc_rate <- numeric(G)

for(g in 1:G){  

	#Initiate counter of accepted particles
	i <- 1
	# Initiate counter of proposed particles
	k <- 1
	while(i <= N){ # While the number of accepted particles is less than N
   		if(g==1){
    	# Sample from prior distribution
			beta_star <- runif(1, min=lm.low[1], max=lm.upp[1]) 
		} else {
			#  Select particle from previous generation
			p <- sample(seq(1,N),1,prob=w.old)		
			sigma <- Sigma[[p]]
			par <- rK(res.old[p,],sigma,lm.low,lm.upp)
			beta_star <- par[1]
		}
      	#  Test if prior non zero
      	if(prior.non.zero(beta_star,lm.low,lm.upp)) {
    			# Set number of accepted simulations to zero
    			m <- 0
    			distance <- matrix(nrow=n,ncol=2)
    			for(j in 1:n){
    				res <- COVID_homeless_model(beta_star,CCMS_data,SF_case_data)     
    				D_S_Star <- res$PCRpos_sx_testing
    				D_T_star <- res$PCRpos
    				# Calculate distances 
    				distance[j,] <- calc_distance(D_S,D_S_Star,D_T,D_T_star) # [ ] - UPDATE to save this in output
    				if((distance[j,1] <= epsilon_S[g]) & (distance[j,2] <= epsilon_T[g])){ # If distances are less than tolerance
    					m <- m+1
    				}
    			}	
    			if (m>0){
    				# Store results
    				res.new[i,] <- beta_star
    				# res_list[[i]] <- res
    				new_infections_mat[i,] <- res$new_infections
    				PCRpos_sx_testing_mat[i,] <- res$PCRpos_sx_testing
    				PCRpos_mat[i,] <- res$PCRpos
    				state_arr[,,i] <- res$state
    				sim_pop_list[[i]] <- res$sim_pop
      			# Calculate weights
      			w1 <- prod(sapply(1:npar, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])))
				if(g==1){
					w2 <- 1
				} else {
					w2 <- sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
				}
      			w.new[i] <- (m/n)*w1/w2
      			print(paste0('Generation: ', g, ", particle: ", i))
      			# Update counter of accepted particles
      			i <- i+1
      			}
      	}
	  # Update counter of proposed particles
	  k <- k+1
    	}
    	Sigma <- list(NA, N)
    for(p in 1:N){
    	Sigma[[p]] <- getSigmaNeighbours(N, M, res.new[p,], res.new, lm.low, lm.upp) 
    }
    	res.old <- res.new
	w.old <- w.new/sum(w.new)
	acc_rate[g] <- i/k
  cat("Acceptance rate for ", g,"th generation = ", round(acc_rate[g],2), "\n", sep = "")

 	write.csv(res.new, file = paste("results_ABC_SMC_MNN_gen_",g,"_9.csv",sep=""), row.names=FALSE)
}

# Calculate ESS
ess <- ESS(w.old)

save(new_infections_mat,PCRpos_sx_testing_mat,PCRpos_mat,state_arr,sim_pop_list,ess,file = "sim_output_8.RData")
