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
# Replace NA's with 0's
CCMS_data[is.na(CCMS_data)] <- 0

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
SF_case_data <- aggregate(Case.Count ~ Date,SF_data,sum)

# Set number of residents and staff in shelter and duration of simulation
N_res <- 350 # Set such that ~255 unique individuals are present in shelter at some point between 3/29 and 4/10  # CCMS_data$Total_Census[1]
N_staff <- 65
N_pop <- N_res + N_staff
start_date <- min(CCMS_data$Date)
end_date <- max(CCMS_data$Date)
T_sim <- as.integer(end_date - start_date + 1) #13 # testing 5 days after 1st cases confirmed on 4/5 + 2 days presymptomatic infectiousness (assumed) + 6 days of symptoms of 1st case  # 30 # Days; Total of 1 month

# Set weights for presence of residents and staff in shelter
w <- c(rep(1,N_res),rep(1/2,N_staff))

# Set natural history parameters
source("set_nat_hist_pars.R")

# Flag for whether to count hospitalisations and deaths
hospitalisation <- F # false as data not available for MSC South outbreak

# Set PCR test parameters
source("set_PCR_test_pars.R")

# PCR testing frequency
testing_dates <- as.Date(c("4/8/2020","4/9/2020"),format="%m/%d/%Y") # testing dates
testing_days <- as.integer(testing_dates - start_date + 1) # days on which testing occurred
N_tested <- rep((144-10)/2,2) # number tested during mass testing = total tested from press release - number tested during first 4 days = 144 - 10, # number tested on each testing day #c(72,72) # assume all testing happened during mass testing and was equally split between 4/8 and 4/9
sx_testing_dates <- as.Date(c("4/4/2020","4/5/2020","4/6/2020","4/7/2020"),format="%m/%d/%Y")
sx_testing_days <- as.integer(sx_testing_dates - start_date + 1)

# Initialise variables
Number <- 1:N_pop
Alive <- rep(1,N_pop)
Resident <- rep(0,N_pop)
Resident[1:N_res] <- 1
res_present0 <- sample(1:N_res,CCMS_data$Total_Census[1])
res_absent0 <- setdiff(1:N_res,res_present0)
staff_present0 <- (N_res+1):N_pop
Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
v <- round(N_res*CCMS_data[1,2:ncol(CCMS_data)]/CCMS_data[1,2],0)
Hi_Risk_60_only_present0 <- sample(res_present0,CCMS_data$Hi_Risk_60_only[1])
Hi_Risk_60_only_absent0 <- sample(res_absent0,v$Hi_Risk_60_only - CCMS_data$Hi_Risk_60_only[1])
Risk[c(Hi_Risk_60_only_present0,Hi_Risk_60_only_absent0)] <- 2
Hi_Risk_Dx_only_present0 <- sample(setdiff(res_present0,Hi_Risk_60_only_present0),CCMS_data$Hi_Risk_Dx_only[1])
Hi_Risk_Dx_only_absent0 <- sample(setdiff(res_absent0,Hi_Risk_60_only_absent0),v$Hi_Risk_Dx_only - CCMS_data$Hi_Risk_Dx_only)
Risk[c(Hi_Risk_Dx_only_present0,Hi_Risk_Dx_only_absent0)] <- 3
Hi_Risk_Both_Age_Dx_present0 <- sample(setdiff(res_present0,c(Hi_Risk_60_only_present0,Hi_Risk_Dx_only_present0)),CCMS_data$Hi_Risk_Both_Age_Dx[1])
Hi_Risk_Both_Age_Dx_absent0 <- sample(setdiff(res_absent0,c(Hi_Risk_60_only_absent0,Hi_Risk_Dx_only_absent0)),v$Hi_Risk_Both_Age_Dx-CCMS_data$Hi_Risk_Both_Age_Dx[1])
Risk[c(Hi_Risk_Both_Age_Dx_present0,Hi_Risk_Both_Age_Dx_absent0)] <- 4
Age <- rep(NA,N_pop)
Age[Risk %in% c(1,3)] <- sample(x=seq(18,59), size=sum(Risk %in% c(1,3)), replace=TRUE)
Age[Risk %in% c(2,4)] <- sample(x=seq(60,69), size=sum(Risk %in% c(2,4)), replace=TRUE) # [ ] - CHECK oldest age (have assumed 69 for now)
TrueState <- rep(1,N_pop)
# "Index" cases:
# 1st case with sx onset on 3/31
i_s_p0 <- sample(res_present0,1) # draw at random from residents who are present
TrueState[i_s_p0] <- 4 # assume initially in severe presymptomatic state
# 2nd case with sx onset on 4/2, assume initially in exposed state
e0 <- sample(setdiff(res_present0,i_s_p0),1) # draw at random from remaining residents who are present
e0ind <- rep(F,N_pop)
e0ind[e0] <- T
TrueState[e0] <- 2 # assume initially in latent state
DayTrueState <- rep(0,N_pop)
DayTrueState[e0] <- 1 # assume 2nd index case has been in latent state for 1 day
WaitingTime <- rep(NA,N_pop)
WaitingTime[i_s_p0] <- 2 # assume 1st index case has presymptomatic duration of 2 days (~mean of presymptomatic duration distribution) # rnbinom(length(i_s_p0),r_p,p_p)+1
WaitingTime[e0] <- 3 # assume 2nd index case has latent duration of 3 days (mean of latent duration distribution)
DaysSinceInfctn <- rep(NA,N_pop) # days since infection
DaysSinceInfctn[i_s_p0] <- rbinom(length(i_s_p0),r_E,p_E)+1 # latent period of 1st index case prior to start of simulation
DaysSinceInfctn[e0] <- DayTrueState[e0]
DaysSinceInfctsnss <- rep(NA,N_pop) # days since start of infectiousness (i.e. start of presymptomatic infectious stage)
DaysSinceInfctsnss[i_s_p0] <- 0
DaysPCRpos <- rep(NA,N_pop) # duration of PCR positivity (viraemia)

# sim_data <- read.csv("sim_data.csv")
# 
# # Length of observed outbreak (weeks)
# n_obs <- nrow(sim_data)

# True (simulated) data
# D_T <- sim_data$new_infections
D_S <- c(2,1,2,5)
D_T <- c(29,33) #70 - sum(D_S) # number of positives out of 134 tested during mass testing on 4/8-4/9

# Number of particles
N <- 10 #1000

# Number of neighbours for covariance matrix calculations
M <- 1 #50

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
infections_mat <- matrix(nrow=N,ncol=T_sim) # daily numbers of new infections
cases_mat <- matrix(nrow=N,ncol=T_sim) # daily numbers of new clinical cases (by symptom onset)
PCRpos_sx_testing_mat <- matrix(nrow=N,ncol=length(D_S)) # daily number PCR positive during symptomatic testing 4/4-4/7
PCRpos_mat <- matrix(nrow=N,ncol=length(D_T)) # number PCR positive in mass testing
state_arr <- array(dim=c(N_res+N_staff,T_sim,N)) # states of individuals over time
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
    				res <- COVID_homeless_model(N_res,N_staff,N_pop,T_sim,w,beta_star,epsilon,r_E,p_E,p_s,h,r_p,p_p,alpha,r_sx,p_sx,
    				                            p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,
    				                            hospitalisation,fit,fit_extrap,spec,testing_days,N_tested,sx_testing_days,
    				                            CCMS_data,Number,Alive,Resident,Present,Risk,Age,TrueState,DayTrueState,
    				                            WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos)
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
    				infections_mat[i,] <- res$infections
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

 	write.csv(res.new, file = paste("results_ABC_SMC_MNN_gen_",g,"_10.csv",sep=""), row.names=FALSE)
}

# Calculate ESS
ess <- ESS(w.old)

save(infections_mat,PCRpos_sx_testing_mat,PCRpos_mat,state_arr,sim_pop_list,ess,file = "sim_output_10.RData")
