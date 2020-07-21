ABC_SMC_MNN <- function(shelter,run_nm){

# Set data and initial conditions for shelter
set_data_and_init_condns(shelter)

# Number of particles
N <- 1000

# Number of neighbours for covariance matrix calculations
M <- 50

# Epsilon tolerance values for number PCR positive in early testing of symptomatics
epsilon_S <- round(seq(sqrt(sum((2*D_S)^2)),sqrt(sum((D_S/1.5)^2)),length.out = 10))

# Epsilon tolerance values for number PCR positive in mass testing
CI <- sapply(1:length(D_T),function(i) binom.test(D_T[i],N_tested[i])$conf.int)
tol <- diff(t(t(CI)*N_tested))/2
epsilon_T <- round(seq(sqrt(sum((2*D_T)^2)),sqrt(sum(tol^2)),length.out = 10))

# Reporting rate of symptomatic cases
sx_reporting_rate <- 1 #55/106 #36/45 # 36 symptomatic cases in DPH data out of 45 in Liz's 4/27 data summary = 80% reporting rate
# Epsilon tolerance values for number of symptom onsets
epsilon_C <- round(seq(sqrt(sum((2*D_C)^2)),sqrt(sum((1.3*D_C)^2)),length.out = 10)) # this is ignored below if D_C=NULL

# Number of generations
G <- length(epsilon_T)

# Number of simulations for each parameter set
n <- 1

# Number of parameters and parameter names
npar <- length(lm.low)
parnames <- c("beta","E0","T")
# parnames <- c("R0","E0","T")
# Indices of continuous and discrete parameters
cpar <- 1
dpar <- 2:3
ncpar <- length(cpar)

# Fixed variance for beta proposals (converted from fixed proposal variance for R0)
sigma_R0 <- 1
sigma <- calc_beta(sqrt(sigma_R0),w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)^2
# # Fixed variance for R0 proposals
# sigma <- 2

# Empty matrices to store results (population x number of model parameters)
res.old <- matrix(nrow=N,ncol=npar)
res.new <- matrix(nrow=N,ncol=npar)
colnames(res.old) <- parnames
colnames(res.new) <- parnames

# Empty vectors to store weights
w.old <- numeric(N)
w.new <- numeric(N)

# # Empty list to store output of each accepted simulation
# res_list <- vector("list",N)

# Empty objects to store outputs for each simulation
# infections_mat <- matrix(nrow=N,ncol=T_sim) # daily numbers of new infections
infections_list <- vector("list",N) # daily numbers of new infections
# cases_mat <- matrix(nrow=N,ncol=T_sim) # daily numbers of new clinical cases (by symptom onset)
cases_list <- vector("list",N) # daily numbers of new clinical cases (by symptom onset)
PCRpos_sx_testing_mat <- matrix(nrow=N,ncol=length(D_S)) # daily number PCR positive during symptomatic testing 4/4-4/7
PCRpos_mat <- matrix(nrow=N,ncol=length(D_T)) # number PCR positive in mass testing
cases_mat <- matrix(nrow=N,ncol=length(D_C)) # daily number of symptomatic cases from XX to XX
# state_arr <- array(dim=c(N_res+N_staff,T_sim,N)) # states of individuals over time
state_list <- vector("list",N) # states of individuals over time
sim_pop_list <- vector("list",N) # individual-level information from end of simulation

# Empty matrix for marginal probabilities of different combinations of discrete parameters
marg_prob <- matrix(nrow = lm.upp[dpar[1]]-lm.low[dpar[1]]+1, ncol = lm.upp[dpar[2]]-lm.low[dpar[2]]+1)

# Empty vector for storing acceptance rates
acc_rate <- numeric(G)

for(g in 1:G){
	# Initialise counter of accepted particles
	i <- 0
	# Initialise counter of proposed particles
	k <- 0
	while(i < N){ # While the number of accepted particles is less than N
	  # start_tot_time <- Sys.time()
   	if(g==1){
    	# Sample from prior distribution
   	  # start_time <- Sys.time()
   	  E0_star <- sample(lm.low[2]:lm.upp[2],1) #round(runif(1, min=lm.low[2], max=lm.upp[2])) #
   	  T_sim_star <- sample(lm.low[3]:lm.upp[3],1) #round(runif(1, min=lm.low[3], max=lm.upp[3])) #
   	  beta_star <- runif(1, min=lm.low[1], max=lm.upp[1])
   	  # R0_star <- runif(1, min=lm.low[1], max=lm.upp[1])
			# beta_star <- calc_beta(R0_star,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
			# end_time <- Sys.time()
			# print(end_time - start_time)
		} else {
			#  Select particle from previous generation
		  # start_time <- Sys.time()
	    E0_star <- sample(lm.low[2]:lm.upp[2],1) #round(par[2]) #
	    T_sim_star <- sample(lm.low[3]:lm.upp[3],1) #round(par[3]) #
	    # p <- sample(seq(1,N),1,prob=w.old)
	    idx <- which(apply(res.old[,dpar],1,function(x) all(x==c(E0_star,T_sim_star))))
	    if (length(idx)!=0){ # if there are any remaining particles with these values of the discrete parameters
	      p <- resample(idx,1,prob=w.old[idx])
	      # sigma <- Sigma[[p]]
	      par <- rK(res.old[p,cpar],sigma,lm.low,lm.upp)
	      beta_star <- par[1]
	      # R0_star <- par[1]
	      # beta_star <- calc_beta(R0_star,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
	    }
			# end_time <- Sys.time()
			# print(end_time - start_time)
		}
	  # Test if there are any remaining particles with the proposed values of the discrete parameters
	  if (g==1 || (g>1 & length(idx)!=0)){
    	#  Test if prior non zero
	    if(prior.non.zero(c(beta_star,E0_star,T_sim_star),lm.low,lm.upp)) {
	    # if(prior.non.zero(c(R0_star,E0_star,T_sim_star),lm.low,lm.upp)) {
  			# Set number of accepted simulations to zero
  			m <- 0
  			if(is.null(D_C)){ # 2 columns if no symptom onset data available, else 3 columns
  			  distance <- matrix(nrow=n,ncol=2)
  			} else {
  			  distance <- matrix(nrow=n,ncol=3)
  			}
  			for(j in 1:n){
  			  if (!is.null(CCMS_data)){
  			    # Backfill CCMS data with numbers in different risk categories on Mar 29
  			    CCMS_data0 <- cbind(Date=CCMS_data$Date[1] - seq(T_sim_star-nrow(CCMS_data),1,by=-1),data.frame(CCMS_data[1,2:ncol(CCMS_data)][rep(1,T_sim_star-nrow(CCMS_data)),]))
  			    tmp <- rbind(CCMS_data0,CCMS_data)
  			  } else {
  			    tmp <- NULL
  			  }

  			  # PCR testing and case record days
  			  intro_date <- end_date - T_sim_star + 1
  			  testing_days <- as.integer(testing_dates - intro_date + 1)
  			  sx_testing_days <- as.integer(sx_testing_dates - intro_date + 1)
  			  if (!is.null(CCMS_data)){
  			    case_days <- as.integer(case_dates - intro_date + 1)
  			  } else {
  			    case_days <- integer()
  			  }

  			  # Initialise infections
  			  TrueState <- rep(1,N_pop)
  			  # "Index" cases:
  			  e0 <- sample(res_present0,E0_star)
  			  # # 1st case with sx onset on 3/31
  			  # i_s_p0 <- sample(res_present0,1) # draw at random from residents who are present
  			  # TrueState[i_s_p0] <- 4 # assume initially in severe presymptomatic state
  			  # # 2nd case with sx onset on 4/2, assume initially in exposed state
  			  # e0 <- sample(setdiff(res_present0,i_s_p0),1) # draw at random from remaining residents who are present
  			  e0ind <- rep(F,N_pop)
  			  # e0ind[e0] <- T
  			  TrueState[e0] <- 2 # assume initially in latent state
  			  DayTrueState <- rep(0,N_pop)
  			  # DayTrueState[e0] <- 1 # assume 2nd index case has been in latent state for 1 day
  			  DayTrueState[e0] <- 0 # assume all index cases are at start of latent stage
  			  WaitingTime <- rep(NA,N_pop)
  			  # WaitingTime[i_s_p0] <- 2 # assume 1st index case has presymptomatic duration of 2 days (~mean of presymptomatic duration distribution) # rnbinom(length(i_s_p0),r_p,p_p)+1
  			  # WaitingTime[e0] <- 3 # assume 2nd index case has latent duration of 3 days (mean of latent duration distribution)
  			  WaitingTime[e0] <- rbinom(E0_star,r_E,p_E) + 1 # draw latent duration of 3 days (mean of latent duration distribution)
  			  DaysSinceInfctn <- rep(NA,N_pop) # days since infection
  			  # DaysSinceInfctn[i_s_p0] <- rbinom(length(i_s_p0),r_E,p_E)+1 # latent period of 1st index case prior to start of simulation
  			  DaysSinceInfctn[e0] <- DayTrueState[e0]
  			  DaysSinceInfctsnss <- rep(NA,N_pop) # days since start of infectiousness (i.e. start of presymptomatic infectious stage)
  			  # DaysSinceInfctsnss[i_s_p0] <- 0
  			  DaysPCRpos <- rep(NA,N_pop) # duration of PCR positivity (viraemia)
  			  # start_time <- Sys.time()
  				res <- COVID_homeless_model(N_res,N_staff,N_pop,T_sim_star,w,beta_star,epsilon,r_E,p_E,p_s,h,r_p,p_p,alpha,r_sx,p_sx,
  				                            p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,
  				                            hospitalisation,fit,fit_extrap,spec,testing_days,N_tested,sx_testing_days,N_sx_tested,
  				                            tmp,Number,Resident,Present,Risk,Age,e0ind,TrueState,DayTrueState,
  				                            WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos,shelter)
  				# end_time <- Sys.time()
  				# print(end_time - start_time)
  				D_S_star <- res$PCRpos_sx_testing
  				D_T_star <- res$PCRpos
  				# D_C_star <- sapply(res$cases[seq(1,T_sim_star) %in% case_days],function(x) rbinom(1,x,sx_reporting_rate))
  				D_C_star <- res$cases[seq(1,T_sim_star) %in% case_days]
  				# Calculate distances
  				distance[j,] <- calc_distance(D_S,D_S_star,D_T,D_T_star,D_C,D_C_star) # [ ] - UPDATE to save this in output
    			if((is.null(D_C) && all(distance[j,] <= c(epsilon_S[g],epsilon_T[g]))) | (!is.null(D_C) && all(distance[j,] <= c(epsilon_S[g],epsilon_T[g],epsilon_C[g])))){ # If distances are less than tolerance
    			  m <- m+1
    			}
  			}
    		if (m>0){
    		  # Update counter of accepted particles
    		  i <- i+1
    		  # Store results
    		  res.new[i,] <- c(beta_star,E0_star,T_sim_star)
    		  # res.new[i,] <- c(R0_star,E0_star,T_sim_star)
  				# res_list[[i]] <- res
  				infections_list[[i]] <- res$infections
  				cases_list[[i]] <- res$cases
  				PCRpos_sx_testing_mat[i,] <- res$PCRpos_sx_testing
  				PCRpos_mat[i,] <- res$PCRpos
  				cases_mat[i,] <- res$cases[seq(1,T_sim_star) %in% case_days]
  				state_list[[i]] <- res$state
  				sim_pop_list[[i]] <- res$sim_pop
    			# Calculate weights
    			w1 <- prod(sapply(1:ncpar, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])))
    			# print(w1)
    			if(g==1){
    			  # start_time <- Sys.time()
    				w2 <- 1
    				# end_time <- Sys.time()
    				# print(end_time - start_time)
    			} else {
    			  # start_time <- Sys.time()
    				# w2 <- sum(sapply(1:N, function(a) w.old[a]*dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
    			  # w2 <- sum(sapply(1:N, function(a) w.old[a]*mydtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
    				# w2 <- sum(sapply(1:N, function(a) w.old[a]*dmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma)))
    				# w2 <- sum(sapply(1:N, function(a) w.old[a]*dmvn(res.new[i,cpar], mu=res.old[a,cpar], sigma=sigma)))
    				w2 <- sum(sapply(idx, function(a) w.old[a]/marg_prob[(lm.low[dpar[1]]:lm.upp[dpar[1]])==res.new[i,dpar[1]],(lm.low[dpar[2]]:lm.upp[dpar[2]])==res.new[i,dpar[2]]]*dmvn(res.new[i,cpar], mu=res.old[a,cpar], sigma=sigma)))
    				# end_time <- Sys.time()
    				# print(end_time - start_time)
    			}
    			# print(w2)
    			w.new[i] <- (m/n)*w1/w2
    			print(paste0('Generation: ', g, ", particle: ", i))
    			# print(D_S)
    			# print(D_S_star)
    			# print(D_T)
    			# print(D_T_star)
    			# print(D_C)
    			# print(D_C_star)
    		}
    	}
	  }
	  # Update counter of proposed particles
	  k <- k+1
	  # end_tot_time <- Sys.time()
	  # print(end_tot_time - start_tot_time)
  }
  # Sigma <- vector("list",N)
  # for(p in 1:N){
  # 	Sigma[[p]] <- getSigmaNeighbours(N, M, res.new[p,], res.new, lm.low, lm.upp)
  # }
  res.old <- res.new
  w.old <- w.new/sum(w.new)
  # Calculate marginal probability for each combination of discrete parameters # [ ] - Think about how this could be generalised for higher dimensions (i.e. more than 2 discrete parameters)
  for (i1 in 1:(lm.upp[dpar[1]]-lm.low[dpar[1]]+1)){
    for (j1 in 1:(lm.upp[dpar[2]]-lm.low[dpar[2]]+1)){
      idx1 <- apply(res.new[,dpar],1,function(x) all(x==c((lm.low[dpar[1]]:lm.upp[dpar[1]])[i1],(lm.low[dpar[2]]:lm.upp[dpar[2]])[j1])))
      marg_prob[i1,j1] <- sum(w.old[idx1])
    }
  }
  # print(w.old)
  acc_rate[g] <- i/k
  cat("Acceptance rate for generation ", g, ": ", round(acc_rate[g],2), "\n", sep = "")

  # Calculate R0 estimates and add to output
  res.new1 <- cbind(res.new,R0=rep(NA,N))
  res.new1[,npar+1] <- sapply(1:N,function(i) calc_R0(res.new[i,1],w,Present,p_s,Risk,h,alpha,mu_p,mu_sx))
  # Save estimates from generation g
  write.csv(res.new1, file = paste("results_ABC_SMC_MNN_gen_",g,run_nm,".csv",sep=""), row.names=FALSE)
}

# Calculate ESS
ess <- ESS(w.old)

# save(infections_mat,cases_mat,PCRpos_sx_testing_mat,PCRpos_mat,state_arr,sim_pop_list,ess,file = "sim_output_10.RData")
save(infections_list,cases_list,PCRpos_sx_testing_mat,PCRpos_mat,cases_mat,state_list,sim_pop_list,ess,acc_rate,file = paste0("sim_output",run_nm,".RData"))

# Plot calibration results
plot_calibration(paste0("results_ABC_SMC_MNN_gen_10",run_nm,".csv"),paste0("sim_output",run_nm,".RData"),run_nm,lm.low,lm.upp,end_date,D_S,D_T,D_C,testing_dates,sx_testing_dates,tol)

# Process calibration results
process_calibration(paste0("results_ABC_SMC_MNN_gen_10",run_nm,".csv"),run_nm,end_date)

}

