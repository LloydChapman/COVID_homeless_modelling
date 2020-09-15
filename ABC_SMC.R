ABC_SMC <- function(shelter,N,run_nm){

# Set data and initial conditions for shelter
set_data_and_init_condns(shelter)

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
parnames <- c("R0","E0","T")
# Indices of continuous and discrete parameters
cpar <- 1
dpar <- 2:3
ncpar <- length(cpar)

# Fixed variance for R0 proposals
sigma <- 1 #2

# Empty matrices to store results (population x number of model parameters)
res.old <- matrix(nrow=N,ncol=npar)
res.new <- matrix(nrow=N,ncol=npar)
colnames(res.old) <- parnames
colnames(res.new) <- parnames
# Empty vector to store beta estimates
res.beta <- numeric(N)

# Empty vectors to store weights
w.old <- numeric(N)
w.new <- numeric(N)

# Empty objects to store outputs for each simulation
infections_list <- vector("list",N) # daily numbers of new infections
cases_list <- vector("list",N) # daily numbers of new clinical cases (by symptom onset)
PCRpos_sx_testing_mat <- matrix(nrow=N,ncol=length(D_S)) # daily number PCR positive during symptomatic testing 4/4-4/7
PCRpos_mat <- matrix(nrow=N,ncol=length(D_T)) # number PCR positive in mass testing
cases_mat <- matrix(nrow=N,ncol=length(D_C)) # number of clinical cases on days with reported clinical cases 
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
   	if(g==1){
    	# Sample from prior distribution for E0 and T_sim
   	  E0_star <- sample(lm.low[2]:lm.upp[2],1)
   	  T_sim_star <- sample(lm.low[3]:lm.upp[3],1)
   	  R0_star <- runif(1, min=lm.low[1], max=lm.upp[1])
			beta_star <- calc_beta(R0_star,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
		} else {
		  # Sample from prior distribution for E0 and T_sim
	    E0_star <- sample(lm.low[2]:lm.upp[2],1)
	    T_sim_star <- sample(lm.low[3]:lm.upp[3],1)
	    #  Find particles in previous generation with (E0,T_sim) combination drawn
	    idx <- which(apply(res.old[,dpar],1,function(x) all(x==c(E0_star,T_sim_star))))
	    if (length(idx)!=0){ # if there are any remaining particles with these values of the discrete parameters
	      # Sample particle from previous generation according to particle weights
	      p <- resample(idx,1,prob=w.old[idx])
	      # Perturb particle
	      par <- rK(res.old[p,cpar],sigma,lm.low,lm.upp)
	      R0_star <- par[1]
	      beta_star <- calc_beta(R0_star,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
	    }
		}
	  # Test if there are any remaining particles with the proposed values of the discrete parameters
	  if (g==1 || (g>1 & length(idx)!=0)){
    	#  Test if prior is non-zero
	    if(prior.non.zero(c(R0_star,E0_star,T_sim_star),lm.low,lm.upp)) {
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
  			  e0 <- sample(res_present0,E0_star) # draw E0 residents to be initially infected individuals
  			  e0ind <- rep(F,N_pop) # don't fix times or states for progression of initial infections
  			  TrueState[e0] <- 2 # assume initially in latent state
  			  DayTrueState <- rep(0,N_pop)
  			  DayTrueState[e0] <- 0 # assume all index cases are at start of latent stage
  			  WaitingTime <- rep(NA,N_pop)
  			  WaitingTime[e0] <- rbinom(E0_star,r_E,p_E) + 1 # draw latent duration
  			  DaysSinceInfctn <- rep(NA,N_pop) # days since infection
  			  DaysSinceInfctn[e0] <- DayTrueState[e0]
  			  DaysSinceInfctsnss <- rep(NA,N_pop) # days since start of infectiousness (i.e. start of presymptomatic infectious stage)
  			  DaysPCRpos <- rep(NA,N_pop) # duration of detectable viral load from start of infectiousness
  				res <- COVID_homeless_model(N_res,N_staff,N_pop,T_sim_star,w,beta_star,epsilon,r_E,p_E,p_s,h,r_p,p_p,alpha,r_sx,p_sx,
  				                            p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,
  				                            hospitalisation,sens,spec,testing_days,N_tested,sx_testing_days,N_sx_tested,
  				                            tmp,Number,Resident,Present,Risk,Age,e0ind,TrueState,DayTrueState,
  				                            WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos,shelter)
  				D_S_star <- res$PCRpos_sx_testing
  				D_T_star <- res$PCRpos
  				case_days_idx <- (seq(1,T_sim_star) %in% case_days)
  				D_C_star <- res$cases_res[case_days_idx] + res$cases_staff[case_days_idx]
  				# Calculate distances
  				distance[j,] <- calc_distance(D_S,D_S_star,D_T,D_T_star,D_C,D_C_star)
    			if((is.null(D_C) && all(distance[j,] <= c(epsilon_S[g],epsilon_T[g]))) | (!is.null(D_C) && all(distance[j,] <= c(epsilon_S[g],epsilon_T[g],epsilon_C[g])))){ # If distances are less than tolerance
    			  m <- m+1
    			}
  			}
    		if (m>0){
    		  # Update counter of accepted particles
    		  i <- i+1
    		  # Store results
    		  res.new[i,] <- c(R0_star,E0_star,T_sim_star)
    		  res.beta[i] <- beta_star
  				infections_list[[i]] <- res$infections_res + res$infections_staff
  				cases_list[[i]] <- res$cases_res + res$cases_staff
  				PCRpos_sx_testing_mat[i,] <- D_S_star
  				PCRpos_mat[i,] <- D_T_star
  				cases_mat[i,] <- D_C_star
  				state_list[[i]] <- res$state
  				sim_pop_list[[i]] <- res$sim_pop
    			# Calculate weights
    			w1 <- prod(sapply(1:ncpar, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])))
    			# print(w1)
    			if(g==1){
    				w2 <- 1
    			} else {
    				w2 <- sum(sapply(idx, function(a) w.old[a]/marg_prob[(lm.low[dpar[1]]:lm.upp[dpar[1]])==res.new[i,dpar[1]],(lm.low[dpar[2]]:lm.upp[dpar[2]])==res.new[i,dpar[2]]]*dmvn(res.new[i,cpar], mu=res.old[a,cpar], sigma=sigma)))
    			}
    			w.new[i] <- (m/n)*w1/w2
    			print(paste0('Generation: ', g, ", particle: ", i))
    		}
    	}
	  }
	  # Update counter of proposed particles
	  k <- k+1
  }
  res.old <- res.new
  w.old <- w.new/sum(w.new)
  # Calculate marginal probability for each combination of discrete parameters
  for (i1 in 1:(lm.upp[dpar[1]]-lm.low[dpar[1]]+1)){
    for (j1 in 1:(lm.upp[dpar[2]]-lm.low[dpar[2]]+1)){
      idx1 <- apply(res.new[,dpar],1,function(x) all(x==c((lm.low[dpar[1]]:lm.upp[dpar[1]])[i1],(lm.low[dpar[2]]:lm.upp[dpar[2]])[j1])))
      marg_prob[i1,j1] <- sum(w.old[idx1])
    }
  }
  acc_rate[g] <- i/k
  cat("Acceptance rate for generation ", g, ": ", round(acc_rate[g],2), "\n", sep = "")

  res.new1 <- cbind(res.new,beta=res.beta)
  # Save estimates from generation g
  write.csv(res.new1, file = paste("results_ABC_SMC_MNN_gen_",g,run_nm,".csv",sep=""), row.names=FALSE)
}

# Calculate ESS
ess <- ESS(w.old)
print(paste0(shelter,": ESS = ",ess))

# save(infections_mat,cases_mat,PCRpos_sx_testing_mat,PCRpos_mat,state_arr,sim_pop_list,ess,file = "sim_output_10.RData")
save(infections_list,cases_list,PCRpos_sx_testing_mat,PCRpos_mat,cases_mat,state_list,sim_pop_list,ess,acc_rate,file = paste0("sim_output",run_nm,".RData"))

# Plot calibration results
ttl <- ifelse(shelter=="SF","San Francisco",gsub("_"," ",shelter))
plot_calibration(paste0("results_ABC_SMC_MNN_gen_10",run_nm,".csv"),paste0("sim_output",run_nm,".RData"),run_nm,lm.low,lm.upp,end_date,D_S,D_T,D_C,testing_dates,sx_testing_dates,tol,ttl)

# Process calibration results
# Overwrite size of simulated population for SF shelter with number of unique individuals that were actually present in shelter between start and end dates
if (shelter=="SF"){
  N_pop <- 255+65
}
process_calibration(paste0("results_ABC_SMC_MNN_gen_10",run_nm,".csv"),paste0("sim_output",run_nm,".RData"),run_nm,end_date,date_first_case,N_pop)

}
