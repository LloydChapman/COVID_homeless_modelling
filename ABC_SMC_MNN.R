# install.packages(c("tmvtnorm","actuar"))
library(tmvtnorm)
library(actuar)
library(splines)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox/Homeless/Code") 
source("ABC_SMC_MNN_functions.R")
source("COVID_homeless_calibration2.R")

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
N_pop <- 320

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
N <- 200

# Number of neighbours for covariance matrix calculations
M <- 50

# Epsilon tolerance values for number PCR positive in early testing of symptomatics
epsilon_S <- round(seq(sqrt(sum((2*D_S)^2)),sqrt(sum((0.2*D_S)^2)),length.out = 10))

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

for(g in 1:G){  

	#Initiate counter
	i <- 1	
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
    				if((distance[j,1] <= epsilon_S[g]) & (distance[j,2] <= epsilon_T[g])){ # If distances is less than tolerance
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
      			# Update counter
      			i <- i+1
      			}
    		} 
    	}
    	Sigma <- list(NA, N)
    for(p in 1:N){
    	Sigma[[p]] <- getSigmaNeighbours(N, M, res.new[p,], res.new, lm.low, lm.upp) 
    }
    	res.old <- res.new
	w.old <- w.new/sum(w.new)

 	write.csv(res.new, file = paste("results_ABC_SMC_MNN_gen_",g,"_7.csv",sep=""), row.names=FALSE)
}

# Calculate ESS
ess <- ESS(w.old)

save(new_infections_mat,PCRpos_sx_testing_mat,PCRpos_mat,state_arr,sim_pop_list,ess,file = "sim_output_7.RData")

# Plot output
pdf("beta_posterior_7.pdf",width = 5, height = 4)
hist(res.new,breaks = 20,freq = F,xlab = "beta",main = "Posterior distribution of beta")
dev.off()

# pdf("num_PCRpos_vs_beta_7.pdf",width = 5,height = 5)
# plot(res.new,PCRpos_mat,xlab = "beta",ylab = "No. PCR positive")
# dev.off()

pdf("onsets_vs_time_7.pdf",width = 6,height = 5)
new_infections_melt <- melt(new_infections_mat)
ggplot(new_infections_melt,aes(x=Var2,y=value,col=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. of onsets") + theme(legend.position = "none")
dev.off()

# Reconstruct the number of PCR positive individuals over time from individual-level simulation data frame
getEventTime <- function(state,i){apply(state,c(1,3),function(x) which(x %in% i)[1])}

tE <- getEventTime(state_arr,2)
tI_m_p <- getEventTime(state_arr,3)
tI_s_p <- getEventTime(state_arr,4)
tI_m_sx <- getEventTime(state_arr,5)
tI_s_sx <- getEventTime(state_arr,6)
tR <- getEventTime(state_arr,7)

StartPresx <- pmin(tI_m_p,tI_s_p,na.rm = T)
Onset <- pmin(tI_m_sx,tI_s_sx,na.rm = T)
EndPCRpos <- matrix(nrow = nrow(Onset), ncol = ncol(Onset))
for (i in 1:N){
  EndPCRpos[,i] <- Onset[,i] + sim_pop_list[[i]]$DaysPCRpos
}

detectable_viral_load_mat <- matrix(nrow = N, ncol = T_sim)
for (i in 1:T_sim){
  detectable_viral_load_mat[,i] <- colSums(StartPresx<=i & EndPCRpos>i,na.rm = T)
}

detectable_viral_load_melt <- melt(detectable_viral_load_mat)
pdf("num_detectabe_viral_load_vs_time_7.pdf",width = 5, height = 4)
ggplot() + geom_line(aes(x=Var2,y=value,col=as.factor(Var1)),detectable_viral_load_melt) + geom_line(aes(x=day,y=PCRpos),data.frame(day=c(12,12),PCRpos=range(PCRpos_mat))) + geom_point(aes(x=12,y=D_T),col="black",size=2) + xlab("Day") + ylab("No. with detectable viral load") + theme(legend.position = "none")
dev.off()

testing_dates <- c(as.Date(c("4/4/2020","4/5/2020","4/6/2020","4/7/2020"),format="%m/%d/%Y"),as.Date(c("4/8/2020","4/9/2020"),format="%m/%d/%Y"))
D <- c(D_S,D_T)
PCRpos_mat1 <- data.frame(cbind(PCRpos_sx_testing_mat,PCRpos_mat))
names(PCRpos_mat1) <- testing_dates
PCRpos_mat1$sim <- row.names(PCRpos_mat1)
PCRpos_melt1 <- melt(PCRpos_mat1,id.vars="sim")
PCRpos_melt1$variable <- as.Date(PCRpos_melt1$variable)
pdf("calibration_7.pdf",width = 6, height = 5)
ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim),col=as.factor(sim)),PCRpos_melt1) + geom_line(aes(x=date,y=PCRpos),data.frame(date=testing_dates,PCRpos=D),col="black") + geom_line(aes(x=date,y=CI),data.frame(date=rep(testing_dates[5],2),CI=c(D_T[1]-tol1,D_T[1]+tol1))) + geom_line(aes(x=date,y=CI),data.frame(date=rep(testing_dates[6],2),CI=c(D_T[2]-tol2,D_T[2]+tol2))) + xlab("Date") + ylab("No. PCR positive") + theme(legend.position = "none")
dev.off()

num_PCRpos_staff <- rep(0,length(sim_pop_list))
for (i in 1:length(sim_pop_list)){
  num_PCRpos_staff[i] <- sum(sim_pop_list[[i]]$HxPCR[sim_pop_list[[i]]$Resident==0],na.rm = T) 
}
pdf("num_PCRpos_staff_7.pdf",width = 5, height = 4)
hist(num_PCRpos_staff,breaks=seq(min(num_PCRpos_staff)-0.5,max(num_PCRpos_staff)+0.5),xlab="No. PCR positive staff",main="Distribution of no. of PCR positive staff")
dev.off()
