rm(list=ls())

library(actuar)
library(splines)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(gsubfn)
library(doParallel)
library(abind)

setwd("~/Dropbox/Homeless/Code")
source("COVID_homeless_interventions.R")
source("COVID_homeless_functions.R")
source("run_simulations.R")
source("process_interventions.R")
source("plot_interventions.R")

# Register doParallel backend with 10 workers
registerDoParallel(10)

# Number of simulations per intervention strategy
nsims <- 1000

## Set intervention strategies
# 1 - Active symptom screening
# 2 - PCR testing once upon entry
# 3 - Routine PCR testing
# 4 - Passive symptom screening
# 5 - Universal masking
# 6 - Removal of high-risk individuals
# 7 - Routine PCR testing of staff only
interventions <- list(NULL,1,c(1,3),c(1,5),c(1,6),c(1,7),c(1,3,5,6))

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

# Start and end date for calculating background incidence
start_date <- as.Date("7/4/2020",format = "%m/%d/%Y")
end_date <- as.Date("7/17/2020",format = "%m/%d/%Y")

underreporting <- 4
homeless_RR <- 2

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition_2020_07_24.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
names(SF_data)[names(SF_data)=="Specimen.Collection.Date"] <- "Date"
names(SF_data)[names(SF_data)=="Case.Count"] <- "Cases"
SF_data$Date <- as.Date(SF_data$Date) 
SF_case_data <- aggregate(Cases ~ Date,SF_data,sum)
epsilon_SF <- calc_epsilon(SF_case_data,start_date,end_date,881549,underreporting,homeless_RR)

# Load Boston case data obtained from https://dashboard.cityofboston.gov/t/Guest_Access_Enabled/views/COVID-19/Dashboard1?:showAppBanner=false&:display_count=n&:showVizHome=n&:origin=viz_share_link&:isGuestRedirectFromVizportal=y&:embed=y
Boston_data <- read.delim("../Data/Boston_Cases_over_time_data.csv",sep = "\t",stringsAsFactors = F, fileEncoding = "UTF-16")
Boston_case_data <- data.frame(Date = as.Date(Boston_data$Day.of.Timestamp,format = "%B %d, %Y"))
Boston_case_data$Cases <- c(1,diff(Boston_data$Cases..Total.Positive,lag = 1))
epsilon_Boston <- calc_epsilon(Boston_case_data,start_date,end_date,692600,underreporting,homeless_RR)

# Load Seattle case data
Seattle_case_data <- read.csv("../Data/Seattle_cases.csv",stringsAsFactors = F)
Seattle_case_data$Date <- as.Date(Seattle_case_data$Date,format = "%m/%d/%y")
epsilon_Seattle <- calc_epsilon(Seattle_case_data,start_date,end_date,753675,underreporting,homeless_RR)

# Set number of residents and staff in shelter and duration of simulation
N_res <- 250 #237 #350
N_staff <- 50 #65
N_pop <- N_res + N_staff
T_sim <- 30

# Set weights for presence of residents and staff in shelter
w <- rep(1,N_pop) # c(rep(1,N_res),rep(1/2,N_staff))

# Set natural history parameters
source("set_nat_hist_pars.R")

# Set background transmission rates
tmp <- c(epsilon_SF,epsilon_Seattle,epsilon_Boston)
epsilons <- c(0,min(tmp),mean(tmp),max(tmp))

# Flag for whether to count hospitalisations and deaths
hospitalisation <- T

# Set PCR test parameters
source("set_PCR_test_pars.R")

# PCR testing frequency
testing_freq <- 2 # testing events per week
testing_days <- floor(seq(1,T_sim,by = 7/testing_freq)) #seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% seq(1,4)] #seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% c(1,4)] # seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% seq(0,6)] # testing twice per week on 1st and 4th day

# Set intervention parameters
max_PCR_tests_per_week <- 2 #10 # maximum number of PCR tests per week
min_days_btw_tests <- 3 # minimum number of days between PCR tests

# PCR testing once upon entry
entry_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing on entry

# Routine PCR testing
routine_PCR_test_compliance <- 0.8 # 1 # 80% compliance with routine PCR testing
sx_pos_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing among those who screen symptom positive

# Mask wearing
mask_compliance <- 0.8 # 1 # 80% compliance with universal masking
mask_eff <- 0.3 # 1 # 30% reduction in transmission from universal masking

# Symptom screening sensitivity and specificity
sens_sx <- c(NA,NA,NA,NA,NA,0.4,NA) # sensitivities for states 1 to 7
spec_sx <- c(0.9,0.9,0.9,0.9,0.9,NA,0.9) # specificities for states 1 to 7
# sens_sx <- c(NA,NA,NA,NA,NA,1,NA) # sensitivities for states 1 to 7
# spec_sx <- c(1,1,1,1,1,NA,1) # specificities for states 1 to 7

# Initialise variables
Number <- 1:N_pop
Resident <- rep(0,N_pop)
Resident[1:N_res] <- 1
res_present0 <- 1:N_res
res_absent0 <- setdiff(1:N_res,res_present0)
staff_present0 <- (N_res+1):N_pop
Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
v <- round(N_res*CCMS_data[1,2:ncol(CCMS_data)]/CCMS_data[1,2],0)
Hi_Risk_60_only_present0 <- sample(res_present0,CCMS_data$Hi_Risk_60_only[1])
Risk[Hi_Risk_60_only_present0] <- 2
Hi_Risk_Dx_only_present0 <- sample(setdiff(res_present0,Hi_Risk_60_only_present0),CCMS_data$Hi_Risk_Dx_only[1])
Risk[Hi_Risk_Dx_only_present0] <- 3
Hi_Risk_Both_Age_Dx_present0 <- sample(setdiff(res_present0,c(Hi_Risk_60_only_present0,Hi_Risk_Dx_only_present0)),CCMS_data$Hi_Risk_Both_Age_Dx[1])
Risk[Hi_Risk_Both_Age_Dx_present0] <- 4
Age <- rep(NA,N_pop)
Age[Risk %in% c(1,3)] <- sample(x=seq(20,59), size=sum(Risk %in% c(1,3)), replace=TRUE)
Age[Risk %in% c(2,4)] <- sample(x=seq(60,80), size=sum(Risk %in% c(2,4)), replace=TRUE) # [ ] - CHECK oldest age (have assumed 69 for now)

# Set number of initial latently infected cases
E0 <- 1

# Names of files with posterior samples
fnms <- sapply(c("_Seattle_A_12","_Boston_11","_SF_22"),function(x) paste0("results_ABC_SMC_MNN_gen_10",x,".csv"))
# Make vector of different estimated R0's
R0s <- numeric(length(fnms))
for (j in 1:length(fnms)){
  pars <- read.csv(fnms[j],stringsAsFactors = F)
  R0s[j] <- median(pars[,1])  
}
# Add a lower R0 value for scenario analyses
R0s <- c(1.5,R0s)
R0lbls <- c("(low-risk)","(Seattle)","(Boston)","(SF)")
# R0s <- R0s[c(2,3)]
ttls <- c("No interventions","Symptom screening","Routine PCR testing","Universal mask wearing","Relocation of high-risk individuals","Routine PCR testing of staff only","Combination strategy")

## Run interventions with median R0's from calibration to different shelter outbreaks
run_nms <- c("_lower_R0_5","_Seattle_A_R0_10","_Boston_R0_10","_SF_R0_10")
for (i in 1:length(epsilons)){
  for (j in 1:length(R0s)){
    run_nm <- paste0(run_nms[j],"_epsilon",i-1)
    run_simulations(R0s[j],w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
                    T_sim,epsilons[i],r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                    min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
                    fit_extrap,spec,testing_days,interventions,max_PCR_tests_per_week,
                    min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
                    sx_pos_PCR_test_compliance,mask_compliance,mask_eff,sens_sx,spec_sx,Number,
                    Resident,Age,res_present0,E0,run_nm)
    process_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm,T_sim)
    plot_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm,ttls)
  }
  # Combine results for each epsilon value
  combine_results("perc_reduction_infections",length(interventions)-1,run_nms,i)
  combine_results("perc_reduction_cases",length(interventions)-1,run_nms,i)
  combine_results("prob_outbreak_averted",length(interventions)-1,run_nms,i)
}

plot_interventions_by_strategy(nsims,T_sim,run_nms,epsilons,interventions,R0s,R0lbls,ttls)

# Run interventions for a range of background incidences for each of the R0 estimates to find relationship
# between background incidence and probability of averting an outbreak
epsilons1 <- seq(0,max(tmp),length.out = 10) #seq(0,max(tmp),length.out = 30)
# epsilons1 <- seq(0,min(tmp),length.out = 8)
# prob_outbreak_averted <- array(NA,dim = c(length(interventions)-1,length(epsilons1),length(run_nms)))

acomb <- function(...){abind(...,along = 3)}

prob_outbreak_averted <- foreach(i=1:length(R0s),.combine = 'acomb',.multicombine=T) %:%
 foreach(j=1:length(epsilons1),.combine = 'cbind') %dopar% {
   run_nm <- paste0(run_nms[i],"_SA_1_epsilon",j-1)
   # run_nm <- paste0(run_nms[i],"_SA_test1_epsilon",j-1)
   run_simulations(R0s[i],w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
                   T_sim,epsilons1[j],r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                   min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
                   fit_extrap,spec,testing_days,interventions,max_PCR_tests_per_week,
                   min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
                   sx_pos_PCR_test_compliance,mask_compliance,mask_eff,sens_sx,spec_sx,Number,
                   Resident,Age,res_present0,E0,run_nm)
   load(paste0("intvntn_sim_output",run_nm,".RData"))
   # prob_outbreak_averted[,j,i] <- calc_prob_outbreak_averted(infections,bckgrnd_infections)
   calc_prob_outbreak_averted(infections,bckgrnd_infections)
 }

# save(interventions,R0s,epsilons1,prob_outbreak_averted,file = "prob_outbreak_averted_SA_epsilon_1.RData")
save(interventions,R0s,epsilons1,prob_outbreak_averted,file = "prob_outbreak_averted_SA_1_epsilon.RData")
# save(interventions,R0s,epsilons1,prob_outbreak_averted,file = "prob_outbreak_averted_SA_test1_epsilon_1.RData")

# plot_epsilon_SA("prob_outbreak_averted_SA_epsilon_1.RData",run_nms)
plot_epsilon_SA("prob_outbreak_averted_SA_1_epsilon.RData",run_nms,homelessRR,R0lbls)

# Run sensitivity analysis for probability of averting an outbreak under different PCR testing frequencies for each R0 value
# testing_freqs <- seq(0.5,2,length.out = 13)
days_btw_tests <- c(1,3) #c(1,3,7,10,14,21,30)

prob_outbreak_averted1 <-
  foreach(i=1:length(R0s),.combine = 'cbind') %:%
    # foreach(j=1:length(testing_freqs),.combine = 'c') %dopar% {
    foreach(j=1:length(days_btw_tests),.combine = 'c') %dopar% {
# # prob_outbreak_averted1 <- matrix(nrow = length(testing_freqs),ncol = length(R0s))
# # for (i in 1:length(R0s)){
# #   for(j in 1:length(testing_freqs)){
# prob_outbreak_averted1 <- matrix(nrow = length(days_btw_tests),ncol = length(R0s))
# for (i in 1:length(R0s)){
#    for(j in 1:length(days_btw_tests)){
      run_nm <- paste0(run_nms[i],"_SA_6_PCR_testing_freq",j)
      # testing_days <- floor(seq(1,T_sim,by = 7/testing_freqs[j]))
      testing_days <- floor(seq(days_btw_tests[j],T_sim,by = days_btw_tests[j]))
      run_simulations(R0s[i],w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
                      T_sim,epsilons[2],r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                      min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
                      fit_extrap,spec,testing_days,interventions[c(1,3)],max_PCR_tests_per_week,
                      min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
                      sx_pos_PCR_test_compliance,mask_compliance,mask_eff,sens_sx,spec_sx,Number,
                      Resident,Age,res_present0,E0,run_nm)
      load(paste0("intvntn_sim_output",run_nm,".RData"))
      calc_prob_outbreak_averted(infections,bckgrnd_infections)
    #   prob_outbreak_averted1[j,i] <- calc_prob_outbreak_averted(infections,bckgrnd_infections)
    # }
}

# save(R0s,testing_freqs,prob_outbreak_averted1,file = "prob_outbreak_averted_SA_5_PCR_testing_freq.RData")
# plot_PCR_testing_freq_SA("prob_outbreak_averted_SA_5_PCR_testing_freq.RData","_5")

save(R0s,days_btw_tests,prob_outbreak_averted1,file = "prob_outbreak_averted_SA_6_PCR_testing_freq.RData")
plot_PCR_testing_freq_SA("prob_outbreak_averted_SA_6_PCR_testing_freq.RData","_6",R0lbls)
