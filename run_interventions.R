library(actuar)
library(splines)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(gsubfn)

setwd("~/Dropbox/Homeless/Code")
source("COVID_homeless_interventions.R")
source("COVID_homeless_functions.R")
source("run_simulations.R")
source("process_interventions.R")
source("plot_interventions.R")

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

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
SF_case_data <- aggregate(Case.Count ~ Date,SF_data,sum)

# Set number of residents and staff in shelter and duration of simulation
N_res <- 250 #237 #350
N_staff <- 50 #65
N_pop <- N_res + N_staff
T_sim <- 30

# Set weights for presence of residents and staff in shelter
w <- rep(1,N_pop) # c(rep(1,N_res),rep(1/2,N_staff))

# Set natural history parameters
source("set_nat_hist_pars.R")

# Set background transmission rate
epsilon <- 0 #0.001

# Flag for whether to count hospitalisations and deaths
hospitalisation <- T

# Set PCR test parameters
source("set_PCR_test_pars.R")

# PCR testing frequency
testing_days <- seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% seq(1,4)] # seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% seq(0,6)] # testing twice per week on 1st and 4th day

# Set intervention parameters
max_PCR_tests_per_week <- 2 #10 # maximum number of PCR tests per week
min_days_btw_tests <- 3 # minimum number of days between PCR tests

# PCR testing once upon entry
entry_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing on entry

# Routine PCR testing
routine_PCR_test_compliance <- 0.8 # 1 # 80% compliance with routine PCR testing

# Mask wearing
mask_compliance <- 0.8 # 1 # 80% compliance with universal masking
mask_eff <- 0.3 # 1 # 30% reduction in transmission from universal masking

# Symptom screening sensitivity and specificity
sens_sx <- c(NA,NA,NA,NA,NA,0.4,NA) # sensitivities for states 1 to 7
spec_sx <- c(0.75,0.75,0.75,0.75,0.75,NA,0.75) # specificities for states 1 to 7
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

## Run interventions with median R0's from calibration to different shelter outbreaks
run_nms <- c("_SF_R0_6","_Seattle_A_R0_6","_Boston_R0_6")
# Names of files with posterior samples
fnms <- sapply(c("_16","_Seattle_A_6","_Boston_5"),function(x) paste0("results_ABC_SMC_MNN_gen_10",x,".csv"))
for (i in 1:length(fnms)){
  pars <- read.csv(fnms[i],stringsAsFactors = F)
  R0 <- median(pars[,1])
  run_simulations(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
                  T_sim,epsilon,r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                  min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
                  fit_extrap,spec,testing_days,interventions,max_PCR_tests_per_week,
                  min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
                  mask_compliance,mask_eff,sens_sx,spec_sx,Number,Resident,
                  Age,res_present0,E0,run_nms[i])
  process_interventions(paste0("intvntn_sim_output",run_nms[i],".RData"),run_nms[i])
  plot_interventions(paste0("intvntn_sim_output",run_nms[i],".RData"),run_nms[i])
}

# # Run interventions for hypothetical shelter with lower R0 and lower background transmission based on SF data
# R0 <- 1.5
# mean_daily_cases <- mean(SF_case_data$Case.Count[SF_case_data$Date>=as.Date("3/28/2020",format="%m/%d/%Y") & SF_case_data$Date<=as.Date("4/10/2020",format="%m/%d/%Y")]) # mean of confirmed cases for period of interest
# mean_daily_inc <- mean_daily_cases/881549 # population estimate from US Census Bureau [https://www.census.gov/quickfacts/sanfranciscocountycalifornia]
# epsilon <- mean_daily_inc/0.14 # transmission rate outside shelter assuming 1/0.14=7.1x as many infections as confirmed cases from Li Science 2020
# run_nm <- "_lower_R0_3"
# run_simulations(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
#                 T_sim,epsilon,r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
#                 min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
#                 fit_extrap,spec,testing_days,interventions,max_PCR_tests_per_week,
#                 min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
#                 mask_compliance,mask_eff,sens_sx,spec_sx,Number,Resident,
#                 Age,e0ind,TrueState,DayTrueState,WaitingTime,DaysSinceInfctn,
#                 DaysSinceInfctsnss,DaysPCRpos,run_nm)
# process_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm)
# plot_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm)
# 
# # Quality checks
# # Run interventions for hypothetical shelter with lower R0 and lower background transmission based on SF data
# R0 <- 2.5
# mean_daily_cases <- mean(SF_case_data$Case.Count[SF_case_data$Date>=as.Date("3/28/2020",format="%m/%d/%Y") & SF_case_data$Date<=as.Date("4/10/2020",format="%m/%d/%Y")]) # mean of confirmed cases for period of interest
# mean_daily_inc <- mean_daily_cases/881549 # population estimate from US Census Bureau [https://www.census.gov/quickfacts/sanfranciscocountycalifornia]
# epsilon <- 6*mean_daily_inc/0.14 # transmission rate outside shelter assuming 1/0.14=7.1x as many infections as confirmed cases from Li Science 2020
# run_nm <- "_high_risk_removal_check"
# run_simulations(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
#                 T_sim,epsilon,r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
#                 min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
#                 fit_extrap,spec,testing_days,list(NULL,6),max_PCR_tests_per_week,
#                 min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
#                 mask_compliance,mask_eff,sens_sx,spec_sx,Number,Resident,
#                 Age,e0ind,TrueState,DayTrueState,WaitingTime,DaysSinceInfctn,
#                 DaysSinceInfctsnss,DaysPCRpos,run_nm)
# process_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm)
# plot_interventions(paste0("intvntn_sim_output",run_nm,".RData"),run_nm)
# 
# # Daily PCR testing with no max number of weekly tests # [ ] - changes hard coded in interventions code, need to UPDATE
# infections_PCR <- numeric(nsims)
# for (j in 1:nsims){
#   res3 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,c(3,4))
#   infections_PCR[j] <- sum(res3$infections)
# }
# 
# barplot(res3$infections)
# View(res3$sim_pop)
# View(res3$state)
# View(res3$presence)
# tE3 <- getEventTime1(res3$state,2)
# tI_m_p3 <- getEventTime1(res3$state,3)
# tI_s_p3 <- getEventTime1(res3$state,4)
# View(cbind(tE3,tI_m_p3,tI_s_p3,res3$sim_pop$DayRemoved))
# mean(infections_PCR)
# 
# # write.table(infections_PCR,"new_onsets_daily_PCR_no_max_per_week.csv",sep = ",",col.names = F, row.names = F)
# 
# ## beta=0.01, epsilon=0, no migration in/out only hospitalisation of some clinical cases 
# # No interventions
# set.seed(1)
# res4 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,NULL)
# View(res4$sim_pop)
# View(res4$state)
# View(res4$presence)
# plot(colSums(res4$presence))
# barplot(res4$infections)
# sum(res4$infections) # 302
# # Passive sx screening and daily PCR testing with 100% compliance and sensitivity=1
# set.seed(1)
# res5 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,c(3,4))
# View(res5$sim_pop)
# View(res5$state)
# View(res5$presence)
# plot(colSums(res5$presence))
# barplot(res5$infections)
# sum(res5$infections) # 295
# 
# ## beta=0.005
# # No interventions
# set.seed(1)
# res6 <- COVID_homeless_intervention_model(0.005,CCMS_data,SF_case_data,NULL)
# View(res6$sim_pop)
# View(res6$state)
# View(res6$presence)
# plot(colSums(res6$presence))
# barplot(res6$infections)
# sum(res6$infections) # 302
# # Passive sx screening and daily PCR testing with 100% compliance and sensitivity=1
# set.seed(1)
# res7 <- COVID_homeless_intervention_model(0.005,CCMS_data,SF_case_data,c(3,4))
# View(res7$sim_pop)
# View(res7$state)
# View(res7$presence)
# plot(colSums(res7$presence))
# barplot(res7$infections)
# sum(res7$infections) # 260
# 
# ## 100% effective masking
# set.seed(1)
# res8 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,5)
# View(res8$sim_pop)
# View(res8$state)
# View(res8$presence)
# plot(colSums(res8$presence))
# barplot(res8$infections)
# sum(res8$infections) # 2 - only the initial cases
# 
# 
