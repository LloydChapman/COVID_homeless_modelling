# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Title: Comparison of public health strategies to reduce COVID-19 in homeless populations across the United States 
# Code author: Lloyd Chapman, Nathan C Lo, MD PhD (UCSF)
# Origin date: 5/19/20
# Last updated: 5/28/20
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Study background  
# 
# Objective: 
# 
# Study outcomes:
# 1) Number of infections
# 2) Number of hospitalizations
# 3) Number of deaths
# 4) Total healthcare spending 

# [ ] - CHECK use of sample() throughout code given how it works for integer vectors of length 1 in first argument! Change to resample() if necessary

COVID_homeless_intervention_model<-function(beta,CCMS_data,SF_case_data,interventions){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Population matrix- microsimulation 
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Columns
  #   1- Person number
  #   2- Alive (1), Dead (0) 
  #   3- Age
  #   4- Risk group: 1- very low; 2- low; 3- medium; 4- high
  #   5- True State: 0- dead; 1- susceptible; 2- exposed; 3- presx infectious mild; 4- presx infectious severe; 5- sx infectious mild ; 6- sx infectious severe; 7- recovered
  #   6- Day in TRUE state: (1,N)
  #   7- Waiting time in days in TRUE state
  #   8- Observed State: 1- non-immune; 2- infected; 3- immune (defined post confirmed infection)
  #   9- Day in OBSERVED state: (1,N)
  #   10- Hx positive PCR: 0- no; 1- yes (detected positive PCR)
  #   11- Hx positive Ab (serology): 0- no; 1- yes
  #   12- True new infection
  
  N_res <- 350
  N_staff <- 65
  N_pop <- N_res + N_staff
  start_date <- min(CCMS_data$Date)
  end_date <- max(CCMS_data$Date)
  T_sim <- 30 #as.integer(end_date - start_date + 1) #13 # testing 5 days after 1st cases confirmed on 4/5 + 2 days presymptomatic infectiousness (assumed) + 6 days of symptoms of 1st case  # 30 # Days; Total of 1 month 
  
  # Backfill CCMS data with assumed population of 255 prior to Mar 29 and same proportions in different risk categories as on Mar 29
  v <- round(N_res*CCMS_data[1,2:ncol(CCMS_data)]/CCMS_data[1,2],0)
  # CCMS_data0 <- cbind(Date=CCMS_data$Date[1] - seq(T_sim-nrow(CCMS_data)+1,1,by=-1),data.frame(v[rep(1,T_sim-nrow(CCMS_data)+1),]))
  # CCMS_data <- rbind(CCMS_data0,CCMS_data)
  
  # Replace NA's with 0's
  CCMS_data[is.na(CCMS_data)] <- 0
  
  # Transmission parameters
  # prob_infection <- 0.05 # daily probability of infection
  # # Wells-Riley parameters
  # q <-  # quantum generation rate
  # Q <-  # ventilation rate
  # p <-  # breathing rate
  # beta <- 0.004 #q*p/Q # transmission rate constant
  h <- 0.55 # relative infectiousness of asymptomatic individuals (individuals with mild symptoms) from Li Science 2020
  mean_daily_cases <- mean(SF_case_data$Case.Count[SF_case_data$Date>=start_date & SF_case_data$Date<=end_date]) # mean of confirmed cases for period of interest
  mean_daily_inc <- mean_daily_cases/881549 # population estimate from US Census Bureau [https://www.census.gov/quickfacts/sanfranciscocountycalifornia]
  epsilon <- 0 #0.001 # mean_daily_inc/0.14 # transmission rate outside shelter assuming 1/0.14=7.1x as many infections as confirmed cases from Li Science 2020 - [ ] consider making this time-dependent
  mu_E <- 3 #5 # mean of negative-binomially-distributed latent period - REDUCED TO 3 DAYS TO ACCOUNT FOR INCLUSION OF PRE-SYMPTOMATIC INFECTIOUS STATE WITH MEAN OF 2.3 DAYS 
  r_E <- 4 # shape parameter of negative-binomially-distributed latent period
  p_E <- r_E/(r_E+mu_E-1) # 'success' probability of negative-binomially-distributed latent period
  # p_s <- 0.2 # probability of developing severe symptoms - [ ] look up age-dependence and check this is representative for a homeless shelter
  
  # Vector of probabilities of developing severe symptoms by age group and co-morbidity status (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
  # p_s <- c((35*0.03 + 10*0.12)/45,(10*0.12+15*0.35)/25,(35*0.06 + 10*0.25)/45,(10*0.25+15*0.76)/25) # Values from Tuite CMAJ 2020 for 15+'s only as only adults in shelter  - [] UPDATE to finer-grained values from supp. mat. 
  p_s <- c((0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30,(0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30) # Values from Davies medRxiv 2020 - no co-morbidity dependence as no information on this in paper
  
  # Vector of probabilities of hospitalisation by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
  # p_h <- c((0.008+0.01+0.019+0.054)/4,(0.151+0.333+0.618)/3,(0.008+0.01+0.019+0.054)/4,(0.151+0.333+0.618)/3) # Values from Davies medRxiv 2020 - no co-morbidity dependence as no information on this
  p_h <- c((0.021+0.025+0.035+0.077)/4,(0.159+0.262+0.446)/3,(0.044+0.054+0.075+0.165)/4,(0.340+0.561+0.954)/3) # Values from Tuite CMAJ 2020
  
  p_ICU <- 0.261 # probability of ICU admission from Wang JAMA 2020
  # Vector of probabilities of death for hospitalised patients admitted to ICU by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
  p_d <- c((20*0.17+10*0.23+10*0.31)/40,(0.41+0.55+0.60)/3,(20*0.45+10*0.60+10*0.81)/40,(1+1+1)/3)
  
  mu_p <- 2.3 # mean of negative binomial presymptomatic period
  r_p <- 4 # shape parameter of negative-binomially-distributed presymptomatic period
  p_p <- r_p/(r_p+mu_p-1) # 'success' probability of negative-binomially-distributed presymptomatic period
  mu_sx <- 7 # mean of negative-binomially-distributed duration of symptoms - perhaps a bit long?
  r_sx <- 4 # shape parameter of negative-binomially-distributed duration of symptoms
  p_sx <- r_sx/(r_sx+mu_sx-1) # 'success' probability of negative-binomially-distributed duration of symptoms
  mean_days_PCR_pos <- 20 # mean duration of PCR positivity following symptom onset based on viral shedding dynamics papers
  # r_PCR_pos <- 10 # shape parameter of negative-binomially-distributed duration of PCR positivity - [ ] just guessed ATM, should be taken from viral shedding dynamics papers
  # p_PCR_pos <- r_PCR_pos/(r_PCR_pos+mean_days_PCR_pos-1)
  days_PCR_pos_presx <- 2 # duration of PCR positivity prior to symptom onset
  
  quaran_period <- 14+1
  hospitalisation <- T # flag for whether hospitalisation was recorded - not available yet for MSC outbreak
  
  alpha <- 0.44/(1-0.44)*mu_sx/mu_p # relative infectiousness of presymptomatic infection (He 2020 Nat Med)
  
  N_staff_cohorts <- 3
  
  # Construct discrete distribution for duration of PCR positivity from reversed log-normal distribution
  # x <- 1:40
  # plot(x,dtrunc(x-1,"nbinom",a=0,b=20,size=r_PCR_pos,prob=p_PCR_pos),xlab="Duration of PCR positivity (days)",ylab="Prob",)
  min_days_PCR_pos <- 5 # He Nat Med 2020, Xiao Jrnl Clin Vir 2020
  max_days_PCR_pos <- 37 # Zhou Lancet 2020 #25
  # discrlnorm <- discretize(plnorm(x,log(max_days_PCR_pos-mean_days_PCR_pos+1),log(2)), from = 0, to = max_days_PCR_pos, step = 1, method = "unbiased" , lev = levlnorm(x,log(max_days_PCR_pos-mean_days_PCR_pos+1),log(2)))
  discrnorm <- discretize(pnorm(x,mean_days_PCR_pos,5), from = min_days_PCR_pos-0.5, to = max_days_PCR_pos+0.5, step = 1, method = "rounding")
  discrnorm <- discrnorm/sum(discrnorm) # normalise
  # plot(min_days_PCR_pos:max_days_PCR_pos,discrnorm)
  # # plot(0:max_days_PCR_pos,discrlnorm)
  # revdiscrlnorm <- rev(discrlnorm[1:(length(discrlnorm)-1)])
  # revdiscrlnorm <- revdiscrlnorm/sum(revdiscrlnorm) # normalise
  # # pdf("viraemia_duration_distn.pdf",width = 5,height = 4)
  # # plot(1:max_days_PCR_pos,revdiscrlnorm,xlab = "Duration of viraemia following symptom onset (days)", ylab = "Prob",pch=19) #, main = "Distribution of duration of viraemia")
  # # dev.off()
  
  # Intervention parameters
  max_PCR_tests_per_week <- 10 # maximum number of PCR tests per week 
  
  # PCR testing once upon entry
  entry_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing on entry
  
  # Routine PCR testing
  routine_PCR_test_compliance <- 1 #0.8 # 80% compliance with routine PCR testing
  
  # Mask wearing
  mask_compliance <- 1 #0.8 # 80% compliance with universal masking
  mask_eff <- 1 #0.3 # 30% reduction in transmission from universal masking
  if (5 %in% interventions) {
    beta <- (1-mask_compliance*mask_eff)*beta
    # print(beta)
  }
  
  # Testing frequency
  # freq_test <- T_sim # test once at the end of the simulation #30 # 5 # test every 5 days
  # testing_days <- seq(1,T_sim)[seq(1,T_sim) %% freq_test==0]
  testing_dates <- as.Date(c("4/8/2020","4/9/2020"),format="%m/%d/%Y") # testing dates
  # testing_dates <- as.Date(c("4/9/2020"),format="%m/%d/%Y") # testing dates - assume all testing happened on 4/9 for now
  # testing_days <- as.integer(testing_dates - start_date + 1) # days on which testing occurred
  testing_days <- seq(1,T_sim)[(seq(1,T_sim) %% 7) %in% seq(0,6)] # testing twice per week on 1st and 4th day
  N_tested <- rep((144-10)/2,2) # number tested during mass testing = total tested from press release - number tested during first 4 days = 144 - 10, # number tested on each testing day #c(72,72) # assume all testing happened during mass testing and was equally split between 4/8 and 4/9
  PCRpos <- rep(0,T_sim) #rep(0,length(testing_days))
  sx_testing_dates <- as.Date(c("4/4/2020","4/5/2020","4/6/2020","4/7/2020"),format="%m/%d/%Y")
  sx_testing_days <- as.integer(sx_testing_dates - start_date + 1)
  PCRpos_sx_testing <- rep(0,length(sx_testing_days))
  
  # # Sensitivity as a function of days since symptom onset
  # sens <- function(x){ifelse(x>=0 & x<=21, 1-x/21, 0)} # linearly decreasing function of time since onset for both mild and severe sx - [ ] Update to exponentially decreasing function based on He 2020 Nat Med paper
  # Sensitivity as a function of days since infection
  FNR <- read.csv("../Data/digitised_sens_graph.csv",header = F, stringsAsFactors = F)
  x <- round(FNR[(round(mu_E)+1):nrow(FNR),1],0)-round(mu_E)
  y <- 1 - FNR[(round(mu_E)+1):nrow(FNR),2]
  fit <- lm(y ~ bs(x,knots = c(3,4,6,7,10,15,17)-round(mu_E)))
  fit_extrap <- approxExtrap(x,y,(x[length(x)]+1):max_days_PCR_pos)
  
  # Symptom screening sensitivity and specificity
  sens_sx <- c(NA,NA,NA,NA,NA,0.9,NA) # sensitivities for states 1 to 7
  spec_sx <- c(0.75,0.75,0.75,0.75,0.75,NA,0.75) # specificities for states 1 to 7
  
  # PCR test specificity 
  spec <- c(1,1,NA,NA,NA,NA,NA) # specificities for states 1 to 7
  
  # [ ] Update Initial conditions 
  # [ ] - UPDATE to include 2nd case 
  # Initialise variables
  Number <- 1:N_pop
  Alive <- rep(1,N_pop)
  Resident <- rep(0,N_pop)
  Resident[1:N_res] <- 1
  res_present0 <- sample(1:N_res,CCMS_data$Total_Census[1])
  res_absent0 <- setdiff(1:N_res,res_present0)
  # StaffCohort <- rep(NA,N_pop)
  # StaffCohort[(N_res+1):N_pop] <- pmin(1 + floor((1:N_staff)/N_staff*N_staff_cohorts),N_staff_cohorts)
  # staff_present0 <- sample((N_res+1):N_pop,floor(N_staff/N_staff_cohorts))
  # staff_absent0 <- setdiff((N_res+1):N_pop,staff_present0)
  staff_present0 <- (N_res+1):N_pop
  # Present <- rep(T,N_pop)
  Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
  Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
  # Hi_Risk_60_only <- sample(1:N_res,CCMS_data$Hi_Risk_60_only[1])
  Hi_Risk_60_only_present0 <- sample(res_present0,CCMS_data$Hi_Risk_60_only[1])
  Hi_Risk_60_only_absent0 <- sample(res_absent0,v$Hi_Risk_60_only - CCMS_data$Hi_Risk_60_only[1])
  Risk[c(Hi_Risk_60_only_present0,Hi_Risk_60_only_absent0)] <- 2
  # Hi_Risk_Dx_only <- sample(setdiff(1:N_res,Hi_Risk_60_only),CCMS_data$Hi_Risk_Dx_only[1])
  Hi_Risk_Dx_only_present0 <- sample(setdiff(res_present0,Hi_Risk_60_only_present0),CCMS_data$Hi_Risk_Dx_only[1])
  Hi_Risk_Dx_only_absent0 <- sample(setdiff(res_absent0,Hi_Risk_60_only_absent0),v$Hi_Risk_Dx_only - CCMS_data$Hi_Risk_Dx_only)
  Risk[c(Hi_Risk_Dx_only_present0,Hi_Risk_Dx_only_absent0)] <- 3
  # Hi_Risk_Both_Age_Dx <- sample(setdiff(1:N_res,union(Hi_Risk_60_only,Hi_Risk_Dx_only)),CCMS_data$Hi_Risk_Both_Age_Dx[1])
  Hi_Risk_Both_Age_Dx_present0 <- sample(setdiff(res_present0,c(Hi_Risk_60_only_present0,Hi_Risk_Dx_only_present0)),CCMS_data$Hi_Risk_Both_Age_Dx[1])
  Hi_Risk_Both_Age_Dx_absent0 <- sample(setdiff(res_absent0,c(Hi_Risk_60_only_absent0,Hi_Risk_Dx_only_absent0)),v$Hi_Risk_Both_Age_Dx-CCMS_data$Hi_Risk_Both_Age_Dx[1])
  Risk[c(Hi_Risk_Both_Age_Dx_present0,Hi_Risk_Both_Age_Dx_absent0)] <- 4
  # print(summary(as.factor(Risk[res_present0])))
  # print(summary(as.factor(Risk[res_absent0])))
  
  # Removal of high-risk individuals and replacement with low-risk individuals
  if (6 %in% interventions){
    Risk[Risk %in% c(2,3,4)] <- 1
  }
  
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
  # DaysPCRpos[i_s_p0] <- sample(1:max_days_PCR_pos,length(i_s_p0),prob=revdiscrlnorm,replace=T)
  Hospitalised <- rep(F,N_pop)
  Dead <- rep(F,N_pop)
  ObsState <- rep(1,N_pop)
  DayObsState <- rep(0,N_pop)
  Tested <- rep(0,N_pop)
  PCRtests <- rep(0,N_pop)
  PCRtestsWeek <- rep(0,N_pop) # counter of number of PCR tests individual has had in current week
  HxPCR <- rep(F,N_pop) # has tested positive on PCR - these individuals are removed and not allowed to return
  DayRemoved <- rep(NA,N_pop) # day individual removed following PCR positive test 
  HxAb <- rep(F,N_pop) # has tested positive for antibodies
  HxSx <- rep(F,N_pop) # has screened positive for symptoms
  NewInfection <- rep(0,N_pop)
  
  sim_pop0 <- data.frame(cbind(Number, Resident, Alive, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection))
  
  # Presence weight
  w <- c(rep(1,N_res),rep(1/2,N_staff))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Decision trees
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Decision tree logic 
  
  # Time step: Daily 
  # Note: Each person has a "true state" (susceptible, exposed, presymptomatic infectious (mild/severe), symptomatic infectious (mild/severe), recovered) and an "observed state" (non-immune, infected or immune)
  
  # Counters
  new_infections <- rep(0,T_sim)
  state <- matrix(NA, nrow = N_pop, ncol = T_sim)
  state[,1] <- TrueState
  presence <- matrix(NA, nrow = N_pop, ncol = T_sim)
  presence[,1] <- Present
  
  # Set marker for whether there were any individuals with severe symptoms on days individuals with clinical symptoms were tested to false
  sx_indvdls_removed <- vector()
  PCRpos_removed <- vector()
    
  for (t in 2:T_sim) {
    
    print(t)
    
    # Advance time by one day - [ ] Think carefully about whether this should be at start or end of loop
    DayTrueState <- DayTrueState + 1
    DayObsState <- DayObsState + 1
    
    # Create index variables
    s <- (TrueState == 1)
    e <- (TrueState == 2)
    i_m_p <- (TrueState == 3)
    i_s_p <- (TrueState == 4)
    i_m_sx <- (TrueState == 5)
    i_s_sx <- (TrueState == 6)
    r <- (TrueState == 7)
    
    # Count numbers in different states inside shelter
    S <- sum(s & Present)
    E <- sum(e & Present)
    I_m_p <- sum(i_m_p & Present)
    I_s_p <- sum(i_s_p & Present)
    I_m_sx <- sum(i_m_sx & Present)
    I_s_sx <- sum(i_s_sx & Present)
    R <- sum(r & Present)
    
    print(paste0("sum(i_m_p)=",sum(i_m_p)))
    print(paste0("I_m_p=",I_m_p))
    print(paste0("sum(i_s_p)=",sum(i_s_p)))
    print(paste0("I_s_p=",I_s_p))
    print(paste0("sum(i_m_sx)=",sum(i_m_sx)))
    print(paste0("I_m_sx=",I_m_sx))
    print(paste0("sum(i_s_sx)=",sum(i_s_sx)))
    print(paste0("I_s_sx=",I_s_sx))
    
    # Increase days since infection by one day for infected individuals
    DaysSinceInfctn[e|i_m_p|i_s_p|i_m_sx|i_s_sx|r] <- DaysSinceInfctn[e|i_m_p|i_s_p|i_m_sx|i_s_sx|r] + 1
    # Increase days since start of infectiousness by one day for infected individuals
    DaysSinceInfctsnss[i_m_p|i_s_p|i_m_sx|i_s_sx|r] <- DaysSinceInfctsnss[i_m_p|i_s_p|i_m_sx|i_s_sx|r] + 1
    
    # Re-set counter for number of PCR tests for each individual in current week to 0 at end of each week
    if (t %% 7 == 0){
      PCRtestsWeek <- rep(0,N_pop) 
    }
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # ALGORITHM 1: TRANSMISSION MODEL # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 
    # Focus on TRUE state
    # 
    # Step 1: S -> E (infection)
    # 
    # Step 2: E -> I_m_p, I_s_p (progression to presymptomatic infectious stage)
    # 
    # Step 3: I_m_p -> I_m_sx, I_s_p -> I_s_sx (symptom onset)
    # 
    # Step 4: I_m_sx, I_s_sx -> R (recovery)
    # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    # Columns
    #   1- Person number
    #   2- Alive (1), Dead (0) 
    #   3- Age
    #   4- Risk group: 1- very low; 2- low; 3- medium; 4- high - [ ] UPDATE defs
    #   5- True State: 0- dead; 1- susceptible; 2- exposed; 3- presx infectious mild; 4- presx infectious severe; 5- sx infectious mild ; 6- sx infectious severe; 7- recovered
    #   6- Day in TRUE state: (1,N)
    #   7- Observed State: 1- non-immune; 2- infected; 3- immune (defined post confirmed infection)
    #   8- Day in OBSERVED state: (1,N)
    #   9- Hx positive PCR: 0- no; 1- yes (detected positive PCR)
    #   10- Hx positive Ab (serology): 0- no; 1- yes
    #   11- True new infection
    
    # Step 1: S -> E (infection)
    lambda <- beta*Present*w*(h*(alpha*I_m_p+I_m_sx)+alpha*I_s_p+I_s_sx) + epsilon # Only individuals present in the shelter can be infected by others in the shelter, those not present can be infected by infectious individuals in the general population
    print(paste0("sum(lambda)=",sum(lambda)))
    prob_infection <- 1-exp(-lambda)
    s2e <- (s & runif(N_pop)<prob_infection)
    
    TrueState[s2e] <- 2 # Update true state based on transmission/natural history
    DayTrueState[s2e] <- 0 # Re-set day in true state
    WaitingTime[s2e] <- rnbinom(sum(s2e),r_E,p_E) + 1 # Draw latent periods
    # DaysPCRpos[s2e] <- WaitingTime[s2e] + sample(1:max_days_PCR_pos,sum(s2e),prob=revdiscrlnorm,replace=T)
    DaysSinceInfctn[s2e] <- 0
    
    # Step 2: E -> I_m_p, I_s_p (progression to presymptomatic infectious stage)
    e2i_p <- (e & (DayTrueState==WaitingTime) & !e0ind)
    
    for (i in 1:4){
      TrueState[e2i_p & Risk==i] <- 3 + (runif(sum(e2i_p & Risk==i)) < p_s[i]) # Update true state based on transmission/natural history
    }
    DayTrueState[e2i_p] <- 0 # Re-set day in true state 
    WaitingTime[e2i_p] <- rnbinom(sum(e2i_p),r_p,p_p) + 1 # Draw infectious periods
    DaysSinceInfctsnss[e2i_p] <- 0
    
    # Treat progression of 2nd index case separately
    e02i_p <- (e & (DayTrueState==WaitingTime) & e0ind)
    TrueState[e02i_p] <- 4 # assume 2nd index case had severe symptoms
    DayTrueState[e02i_p] <- 0 # Re-set day in true state
    WaitingTime[e02i_p] <- 2 # assume 2nd index case has presymptomatic duration of 2 days
    DaysSinceInfctsnss[e02i_p] <- 0
    
    # Step 3: I_m_p -> I_m_sx, I_s_p -> I_s_sx (symptom onset)
    i_p2i_sx <- ((i_m_p | i_s_p) & (DayTrueState==WaitingTime))
    
    TrueState[i_p2i_sx] <- TrueState[i_p2i_sx] + 2 # Update true state based on transmission/natural history: just add 2 for both transitions 3->5, 4->6
    DayTrueState[i_p2i_sx] <- 0 # Re-set day in true state
    NewInfection[i_p2i_sx] <- 1 # Count new infections 
    WaitingTime[i_p2i_sx] <- rnbinom(sum(i_p2i_sx),r_sx,p_sx) + 1
    DaysPCRpos[i_p2i_sx] <- sample(min_days_PCR_pos:max_days_PCR_pos,sum(i_p2i_sx),prob=discrnorm,replace=T)
    
    # Step 4: I_m_sx, I_s_sx -> R (recovery)
    i_s2r <- ((i_m_sx | i_s_sx) & (DayTrueState==WaitingTime))
    
    TrueState[i_s2r] <- 7 # Update true state based on transmission/natural history
    DayTrueState[i_s2r] <- 0 # Re-set day in true state
    
    if (hospitalisation){
      for (j in 1:4){
        i_s_sx_Risk_j <- (Present & i_s_sx & i_s2r & Risk==j)
        i_hosp <- (runif(sum(i_s_sx_Risk_j)) < p_h[j])
        if (any(i_hosp)){
          Hospitalised[i_s_sx_Risk_j] <- i_hosp
          Present[i_s_sx_Risk_j] <- !i_hosp
          i_d <- ((runif(sum(i_hosp)) < p_ICU) & (runif(sum(i_hosp)) < p_d[j]))
          Dead[i_s_sx_Risk_j][i_hosp] <- i_d
        }
      }    
    }
    
    # Count number of new infections on day t
    new_infections[t] <- sum(i_p2i_sx)
    
    # Store true states of individuals on day t
    state[,t] <- TrueState
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 2: REMOVAL OF INDIVIDUALS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    Present[c(sx_indvdls_removed,PCRpos_removed)] <- F # All symptomatic individuals who tested PCR positive on the previous day are removed
    DayRemoved[c(sx_indvdls_removed,PCRpos_removed)] <- t
    # # Present[1:N_res] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==1,sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==1) + sample(-5:5,1),0)
    # # Present[1:N_res] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==2,sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==2) + sample(-5:5,1),0)
    # # Present[1:N_res] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==3,sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==3) + sample(-5:5,1),0)
    # # Present[1:N_res] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==4,sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==4) + sample(-5:5,1),0)
    # list[Present[1:N_res],,Numbers_add1] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==1,sample(-5:5,1),0)
    # list[Present[1:N_res],,Numbers_add2] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==2,sample(-5:5,1),0)
    # list[Present[1:N_res],,Numbers_add3] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==3,sample(-5:5,1),0)
    # list[Present[1:N_res],,Numbers_add4] <- update_present(Present[1:N_res],HxPCR[1:N_res],Risk[1:N_res]==4,sample(-5:5,1),0)
    # # print("Indvdls added")
    # # print(c(Numbers_add1,Numbers_add2,Numbers_add3,Numbers_add4))
    
    presence[,t] <- Present
    sx_indvdls_removed <- vector() # reset list of symptomatic individuals removed to empty
    PCRpos_removed <- vector() # reset list of PCR positive individuals removed to empty
    # print(c(sx_indvdls_removed,PCRpos_removed))
    # print(sum(Present[1:N_res]))
    # print(summary(as.factor(Risk[Resident==1 & Present])))
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 3: INTERVENTIONS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # PCR testing once upon entry
    if (2 %in% interventions){
      Numbers_add <- c(Numbers_add1,Numbers_add2,Numbers_add3,Numbers_add4)
      # print(Numbers_add)
      for (i in 1:7){
        list[PCRtests,Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_on_entry_update(i,TrueState,Number,Numbers_add,entry_PCR_test_compliance,PCRtests,Tested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed) 
      }
      # print(PCRpos_removed)
    }
        
    # Active symptom screening
    if (1 %in% interventions){
      for (i in 1:7){
        list[HxSx,Tested,PCRtests,PCRtestsWeek,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- sx_screening_update(i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,sens_sx[i],spec_sx[i],HxSx,PCRtests,Tested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
      }
      # print(PCRpos_removed)
    }
    
    # Passive symptom screening
    if (4 %in% interventions){
      idx_sx <- which(TrueState==6 & ObsState==1 & Present & DayTrueState==2) # severe symptomatic individuals present in shelter and not known to be infected self present for PCR testing on 3rd day of symptoms
      # print(idx_sx)
      PCRtests[idx_sx] <- PCRtests[idx_sx] + 1
      PCRtestsWeek[idx_sx] <- PCRtestsWeek[idx_sx] + 1
      list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx_sx,6,spec[6],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
    }
    
    # Routine PCR testing
    if (3 %in% interventions){
      for (i in 1:7){
        list[PCRtests,PCRtestsWeek,Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- routine_PCR_testing_update(t,testing_days,i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,routine_PCR_test_compliance,PCRtests,Tested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,PCRpos_removed)
      }
      # print(PCRpos_removed)
    }
  }
  
  sim_pop <- data.frame(cbind(Number, Resident, Alive, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection))
  return(res=list(new_infections=new_infections,PCRpos_sx_testing=PCRpos_sx_testing,PCRpos=PCRpos,sim_pop=sim_pop,state=state,presence=presence))
}
