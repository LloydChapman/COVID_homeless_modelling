run_simulations <- function(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,nsims,N_res,N_staff,N_pop,
                            T_sim,epsilon,r_E,p_E,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                            min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,sens,
                            spec,testing_days,interventions,max_PCR_tests_per_week,
                            min_days_btw_tests,entry_PCR_test_compliance,routine_PCR_test_compliance,
                            sx_pos_PCR_test_compliance,mask_compliance,mask_eff_susc,mask_eff_inf,sens_sx,spec_sx,Number,
                            Resident,Age,res_present0,E0,run_nm,dir=""){
  # Calculate beta
  beta <- calc_beta(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
  
  # Set initial conditions
  TrueState <- rep(1,N_pop)
  # "Index" case(s):
  e0 <- sample(res_present0,E0) # draw E0 residents to be initially infected individuals
  e0ind <- rep(F,N_pop) # don't fix times or states for progression of initial infections
  TrueState[e0] <- 2 # assume initially in latent state
  DayTrueState <- rep(0,N_pop)
  DayTrueState[e0] <- 0 # assume all index cases are at start of latent stage
  WaitingTime <- rep(NA,N_pop)
  WaitingTime[e0] <- rbinom(E0,r_E,p_E) + 1 # draw latent duration
  DaysSinceInfctn <- rep(NA,N_pop) # days since infection
  DaysSinceInfctn[e0] <- DayTrueState[e0]
  DaysSinceInfctsnss <- rep(NA,N_pop) # days since start of infectiousness (i.e. start of early infectious stage)
  DaysPCRpos <- rep(NA,N_pop) # duration of detectable viral load from start of infectiousness
  
  # Arrays to store number of infections and clinical cases over time from each simulation
  infections_res <- array(NA,c(nsims,T_sim,length(interventions)))
  cases_res <- array(NA,c(nsims,T_sim,length(interventions)))
  infections_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  cases_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  
  # Arrays to store number of infections and clinical cases over time from outside shelter from each simulation
  bckgrnd_infections_res <- array(NA,c(nsims,T_sim,length(interventions)))
  bckgrnd_cases_res <- array(NA,c(nsims,T_sim,length(interventions)))
  bckgrnd_infections_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  bckgrnd_cases_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  
  # Matrices to store total number of infections, clinical cases, hospitalisations and deaths from each simulation
  total_infections_res <- matrix(NA,nsims,length(interventions))
  total_cases_res <- matrix(NA,nsims,length(interventions))
  total_infections_staff <- matrix(NA,nsims,length(interventions))
  total_cases_staff <- matrix(NA,nsims,length(interventions))
  total_hospitalisations <- matrix(NA,nsims,length(interventions))
  total_deaths <- matrix(NA,nsims,length(interventions))
  
  # Array to store number of PCR tests per individual from each simulation
  PCRtests <- array(NA,c(N_pop,nsims,length(interventions)))
  
  # Run simulations for each intervention strategy
  for (j in 1:length(interventions)){
    for (i in 1:nsims){
      set.seed(i)
      res <- COVID_homeless_intervention_model(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,
                                               alpha,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,
                                               max_days_PCR_pos,discrnorm,hospitalisation,sens,spec,
                                               testing_days,interventions[[j]],max_PCR_tests_per_week,
                                               min_days_btw_tests,entry_PCR_test_compliance,
                                               routine_PCR_test_compliance,sx_pos_PCR_test_compliance,
                                               mask_compliance,mask_eff_susc,mask_eff_inf,sens_sx,spec_sx,Number,Resident,Present,
                                               Risk,Age,e0ind,TrueState,DayTrueState,WaitingTime,DaysSinceInfctn,
                                               DaysSinceInfctsnss,DaysPCRpos)
      infections_res[i,,j] <- res$infections_res
      cases_res[i,,j] <- res$cases_res
      bckgrnd_infections_res[i,,j] <- res$bckgrnd_infections_res
      bckgrnd_cases_res[i,,j] <- res$bckgrnd_cases_res
      total_infections_res[i,j] <- sum(res$infections_res)
      total_cases_res[i,j] <- sum(res$cases_res)
      infections_staff[i,,j] <- res$infections_staff
      cases_staff[i,,j] <- res$cases_staff
      bckgrnd_infections_staff[i,,j] <- res$bckgrnd_infections_staff
      bckgrnd_cases_staff[i,,j] <- res$bckgrnd_cases_staff
      total_infections_staff[i,j] <- sum(res$infections_staff)
      total_cases_staff[i,j] <- sum(res$cases_staff)
      total_hospitalisations[i,j] <- sum(res$sim_pop$Hospitalised)
      total_deaths[i,j] <- sum(res$sim_pop$Dead)
      PCRtests[,i,j] <- res$sim_pop$PCRtests
    }  
  }

  infections <- infections_res + infections_staff
  cases <- cases_res + cases_staff
  bckgrnd_infections <- bckgrnd_infections_res + bckgrnd_infections_staff
  bckgrnd_cases <- bckgrnd_cases_res + bckgrnd_cases_staff
  total_infections <- total_infections_res + total_infections_staff
  total_cases <- total_cases_res + total_cases_staff
  
  print(colMeans(total_infections))
  
  # Save total number of new infections and new cases per simulation
  total_infections_df <- data.frame(total_infections)
  names(total_infections_df) <- c("NoIntervention",sapply(1:(length(interventions)-1),function(x) paste0("Strategy",x)))
  if (dir!=""){
    dir.create(dir,recursive = T)
  }
  write.csv(total_infections_df,paste0(dir,"total_infections_interventions",run_nm,".csv"))
  total_cases_df <- data.frame(total_cases)
  names(total_cases_df) <- c("NoIntervention",sapply(1:(length(interventions)-1),function(x) paste0("Strategy",x)))
  write.csv(total_cases_df,paste0(dir,"total_cases_interventions",run_nm,".csv"))
  
  save(infections_res,cases_res,infections_staff,cases_staff,bckgrnd_infections_res,bckgrnd_cases_res,bckgrnd_infections_staff,bckgrnd_cases_staff,infections,cases,bckgrnd_infections,bckgrnd_cases,total_infections_res,total_cases_res,total_infections_staff,total_cases_staff,total_infections,total_cases,total_hospitalisations,total_deaths,PCRtests,file=paste0(dir,"intvntn_sim_output",run_nm,".RData"))
}
