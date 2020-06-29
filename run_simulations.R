run_simulations <- function(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx,N_res,N_staff,N_pop,T_sim,
                            epsilon,r_E,p_E,h,r_p,p_p,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,
                            min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,fit,
                            fit_extrap,spec,testing_days,interventions,max_PCR_tests_per_week,
                            entry_PCR_test_compliance,routine_PCR_test_compliance,
                            mask_compliance,mask_eff,sens_sx,spec_sx,Number,Alive,Resident,
                            Age,e0ind,TrueState,DayTrueState,WaitingTime,DaysSinceInfctn,
                            DaysSinceInfctsnss,DaysPCRpos,run_nm){
  # Calculate beta
  beta <- calc_beta(R0,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx)
  
  # Arrays to store number of infections and clinical cases over time from each simulation
  infections_res <- array(NA,c(nsims,T_sim,length(interventions)))
  cases_res <- array(NA,c(nsims,T_sim,length(interventions)))
  infections_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  cases_staff <- array(NA,c(nsims,T_sim,length(interventions)))
  
  # Matrices to store total number of infections and clinical cases from each simulation
  total_infections_res <- matrix(NA,nsims,length(interventions))
  total_cases_res <- matrix(NA,nsims,length(interventions))
  total_infections_staff <- matrix(NA,nsims,length(interventions))
  total_cases_staff <- matrix(NA,nsims,length(interventions))
  
  # Run simulations for each intervention strategy
  for (j in 1:length(interventions)){
    for (i in 1:nsims){
      res <- COVID_homeless_intervention_model(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,
                                               alpha,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,
                                               max_days_PCR_pos,discrnorm,hospitalisation,fit,fit_extrap,spec,
                                               testing_days,interventions[[j]],max_PCR_tests_per_week,
                                               entry_PCR_test_compliance,routine_PCR_test_compliance,
                                               mask_compliance,mask_eff,sens_sx,spec_sx,Number,Alive,Resident,
                                               Present,Risk,Age,e0ind,TrueState,DayTrueState,WaitingTime,
                                               DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos)
      infections_res[i,,j] <- res$infections
      cases_res[i,,j] <- res$cases
      total_infections_res[i,j] <- sum(res$infections)
      total_cases_res[i,j] <- sum(res$cases)
      infections_staff[i,,j] <- res$infections
      cases_staff[i,,j] <- res$cases
      total_infections_staff[i,j] <- sum(res$infections_staff)
      total_cases_staff[i,j] <- sum(res$cases_staff)
    }  
  }
  
  total_infections <- total_infections_res + total_infections_staff
  total_cases <- total_cases_res + total_cases_staff
  infections <- infections_res + infections_staff
  cases <- cases_res + cases_staff
  
  print(colMeans(total_infections))
  
  # Save total number of new infections per simulation
  # total_infections_df <- data.frame(rbind(total_infections,colMeans(total_infections)))
  total_infections_df <- data.frame(total_infections)
  names(total_infections_df) <- c("NoIntervention",sapply(1:(length(interventions)-1),function(x) paste0("Strategy",x)))
  # row.names(total_infections_df)[nsims+1] <- "mean"
  write.csv(total_infections_df,paste0("total_infections_interventions",run_nm,".csv"))
  
}
