# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Title: Comparison of public health strategies to reduce COVID-19 in homeless populations across the United States 
# Code author: Lloyd Chapman, Nathan C Lo, MD PhD (UCSF)
# Origin date: 5/19/20
# Last updated: 6/9/20
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

COVID_homeless_intervention_model<-function(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,
                                            alpha,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,
                                            max_days_PCR_pos,discrnorm,hospitalisation,fit,fit_extrap,spec,
                                            testing_days,interventions,max_PCR_tests_per_week,
                                            entry_PCR_test_compliance,routine_PCR_test_compliance,
                                            mask_compliance,mask_eff,sens_sx,spec_sx,Number,Alive,Resident,
                                            Present,Risk,Age,e0ind,TrueState,DayTrueState,WaitingTime,
                                            DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Population microsimulation 
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
  
  # Mask wearing
  if (5 %in% interventions) {
    beta <- (1-mask_compliance*mask_eff)*beta
    # print(beta)
  }
  
  # Removal of high-risk individuals and replacement with low-risk individuals
  if (6 %in% interventions){
    Risk[Risk %in% c(2,3,4)] <- 1
  }
  
  # Initialise variables 
  Hospitalised <- rep(F,N_pop)
  Dead <- rep(F,N_pop)
  ObsState <- rep(1,N_pop)
  DayObsState <- rep(0,N_pop)
  Tested <- rep(0,N_pop)
  DayTested <- rep(NA,N_pop)
  PCRtests <- rep(0,N_pop)
  PCRtestsWeek <- rep(0,N_pop) # counter of number of PCR tests individual has had in current week
  HxPCR <- rep(F,N_pop) # has tested positive on PCR - these individuals are removed and not allowed to return
  DayRemoved <- rep(NA,N_pop) # day individual removed following PCR positive test 
  HxAb <- rep(F,N_pop) # has tested positive for antibodies
  HxSx <- rep(F,N_pop) # has screened positive for symptoms
  NewInfection <- rep(0,N_pop)
  
  sim_pop0 <- data.frame(cbind(Number, Resident, Alive, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, DayTested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Decision trees
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Decision tree logic 
  
  # Time step: Daily 
  # Note: Each person has a "true state" (susceptible, exposed, presymptomatic infectious (mild/severe), symptomatic infectious (mild/severe), recovered) and an "observed state" (non-immune, infected or immune)
  
  # Create objects for storing simulation output
  infections <- rep(0,T_sim)
  cases <- rep(0,T_sim)
  infections_staff <- rep(0,T_sim)
  cases_staff <- rep(0,T_sim)
  PCRpos <- rep(0,T_sim)
  state <- matrix(NA, nrow = N_pop, ncol = T_sim)
  state[,1] <- TrueState
  presence <- matrix(NA, nrow = N_pop, ncol = T_sim)
  presence[,1] <- Present
  
  # Initialise vectors for storing individuals removed due to testing PCR positive
  sx_indvdls_removed <- integer()
  PCRpos_removed <- integer()
  # Initialise vector for individuals who are isolated due to screening positive for symptoms
  sx_indvdls_isolated <- integer()
    
  for (t in 2:T_sim) {
    
    # print(t)
    
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

    # Update states of individuals based on transmission on previous day
    list[TrueState,DayTrueState,DayObsState,WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,NewInfection,DaysPCRpos,infections[t],cases[t],infections_staff[t],cases_staff[t]] <- iterate(TrueState,Present,DayTrueState,DayObsState,DaysSinceInfctn,DaysSinceInfctsnss,beta,w,h,alpha,epsilon,N_pop,WaitingTime,r_E,p_E,e0ind,Risk,p_s,r_p,p_p,NewInfection,DaysPCRpos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,Resident,hospitalisation,Hospitalised,Dead,p_h,p_ICU,p_d)
    
    # Store true states of individuals on day t
    state[,t] <- TrueState
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 2: REMOVAL OF INDIVIDUALS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    # Remove all individuals who tested PCR positive on the previous day
    Present[c(sx_indvdls_removed,PCRpos_removed)] <- F
    # Return symptom positive individuals who tested PCR negative
    # print(intersect(sx_indvdls_isolated,PCRpos_removed))
    # print(setdiff(sx_indvdls_isolated,PCRpos_removed))
    Present[setdiff(sx_indvdls_isolated,PCRpos_removed)] <- T
    DayRemoved[c(sx_indvdls_removed,PCRpos_removed)] <- t

    if (2 %in% interventions){
      idx <- (Resident==1)
      list[Present[idx],,Numbers_add1] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==1,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add2] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==2,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add3] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==3,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add4] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==4,sample(-5:5,1),0)
      # print("Indvdls added")
      # print(c(Numbers_add1,Numbers_add2,Numbers_add3,Numbers_add4))
    }
    
    presence[,t] <- Present
    sx_indvdls_removed <- integer() # reset list of symptomatic individuals removed to empty
    PCRpos_removed <- integer() # reset list of PCR positive individuals removed to empty
    # print(c(sx_indvdls_removed,PCRpos_removed))
    # print(sum(Present[idx]))
    # print(summary(as.factor(Risk[Resident==1 & Present])))
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 3: INTERVENTIONS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Re-set counter for number of PCR tests for each individual in current week to 0 at end of each week
    if (t %% 7 == 0){
      PCRtestsWeek <- rep(0,N_pop) 
    }
    
    # PCR testing once upon entry
    if (2 %in% interventions){
      Numbers_add <- c(Numbers_add1,Numbers_add2,Numbers_add3,Numbers_add4)
      # print(Numbers_add)
      for (i in 1:7){
        list[PCRtests,Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_on_entry_update(i,TrueState,Number,Numbers_add,entry_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed) 
      }
      # print(PCRpos_removed)
    }
        
    # Active symptom screening
    if (1 %in% interventions){
      HxSx_prev <- HxSx 
      for (i in 1:7){
        list[HxSx,Tested,DayTested,PCRtests,PCRtestsWeek,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- sx_screening_update(i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,sens_sx[i],spec_sx[i],HxSx,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
      }
      sx_indvdls_isolated <- which(HxSx & !HxSx_prev) # individuals who screen positive for symptoms on this day 
      Present[sx_indvdls_isolated] <- F # immediately isolate individuals who screen symptom positive until test results are returned the next day
      # print(sx_indvdls_isolated)
      # print(PCRpos_removed)
    }
    
    # Passive symptom screening
    if (4 %in% interventions){
      idx_sx <- which(TrueState==6 & ObsState==1 & Present & DayTrueState==2) # severe symptomatic individuals present in shelter and not known to be infected self present for PCR testing on 3rd day of symptoms
      # print(idx_sx)
      PCRtests[idx_sx] <- PCRtests[idx_sx] + 1
      PCRtestsWeek[idx_sx] <- PCRtestsWeek[idx_sx] + 1
      list[Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,DayTested,idx_sx,6,spec[6],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
    }
    
    # Routine PCR testing
    if (3 %in% interventions){
      if (t %in% testing_days){
        for (i in 1:7){
          list[PCRtests,PCRtestsWeek,Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- routine_PCR_testing_update(i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,routine_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
        }        
      }
      # print(PCRpos_removed)
    }
  }
  
  sim_pop <- data.frame(cbind(Number, Resident, Alive, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, DayTested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection))
  return(res=list(infections=infections,cases=cases,infections_staff=infections_staff,cases_staff=cases_staff,PCRpos=PCRpos,sim_pop=sim_pop,state=state,presence=presence))
}
