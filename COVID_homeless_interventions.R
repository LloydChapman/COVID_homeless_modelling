COVID_homeless_intervention_model<-function(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,
                                            alpha,r_sx,p_sx,p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,
                                            max_days_PCR_pos,discrnorm,hospitalisation,sens,spec,
                                            testing_days,interventions,max_PCR_tests_per_week,min_days_btw_tests,
                                            entry_PCR_test_compliance,routine_PCR_test_compliance,
                                            sx_pos_PCR_test_compliance,mask_compliance,mask_eff,sens_sx,spec_sx,
                                            Number,Resident,Present,Risk,Age,e0ind,TrueState,DayTrueState,
                                            WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Population microsimulation 
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Note: Each person has a "true state" (TrueState): susceptible, exposed, early infectious (subclinical/clinical), late infectious (subclinical/clinical), recovered; and an "observed state" (ObsState): never infected or ever infected
  
  # Mask wearing
  if (5 %in% interventions) {
    beta <- (1-mask_compliance*mask_eff)*beta
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
  Tested <- rep(0,N_pop) # ever PCR tested
  DayTested <- rep(NA,N_pop) # day individual last PCR tested
  PCRtests <- rep(0,N_pop) # total number of PCR tests individual has had
  PCRtestsWeek <- rep(0,N_pop) # counter of number of PCR tests individual has had in current week
  HxPCR <- rep(F,N_pop) # has tested positive on PCR - these individuals are removed and not allowed to return
  DayRemoved <- rep(NA,N_pop) # day individual removed following PCR positive test 
  HxAb <- rep(F,N_pop) # has tested positive for antibodies
  HxSx <- rep(F,N_pop) # has screened positive for symptoms
  NewInfection <- rep(0,N_pop) # true new infection
  Background <- rep(0,N_pop) # was infected outside the shelter
  
  # Make data frame of initial values of variables for all individuals
  # Columns
  #   1- Person number
  #   2- Resident (1), Staff (0)
  #   3- Present (T), Absent (F)
  #   4- Age
  #   5- Risk group: (1) under-60 & no co-morbidities, (2) 60+ & no co-morbidities, (3) under-60 & co-morbidities, (4) 60+ & co-morbidities
  #   6- True State: (1) susceptible; (2) exposed; (3) subclinical early infectious; (4) clinical early infectious; (5) subclinical late infectious; (6) clinical late infectious; (7) recovered
  #   7- Days elapsed in TRUE state
  #   8- Waiting time in days in TRUE state
  #   9- Number of days since infection
  #   10- Number of days since becoming infectious (entering early infectious stage)
  #   11- Days from start of infectiousness for which individual has detectable viral load
  #   12- Hospitalised (T) or Not-hospitalised (F) if clinical symptoms
  #   13- Dead (T), Alive (F)
  #   14- Observed State: Uninfected (1); Infected (2)
  #   15- Days elapsed in OBSERVED state
  #   16- Tested (1), Untested (0) 
  #   17- Day last PCR tested
  #   18- Total number of PCR tests
  #   19- Number of PCR tests in current week
  #   20- History of positive PCR test: Yes (T), No (F)
  #   21- Day removed from shelter following positive PCR test
  #   22- History of positive Ab test (serology): Yes (T), No (F)
  #   23- History of positive symptom screen: Yes (T), No (F)
  #   23- True new infection
  #   24- Infected outside shelter: Yes (1), No (0)
  sim_pop0 <- data.frame(cbind(Number, Resident, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, DayTested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection, Background))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Decision trees
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Decision tree logic 
  
  # Time step: Daily
  
  # Create objects for storing simulation output
  infections_res <- rep(0,T_sim)
  cases_res <- rep(0,T_sim)
  infections_staff <- rep(0,T_sim)
  cases_staff <- rep(0,T_sim)
  bckgrnd_infections_res <- rep(0,T_sim)
  bckgrnd_cases_res <- rep(0,T_sim)
  bckgrnd_infections_staff <- rep(0,T_sim)
  bckgrnd_cases_staff <- rep(0,T_sim)
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
    # Transition 1: S -> E (infection)
    #
    # Transition 2: E -> I_s1, I_c1 (progression to early infectious stage)
    #
    # Transition 3: I_s1 -> I_s2, I_c1 -> I_c2 (progression to late infectious stage)
    #
    # Transition 4: I_s2, I_c2 -> R (recovery)
    #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Update states of individuals based on transmission on previous day
    list[TrueState,DayTrueState,DayObsState,WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,NewInfection,Background,DaysPCRpos,infections_res[t],cases_res[t],infections_staff[t],cases_staff[t],Hospitalised,Dead,bckgrnd_infections_res[t],bckgrnd_cases_res[t],bckgrnd_infections_staff[t],bckgrnd_cases_staff[t]] <- iterate(TrueState,Present,DayTrueState,DayObsState,DaysSinceInfctn,DaysSinceInfctsnss,beta,w,h,alpha,epsilon,N_pop,WaitingTime,r_E,p_E,e0ind,Risk,p_s,r_p,p_p,NewInfection,Background,DaysPCRpos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,Resident,hospitalisation,Hospitalised,Dead,p_h,p_ICU,p_d)
    
    # Store true states of individuals on day t
    state[,t] <- TrueState
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 2: REMOVAL OF INDIVIDUALS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    # Remove all individuals who tested PCR positive on the previous day
    Present[c(sx_indvdls_removed,PCRpos_removed)] <- F
    Present[setdiff(sx_indvdls_isolated,PCRpos_removed)] <- T
    DayRemoved[c(sx_indvdls_removed,PCRpos_removed)] <- t

    if (2 %in% interventions){
      idx <- (Resident==1)
      list[Present[idx],,Numbers_add1] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==1,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add2] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==2,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add3] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==3,sample(-5:5,1),0)
      list[Present[idx],,Numbers_add4] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==4,sample(-5:5,1),0)
    }
    
    presence[,t] <- Present
    sx_indvdls_removed <- integer() # reset list of symptomatic individuals removed to empty
    PCRpos_removed <- integer() # reset list of PCR positive individuals removed to empty
    
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
      for (i in 1:7){
        list[PCRtests,Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_on_entry_update(i,TrueState,Number,Numbers_add,entry_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed) 
      }
    }
        
    # Active symptom screening
    if (1 %in% interventions){
      HxSx_prev <- HxSx 
      for (i in 1:7){
        list[HxSx,Tested,DayTested,PCRtests,PCRtestsWeek,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- sx_screening_update(i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,min_days_btw_tests,sens_sx[i],spec_sx[i],HxSx,sx_pos_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
      }
      sx_indvdls_isolated <- which(HxSx & !HxSx_prev) # individuals who screen positive for symptoms on this day 
      Present[sx_indvdls_isolated] <- F # immediately isolate individuals who screen symptom positive until test results are returned the next day
    }
    
    # Passive symptom screening
    if (4 %in% interventions){
      idx_sx <- which(TrueState==6 & ObsState==1 & Present & DayTrueState==2) # severe symptomatic individuals present in shelter and not known to be infected self present for PCR testing on 3rd day of symptoms
      PCRtests[idx_sx] <- PCRtests[idx_sx] + 1
      PCRtestsWeek[idx_sx] <- PCRtestsWeek[idx_sx] + 1
      list[Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,DayTested,idx_sx,6,spec[6],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
    }
    
    # Routine PCR testing
    if (3 %in% interventions){
      if (t %in% testing_days){
        for (i in 1:7){
          list[PCRtests,PCRtestsWeek,Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- routine_PCR_testing_update(i,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,min_days_btw_tests,routine_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
        }        
      }
    }
    
    # Routine PCR testing of staff only
    if (7 %in% interventions){
      if (t %in% testing_days){
        for (i in 1:7){
          list[PCRtests,PCRtestsWeek,Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- routine_PCR_testing_update(i,TrueState,ObsState,(Present & Resident==0),PCRtestsWeek,max_PCR_tests_per_week,min_days_btw_tests,routine_PCR_test_compliance,PCRtests,Tested,DayTested,spec[i],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
        }        
      }
    }
  }
  
  # Make data frame of final values of variables for all individuals (columns as above)
  sim_pop <- data.frame(cbind(Number, Resident, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, Dead, ObsState, DayObsState, Tested, DayTested, PCRtests, PCRtestsWeek, HxPCR, DayRemoved, HxAb, HxSx, NewInfection, Background))
  return(res=list(infections_res=infections_res,cases_res=cases_res,infections_staff=infections_staff,cases_staff=cases_staff,bckgrnd_infections_res=bckgrnd_infections_res,bckgrnd_cases_res=bckgrnd_cases_res,bckgrnd_infections_staff=bckgrnd_infections_staff,bckgrnd_cases_staff=bckgrnd_cases_staff,PCRpos=PCRpos,sim_pop=sim_pop,state=state,presence=presence))
}
