COVID_homeless_model<-function(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,alpha,r_sx,p_sx,
                               p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,
                               hospitalisation,sens,spec,testing_days,N_tested,sx_testing_days,N_sx_tested,
                               CCMS_data,Number,Resident,Present,Risk,Age,e0ind,TrueState,DayTrueState,
                               WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos,shelter){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Population microsimulation 
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Note: Each person has a "true state" (TrueState): susceptible, exposed, early infectious (subclinical/clinical), late infectious (subclinical/clinical), recovered; and an "observed state" (ObsState): never infected or ever infected
  
  # Initialise variables
  Hospitalised <- rep(F,N_pop)
  Dead <- rep(F,N_pop)
  ObsState <- rep(1,N_pop)
  DayObsState <- rep(0,N_pop)
  Tested <- rep(0,N_pop) # ever PCR tested
  DayTested <- rep(NA,N_pop) # day individual last PCR tested
  HxPCR <- rep(F,N_pop) # has tested positive on PCR - these individuals are removed and not allowed to return
  DayRemoved <- rep(NA,N_pop) # day individual removed following PCR positive test 
  HxAb <- rep(F,N_pop) # has tested positive for antibodies
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
  #   13- Observed State: Uninfected (1); Infected (2)
  #   14- Days elapsed in OBSERVED state
  #   15- Tested (1), Untested (0) 
  #   16- Day last PCR tested
  #   17- History of positive PCR test: Yes (T), No (F)
  #   18- Day removed from shelter following positive PCR test
  #   19- History of positive Ab test (serology): Yes (T), No (F)
  #   20- True new infection
  #   21- Infected outside shelter: Yes (1), No (0)
  sim_pop0 <- data.frame(cbind(Number, Resident, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, ObsState, DayObsState, Tested, DayTested, HxPCR, DayRemoved, HxAb, NewInfection, Background))
  
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
  PCRpos_sx_testing <- rep(0,length(sx_testing_days))
  PCRpos <- rep(0,length(testing_days))
  state <- matrix(NA, nrow = N_pop, ncol = T_sim)
  state[,1] <- TrueState
  presence <- matrix(F, nrow = N_pop, ncol = T_sim)
  presence[,1] <- Present
  
  # Initialise vectors for storing individuals removed due to testing PCR positive
  sx_indvdls_removed <- integer()
  PCRpos_removed <- integer()
   
  for (t in 2:T_sim) {
    
    # print(t)
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # STEP 1: TRANSMISSION MODEL # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
    list[TrueState,DayTrueState,DayObsState,WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,NewInfection,Background,DaysPCRpos,infections_res[t],cases_res[t],infections_staff[t],cases_staff[t]] <- iterate(TrueState,Present,DayTrueState,DayObsState,DaysSinceInfctn,DaysSinceInfctsnss,beta,w,h,alpha,epsilon,N_pop,WaitingTime,r_E,p_E,e0ind,Risk,p_s,r_p,p_p,NewInfection,Background,DaysPCRpos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,Resident,hospitalisation,Hospitalised,Dead,p_h,p_ICU,p_d)
    
    # Store true states of individuals on day t
    state[,t] <- TrueState
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # STEP 2: REMOVAL OF INDIVIDUALS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    # Remove all individuals who tested PCR positive on the previous day
    Present[c(sx_indvdls_removed,PCRpos_removed)] <- F
    DayRemoved[c(sx_indvdls_removed,PCRpos_removed)] <- t
    # Remove individuals at random according to daily numbers in each risk category in CCMS data
    if (!is.null(CCMS_data)){
      idx <- which(Resident==1)
      list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==1,CCMS_data$Total_Low_Risk[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed %in% idx],PCRpos_removed[PCRpos_removed %in% idx])]==1),CCMS_data$Total_Low_Risk[t])
      list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==2,CCMS_data$Hi_Risk_60_only[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed %in% idx],PCRpos_removed[PCRpos_removed %in% idx])]==2),CCMS_data$Hi_Risk_60_only[t])
      list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==3,CCMS_data$Hi_Risk_Dx_only[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed %in% idx],PCRpos_removed[PCRpos_removed %in% idx])]==3),CCMS_data$Hi_Risk_Dx_only[t])
      list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==4,CCMS_data$Hi_Risk_Both_Age_Dx[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed %in% idx],PCRpos_removed[PCRpos_removed %in% idx])]==4),CCMS_data$Hi_Risk_Both_Age_Dx[t])
    }
    
    presence[,t] <- Present
    sx_indvdls_removed <- integer() # reset list of symptomatic individuals removed to empty
    PCRpos_removed <- integer() # reset list of PCR positive individuals removed to empty
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # STEP 3: DETECTION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if (t %in% sx_testing_days){ # Step 3.1: Selects for time steps where testing of cases with clinical symptoms takes place
      # Testing of TRUE clinical late infectious individuals
      i_c2 <- which((TrueState==6) & (ObsState==1) & (Tested==0) & Present) # Individuals have to be untested and present 
      if (length(i_c2)>0){
        sx_indvdls_tested <- resample(i_c2,min(N_sx_tested[sx_testing_days==t],length(i_c2)),prob=w[i_c2])
        list[Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos_sx_testing,sx_indvdls_removed] <- PCR_testing_update(Tested,DayTested,sx_indvdls_tested,6,spec[6],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos_sx_testing,which(sx_testing_days==t),sx_indvdls_removed)
      }
    }
    
    if (t %in% testing_days) { # Step 3.2: Selects for time steps where mass testing takes place
  
      # Allow repeat testing of individuals for Seattle shelters (which had 2 mass testings) but not for other shelters
      if (shelter %in% c("Seattle","Seattle_A","Seattle_B","Seattle_C")){
        idx1 <- which(Present) # & (Resident==1)
      } else {
        idx1 <- which((Tested==0) & Present) # & (Resident==1)  
      }
      if (length(idx1)>0){
        indvdls_tested <- resample(idx1,min(N_tested[testing_days==t],length(idx1))) #,prob=w[idx1] # randomly pick untested individuals to test, residents and staff have same chance of being tested
        for (i in 1:7){
          idx2 <- which((TrueState==i) & (ObsState==1) & (Number %in% indvdls_tested))
          list[Tested,DayTested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,DayTested,idx2,i,spec[i],DaysSinceInfctsnss,sens,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,which(testing_days==t),PCRpos_removed)
        }        
      }
    }
    
    # Overwrite individuals removed on 4/10 in MSC South with negatives from testing on 4/9 as negatives were removed on last day
    if (shelter=="SF"){
      if (t==testing_days[length(testing_days)-1]){
        PCRpos_removed <- setdiff(indvdls_tested,PCRpos_removed)
      }      
    }
}

# Make data frame of final values of variables for all individuals (columns as above)
sim_pop <- data.frame(cbind(Number, Resident, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, ObsState, DayObsState, Tested, DayTested, HxPCR, DayRemoved, HxAb, NewInfection, Background))
return(res=list(infections_res=infections_res,cases_res=cases_res,infections_staff=infections_staff,cases_staff=cases_staff,PCRpos_sx_testing=PCRpos_sx_testing,PCRpos=PCRpos,sim_pop=sim_pop,state=state,presence=presence))

}
