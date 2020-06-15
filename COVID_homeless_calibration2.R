# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Title: Comparison of public health strategies to reduce COVID-19 in homeless populations across the United States 
# Code author: Lloyd Chapman, Nathan C Lo, MD PhD (UCSF)
# Origin date: 4/21/20
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

COVID_homeless_model<-function(N_res,N_staff,N_pop,T_sim,w,beta,epsilon,r_E,p_E,p_s,h,r_p,p_p,alpha,r_sx,p_sx,
                               p_h,p_ICU,p_d,mean_days_PCR_pos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,
                               hospitalisation,fit,fit_extrap,spec,testing_days,N_tested,sx_testing_days,
                               CCMS_data,Number,Alive,Resident,Present,Risk,Age,TrueState,DayTrueState,
                               WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,DaysPCRpos){

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
  
  # Initialise variables
  Hospitalised <- rep(F,N_pop)
  Dead <- rep(F,N_pop)
  ObsState <- rep(1,N_pop)
  DayObsState <- rep(0,N_pop)
  Tested <- rep(0,N_pop)
  HxPCR <- rep(F,N_pop) # has tested positive on PCR - these individuals are removed and not allowed to return
  DayRemoved <- rep(NA,N_pop) # day individual removed following PCR positive test 
  HxAb <- rep(F,N_pop) # has tested positive for antibodies
  NewInfection <- rep(0,N_pop)
  
  sim_pop0 <- data.frame(cbind(Number, Resident, Alive, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, ObsState, DayObsState, Tested, HxPCR, DayRemoved, HxAb, NewInfection))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Decision trees
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Decision tree logic 
  
  # Time step: Daily 
  # Note: Each person has a "true state" (susceptible, exposed, presymptomatic infectious (mild/severe), symptomatic infectious (mild/severe), recovered) and an "observed state" (non-immune, infected or immune)
  
  # Create objects for storing simulation output
  infections <- rep(0,T_sim)
  cases <- rep(0,T_sim)
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
    list[TrueState,DayTrueState,DayObsState,WaitingTime,DaysSinceInfctn,DaysSinceInfctsnss,NewInfection,DaysPCRpos,infections[t],cases[t]] <- iterate(TrueState,Present,DayTrueState,DayObsState,DaysSinceInfctn,DaysSinceInfctsnss,beta,w,h,alpha,epsilon,N_pop,WaitingTime,r_E,p_E,e0ind,Risk,p_s,r_p,p_p,NewInfection,DaysPCRpos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation,Hospitalised,Dead,p_h,p_ICU,p_d)
    
    # Store true states of individuals on day t
    state[,t] <- TrueState
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 2: REMOVAL OF INDIVIDUALS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    # Remove all individuals who tested PCR positive on the previous day
    Present[c(sx_indvdls_removed,PCRpos_removed)] <- F
    DayRemoved[c(sx_indvdls_removed,PCRpos_removed)] <- t
    # Remove individuals at random according to daily numbers in each risk category in CCMS data
    # Create index for individuals who can be removed
    if (t<sx_testing_days[1]){ # if before first symptomatic testing date, first cases can't be removed
      idx <- setdiff(1:N_res,c(e0,i_s_p0))
    } else {
      idx <- 1:N_res
    }
    list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==1,CCMS_data$Total_Low_Risk[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==1),CCMS_data$Total_Low_Risk[t])
    list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==2,CCMS_data$Hi_Risk_60_only[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==2),CCMS_data$Hi_Risk_60_only[t])
    list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==3,CCMS_data$Hi_Risk_Dx_only[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==3),CCMS_data$Hi_Risk_Dx_only[t])
    list[Present[idx]] <- presence_update(Present[idx],HxPCR[idx],Risk[idx]==4,CCMS_data$Hi_Risk_Both_Age_Dx[t-1] - sum(Risk[c(sx_indvdls_removed[sx_indvdls_removed<=N_res],PCRpos_removed[PCRpos_removed<=N_res])]==4),CCMS_data$Hi_Risk_Both_Age_Dx[t])
    
    presence[,t] <- Present
    # print(c(sx_indvdls_removed,PCRpos_removed))
    sx_indvdls_removed <- integer() # reset list of symptomatic individuals removed to empty
    PCRpos_removed <- integer() # reset list of PCR positive individuals removed to empty
    # print(sum(Present[1:N_res]))
    # print(summary(as.factor(Risk[Resident==1 & Present])))
    # print(summary(as.factor(Risk[Resident==1 & !Present])))
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ALGORITHM 3: DETECTION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #
    # Focus on TRUE state
    #
    # Step 1: Is this person at work?
    #         - Yes: Present at work
    #         - No: Not present at work (due to having confirmed PCR infection req 14 day self-quarantine, hospitalization due to COVID)
    # Step 2: If YES (Step 1), What is the OBSERVED state (i.e. confirmed immune status)?
    #         - Non-immuune: No hx confirmed PCR/Ab positivity
    #         - Immune: hx confirmed PCR/Ab positivity
    # Step 3: If NON-IMMUNE (Step 2), is this a test day based on the frequency of routine testing?
    #         - Yes: Testing will take place, proceed to PCR testing (sens/spec ~ true state)
    #         - No: Testing will not take place
    # Step 4: YES (Step 3), is the tested person positive?
    #         - Yes: Confirmed new infection, will req 14 day self-quarantine during Step 1 the subsequent dayy.
    #         - No: Goes to work, no change
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
    if (t %in% sx_testing_days){ # Step 3.1: Selects for time steps where testing of cases with clinical symptoms takes place
      # Testing of TRUE sx infectious severe
      i_s_sx1 <- which((TrueState==6) & (ObsState==1) & (Tested==0) & Present) # Individuals have to be untested and present 
      # print(length(i_s_sx1))
      # print(DayTrueState[i_s_sx1])
      # print(sens(DaysSinceInfctsnss[i_s_sx1],fit,fit_extrap,max_days_PCR_pos))
      list[Tested,ObsState,DayObsState,HxPCR,PCRpos_sx_testing,sx_indvdls_removed] <- PCR_testing_update(Tested,i_s_sx1,6,spec[6],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos_sx_testing,which(sx_testing_days==t),sx_indvdls_removed)
      print(PCRpos_sx_testing)
      # print(sx_indvdls_removed)
    }
    
    if (t %in% testing_days) { # Step 3.2: Selects for time steps where mass testing takes place
      # STEP 1 (remove ObsState=2)
      # STEP 2 (remove ObsState=3)
      
      # Columns
      #   1- Person number
      #   2- Alive (1), Dead (0)
      #   3- Age
      #   4- Risk group: 1- very low; 2- low; 3- medium; 4- high - [ ] UPDATE defs
      #   5- True State: 0- dead; 1- susceptible; 2- exposed; 3- infected & asx; 4- infected & sx; 5- infected & hospitalized ; 6- recovered
      #   6- Day in TRUE state: (1,N)
      #   7- Observed State: 1- non-immune; 2- infected; 3- immune (defined post confirmed infection)
      #   8- Day in OBSERVED state: (1,N)
      #   9- Hx positive PCR: 0- no; 1- yes (detected positive PCR)
      #   10- Hx positive Ab (serology): 0- no; 1- yes
      #   11- True new infection
  
      indvdls_tested <- sample(which(Tested==0 & Present),N_tested[testing_days==t],prob=w[Tested==0 & Present]) # randomly pick untested individuals to test
      for (i in 1:7){
        idx1 <- which((TrueState==i) & (ObsState==1) & (Number %in% indvdls_tested))
        list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx1,i,spec[i],DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,which(testing_days==t),PCRpos_removed)
      }
    }
}

sim_pop <- data.frame(cbind(Number, Resident, Alive, Present, Age, Risk, TrueState, DayTrueState, WaitingTime, DaysSinceInfctn, DaysSinceInfctsnss, DaysPCRpos, Hospitalised, ObsState, DayObsState, Tested, HxPCR, DayRemoved, HxAb, NewInfection))
return(res=list(infections=infections,cases=cases,PCRpos_sx_testing=PCRpos_sx_testing,PCRpos=PCRpos,sim_pop=sim_pop,state=state,presence=presence))

}
