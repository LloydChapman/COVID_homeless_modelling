calc_R0 <- function(beta,w,Present,p_s,Risk,h,alpha,mu_p,mu_sx){
  R0 <- beta*sum(w*Present*(((1-p_s[Risk])*h+p_s[Risk])*(alpha*mu_p+mu_sx)))
  return(R0)
}

# PCR sensitivity as a function of time since start of infectiousness
sens <- function(x,fit,fit_extrap,max_days_PCR_pos){
  # res <- rep(1,length(x))
  # # res <- rep(0.75,length(x))
  res <- rep(0,length(x))
  x1 <- length(fit$fitted.values)
  idx <- (x>0 & x<=x1)
  if (any(idx)){
    res[idx] <- suppressWarnings(predict(fit,data.frame(x=x[idx])))
  }
  # idx1 <- which(x>x1 & x<=max_days_PCR_pos)
  # for (i in 1:length(idx1)){
  #   res[idx1[i]] <- max(fit_extrap$y[fit_extrap$x==x[idx1[i]]],0)
  # }
  idx1 <- (x>x1 & x<=max_days_PCR_pos)
  res[idx1] <- sapply(x[idx1],function(x) fit_extrap$y[fit_extrap$x==x])
  return(res)
}

resample <- function(x,...) x[sample.int(length(x),...)]

iterate <- function(TrueState,Present,DayTrueState,DayObsState,DaysSinceInfctn,DaysSinceInfctsnss,beta,w,h,alpha,epsilon,N_pop,WaitingTime,r_E,p_E,e0ind,Risk,p_s,r_p,p_p,NewInfection,DaysPCRpos,min_days_PCR_pos,max_days_PCR_pos,discrnorm,hospitalisation = F,Hospitalised = NULL,Dead = NULL,p_h = NULL,p_ICU = NULL,p_d = NULL){
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
  
  # Advance time by one day - [ ] Think carefully about whether this should be at start or end of loop
  DayTrueState <- DayTrueState + 1
  DayObsState <- DayObsState + 1
  
  # Increase days since infection by one day for infected individuals
  DaysSinceInfctn[e|i_m_p|i_s_p|i_m_sx|i_s_sx|r] <- DaysSinceInfctn[e|i_m_p|i_s_p|i_m_sx|i_s_sx|r] + 1
  # Increase days since start of infectiousness by one day for infected individuals
  DaysSinceInfctsnss[i_m_p|i_s_p|i_m_sx|i_s_sx|r] <- DaysSinceInfctsnss[i_m_p|i_s_p|i_m_sx|i_s_sx|r] + 1
  
  # Step 1: S -> E (infection)
  lambda <- beta*Present*w*(h*(alpha*I_m_p+I_m_sx)+alpha*I_s_p+I_s_sx) + epsilon # Only individuals present in the shelter can be infected by others in the shelter, those not present can be infected by infectious individuals in the general population
  # print(paste0("sum(lambda)=",sum(lambda)))
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
  
  # Count number of new infections
  infections <- sum(s2e)
  
  # Count number of new clinical cases
  cases <- sum(i_p2i_sx & i_s_p)
  
  return(list(TrueState=TrueState,DayTrueState=DayTrueState,DayObsState=DayObsState,WaitingTime=WaitingTime,DaysSinceInfctn=DaysSinceInfctn,DaysSinceInfctsnss=DaysSinceInfctsnss,NewInfection=NewInfection,DaysPCRpos=DaysPCRpos,infections=infections,cases=cases,Hospitalised=Hospitalised,Dead=Dead))

}

# Function for updating presence in shelter
presence_update <- function(Present,HxPCR,ind,x,y){
  if (y<x){
    ind_present <- which(ind & Present)
    Numbers_rm <- resample(ind_present,min(x-y,length(ind_present))) # remove number of individuals to match register if possible or all present in chosen risk group if not
    Present[Numbers_rm] <- F
    Numbers_add <- integer()
  } else if (y>x){
    # print(x)
    # print(y)
    # # print(Risk[ind][1])
    # print(length(which(ind & !Present & !HxPCR)))
    ind_absent_noHxPCR <- which(ind & !Present & !HxPCR) # individuals that have tested PCR positive cannot return to shelter
    Numbers_add <- resample(ind_absent_noHxPCR,min(y-x,length(ind_absent_noHxPCR))) # add number of individuals to match register if possible or all those in chosen risk group absent without positive PCR result if not
    Present[Numbers_add] <- T
    Numbers_rm <- integer()
  } else {
    Numbers_rm <- integer()
    Numbers_add <- integer()
  }
  return(list(Present=Present,Numbers_rm=Numbers_rm,Numbers_add=Numbers_add))
}

# Function for updates from PCR testing
PCR_testing_update <- function(Tested,idx,state,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed){
  Tested[idx] <- 1
  if (state %in% c(1,2)){
    ni2i <- (runif(length(idx))<1-spec)
  } else if (state %in% c(3,4)){
    ni2i <- (runif(length(idx))<sens(DaysSinceInfctsnss[idx],fit,fit_extrap,max_days_PCR_pos)) # Infected and test is positive
  } else if (state %in% c(5,6)){
    ni2i <- (DayTrueState[idx]<=DaysPCRpos[idx] & runif(length(idx))<sens(DaysSinceInfctsnss[idx],fit,fit_extrap,max_days_PCR_pos)) # Time since sx onset < duration of PCR positivity (viraemia) and test is positive
  } else {
    ni2i <- (WaitingTime[idx]+DayTrueState[idx]<=DaysPCRpos[idx] & runif(length(idx))<sens(DaysSinceInfctsnss[idx],fit,fit_extrap,max_days_PCR_pos)) # Time since sx onset = waiting time in sx state + time in recovered state < duration of PCR positivity and test is positive 
  }
  ObsState[idx] <- 1 + ni2i # update observed data if tested positive
  DayObsState[idx][ni2i] <- 0 # re-set day in observed state
  HxPCR[idx][ni2i] <- 1 # assign yes for hx positive PCR
  HxPCR[idx][!ni2i] <- 0 # assign no for hx positive PCR
  PCRpos[t] <- PCRpos[t] + sum(ni2i) # add number of PCR positives in this state to number PCR positive on day t
  PCRpos_removed <- c(PCRpos_removed,idx[ni2i]) # add PCR positives to list of PCR positive individuals to be removed the next day
  
  return(list(Tested=Tested,ObsState=ObsState,DayObsState=DayObsState,HxPCR=HxPCR,PCRpos=PCRpos,PCRpos_removed=PCRpos_removed))
}

# Function for updates resulting from symptom screening
sx_screening_update <- function(state,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,sens_sx,spec_sx,HxSx,PCRtests,Tested,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed){
  idx <- ((TrueState==state) & (ObsState==1) & Present & (PCRtestsWeek<max_PCR_tests_per_week))
  if (any(idx)){
    if (state %in% c(1,2,3,4,5,7)){
      sx_pos <- (runif(sum(idx))<1-spec_sx) # draw which individuals are positive for symptoms
    } else {
      sx_pos <- (runif(sum(idx))<sens_sx)
    }
    idx_pos <- which(idx)[sx_pos] # get indices of positives
    HxSx[idx_pos] <- T # mark as having symptoms
    
    # PCR test those positive for symptoms
    PCRtests[idx_pos] <- PCRtests[idx_pos] + 1 # increment counter for total number of PCR tests
    PCRtestsWeek[idx_pos] <- PCRtestsWeek[idx_pos] + 1 # increment counter for number of PCR tests this week
    list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx_pos,state,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
  }
  return(list(HxSx=HxSx,Tested=Tested,PCRtests=PCRtests,PCRtestsWeek=PCRtestsWeek,ObsState=ObsState,DayObsState=DayObsState,HxPCR=HxPCR,PCRpos=PCRpos,PCRpos_removed=PCRpos_removed))
}

# Function for updates resulting from PCR testing once upon entry
PCR_testing_on_entry_update <- function(state,TrueState,Number,Numbers_add,entry_PCR_test_compliance,PCRtests,Tested,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed){
  idx <- which((TrueState==state) & (Number %in% Numbers_add) & runif(length(TrueState))<entry_PCR_test_compliance) # test anybody entering shelter regardless of whether they've been tested previously and what the result was
  PCRtests[idx] <- PCRtests[idx] + 1
  list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx,state,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
  
  return(list(PCRtests=PCRtests,Tested=Tested,ObsState=ObsState,DayObsState=DayObsState,HxPCR=HxPCR,PCRpos=PCRpos,PCRpos_removed=PCRpos_removed))
}

# Function for updates resulting from routine PCR testing
routine_PCR_testing_update <- function(state,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,routine_PCR_test_compliance,PCRtests,Tested,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,t,PCRpos_removed){
  idx <- which(TrueState==state & ObsState==1 & Present & PCRtestsWeek<max_PCR_tests_per_week & runif(length(TrueState))<routine_PCR_test_compliance) # individuals present in the shelter not known to be infected who have had less than max_PCR_tests_per_week PCR tests in the current week and comply
  PCRtests[idx] <- PCRtests[idx] + 1
  PCRtestsWeek[idx] <- PCRtestsWeek[idx] + 1
  list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx,state,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
  return(list(PCRtests=PCRtests,PCRtestsWeek=PCRtestsWeek,Tested=Tested,ObsState=ObsState,DayObsState=DayObsState,HxPCR=HxPCR,PCRpos=PCRpos,PCRpos_removed=PCRpos_removed))
}