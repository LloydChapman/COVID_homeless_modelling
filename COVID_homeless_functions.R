# PCR sensitivity as a function of time since start of infectiousness
sens <- function(x,fit,fit_extrap,max_days_PCR_pos){
  # res <- rep(1,length(x))
  # res <- rep(0.75,length(x))
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

# Function for updating presence in shelter
update_present <- function(Present,HxPCR,ind,x,y){
  if (y<x){
    ind_present <- which(ind & Present)
    Numbers_rm <- resample(ind_present,min(x-y,length(ind_present))) # remove number of individuals to match register if possible or all present in chosen risk group if not
    Present[Numbers_rm] <- F
    Numbers_add <- vector()
  } else if (y>x){
    # print(x)
    # print(y)
    # # print(Risk[ind][1])
    # print(length(which(ind & !Present & !HxPCR)))
    ind_absent_noHxPCR <- which(ind & !Present & !HxPCR) # individuals that have tested PCR positive cannot return to shelter
    Numbers_add <- resample(ind_absent_noHxPCR,min(y-x,length(ind_absent_noHxPCR))) # add number of individuals to match register if possible or all those in chosen risk group absent without positive PCR result if not
    Present[Numbers_add] <- T
    Numbers_rm <- vector()
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
routine_PCR_testing_update <- function(t,testing_days,state,TrueState,ObsState,Present,PCRtestsWeek,max_PCR_tests_per_week,routine_PCR_test_compliance,PCRtests,Tested,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,DayObsState,HxPCR,PCRpos,PCRpos_removed){
  # Routine PCR testing on testing days
  if (t %in% testing_days){
    idx <- which(TrueState==state & ObsState==1 & Present & PCRtestsWeek<max_PCR_tests_per_week & runif(length(TrueState))<routine_PCR_test_compliance) # individuals present in the shelter not known to be infected who have had less than max_PCR_tests_per_week PCR tests in the current week and comply
    PCRtests[idx] <- PCRtests[idx] + 1
    PCRtestsWeek[idx] <- PCRtestsWeek[idx] + 1
    list[Tested,ObsState,DayObsState,HxPCR,PCRpos,PCRpos_removed] <- PCR_testing_update(Tested,idx,state,spec,DaysSinceInfctsnss,fit,fit_extrap,max_days_PCR_pos,DayTrueState,DaysPCRpos,WaitingTime,ObsState,DayObsState,HxPCR,PCRpos,t,PCRpos_removed)
  }
  return(list(PCRtests=PCRtests,PCRtestsWeek=PCRtestsWeek,Tested=Tested,ObsState=ObsState,DayObsState=DayObsState,HxPCR=HxPCR,PCRpos=PCRpos,PCRpos_removed=PCRpos_removed))
}