process_interventions <- function(fnm,run_nm,T_sim){
  load(fnm)
  
  probs <- c(0.5,0.025,0.975)

  # For each intervention strategy:
  # Calculate percentage reduction in total number of infections
  list[perc_reduction_infections,q_reduction_infections] <- calc_perc_reduction(total_infections,probs)
  write.table(q_reduction_infections,paste0("perc_reduction_infections",run_nm,".csv"),sep = ",",col.names = NA)

  # Calculate percentage reduction in total number of clinical cases
  list[perc_reduction_cases,q_reduction_cases] <- calc_perc_reduction(total_cases,probs)
  write.table(q_reduction_cases,paste0("perc_reduction_cases",run_nm,".csv"),sep = ",",col.names = NA)

  # Calculate percentage reduction in total number of hospitalisations
  list[perc_reduction_hospitalisations,q_reduction_hospitalisations] <- calc_perc_reduction(total_hospitalisations,probs)
  write.table(q_reduction_hospitalisations,paste0("perc_reduction_hospitalisations",run_nm,".csv"),sep = ",",col.names = NA)

  # Calculate median and IQR of reduction in total number of hospitalisations
  q_hospitalisations <- t(apply(total_hospitalisations,2,function(x) quantile(x,probs=c(0.5,0.25,0.75))))
  write.table(q_hospitalisations,paste0("hospitalisations",run_nm,".csv"),sep = ",",col.names = NA)

  # Calculate percentage reduction in total number of deaths
  list[perc_reduction_deaths,q_reduction_deaths] <- calc_perc_reduction(total_deaths,probs)
  write.table(q_reduction_deaths,paste0("perc_reduction_deaths",run_nm,".csv"),sep = ",",col.names = NA)

  # Calculate median and IQR of reduction in total number of deaths
  q_deaths <- t(apply(total_deaths,2,function(x) quantile(x,probs=c(0.5,0.25,0.75))))
  write.table(q_deaths,paste0("deaths",run_nm,".csv"),sep = ",",col.names = NA)

  # Save results
  save(perc_reduction_infections,perc_reduction_cases,perc_reduction_hospitalisations,perc_reduction_deaths,q_reduction_infections,q_reduction_cases,q_reduction_hospitalisations,q_reduction_deaths,file=paste0("perc_reductions",run_nm,".RData"))

  # Calculate probability of averting outbreak
  prob_outbreak_averted <- calc_prob_outbreak_averted(infections,bckgrnd_infections)
  write.table(prob_outbreak_averted,paste0("prob_outbreak_averted",run_nm,".csv"),sep = ",",col.names = NA)
  
  # Calculate number of PCR tests used
  PCR_tests_used <- calc_PCR_tests_used(PCRtests,T_sim)
  write.csv(PCR_tests_used,paste0("PCR_tests_used",run_nm,".csv"),row.names = F)
}

# Function to combine results with same background infection rate
combine_results <- function(fnm,nintvntns,run_nms,i){
  df <- data.frame(matrix(nrow = nintvntns,ncol = 0))
  for (j in 1:length(run_nms)){
    tmp <- read.csv(paste0(fnm,run_nms[j],"_epsilon",i-1,".csv"),stringsAsFactors = F)
    df <- cbind(df,tmp[,2])
  }
  names(df) <- rep(fnm,length(run_nms))
  write.csv(df,paste0(fnm,"_epsilon",i-1,".csv"))
}
