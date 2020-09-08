process_interventions <- function(fnm,run_nm,T_sim){
  load(fnm)
  
  probs <- c(0.5,0.025,0.975)

  list[perc_reduction_infections,q_reduction_infections] <- calc_perc_reduction(total_infections,probs)
  write.table(q_reduction_infections,paste0("perc_reduction_infections",run_nm,".csv"),sep = ",",col.names = NA)
  # # perc_reduction_infections <- calc_perc_reduction(total_infections)
  # # total_infections_sort <- apply(total_infections,2,sort)
  # # perc_reduction_infections <- calc_perc_reduction(total_infections_sort)
  # q_reduction_infections <- t(apply(perc_reduction_infections,2,function(x){quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)}))
  # # rownames(q_reduction_infections) <- sapply(1:nrow(q_reduction_infections),function(x) paste0("Strategy ",x))
  # write.table(q_reduction_infections,paste0("perc_reduction_infections",run_nm,".csv"),sep = ",",col.names = NA)

  list[perc_reduction_cases,q_reduction_cases] <- calc_perc_reduction(total_cases,probs)
  write.table(q_reduction_cases,paste0("perc_reduction_cases",run_nm,".csv"),sep = ",",col.names = NA)
  # perc_reduction_cases <- calc_perc_reduction(total_cases)
  # # total_cases_sort <- apply(total_cases,2,sort)
  # # perc_reduction_cases <- calc_perc_reduction(total_cases_sort)
  # q_reduction_cases <- t(apply(perc_reduction_cases,2,function(x){quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)}))
  # # rownames(q_reduction_cases) <- sapply(1:nrow(q_reduction_cases),function(x) paste0("Strategy ",x))
  # write.table(q_reduction_cases,paste0("perc_reduction_cases",run_nm,".csv"),sep = ",",col.names = NA)

  list[perc_reduction_hospitalisations,q_reduction_hospitalisations] <- calc_perc_reduction(total_hospitalisations,probs)
  write.table(q_reduction_hospitalisations,paste0("perc_reduction_hospitalisations",run_nm,".csv"),sep = ",",col.names = NA)
  # perc_reduction_hospitalisations <- calc_perc_reduction(total_hospitalisations)
  # q_reduction_hospitalisations <- t(apply(perc_reduction_hospitalisations,2,function(x){quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)}))
  # write.table(q_reduction_hospitalisations,paste0("perc_reduction_hospitalisations",run_nm,".csv"),sep = ",",col.names = NA)

  q_hospitalisations <- t(apply(total_hospitalisations,2,function(x) quantile(x,probs=c(0.5,0.25,0.75))))
  write.table(q_hospitalisations,paste0("hospitalisations",run_nm,".csv"),sep = ",",col.names = NA)

  list[perc_reduction_deaths,q_reduction_deaths] <- calc_perc_reduction(total_deaths,probs)
  write.table(q_reduction_deaths,paste0("perc_reduction_deaths",run_nm,".csv"),sep = ",",col.names = NA)

  q_deaths <- t(apply(total_deaths,2,function(x) quantile(x,probs=c(0.5,0.25,0.75))))
  write.table(q_deaths,paste0("deaths",run_nm,".csv"),sep = ",",col.names = NA)

  save(perc_reduction_infections,perc_reduction_cases,perc_reduction_hospitalisations,perc_reduction_deaths,q_reduction_infections,q_reduction_cases,q_reduction_hospitalisations,q_reduction_deaths,file=paste0("perc_reductions",run_nm,".RData"))

  # perc_outbreaks_averted <- calc_perc_outbreaks_averted(total_infections)
  # write.table(perc_outbreaks_averted,paste0("perc_outbreaks_averted",run_nm,".csv"),sep = ",",col.names = NA)

  prob_outbreak_averted <- calc_prob_outbreak_averted(infections,bckgrnd_infections)
  write.table(prob_outbreak_averted,paste0("prob_outbreak_averted",run_nm,".csv"),sep = ",",col.names = NA)
  
  PCR_tests_used <- calc_PCR_tests_used(PCRtests,T_sim)
  write.csv(PCR_tests_used,paste0("PCR_tests_used",run_nm,".csv"),row.names = F)
}

combine_results <- function(fnm,nintvntns,run_nms,i){
  df <- data.frame(matrix(nrow = nintvntns,ncol = 0))
  for (j in 1:length(run_nms)){
    tmp <- read.csv(paste0(fnm,run_nms[j],"_epsilon",i-1,".csv"),stringsAsFactors = F)
    df <- cbind(df,tmp[,2])
  }
  names(df) <- rep(fnm,length(run_nms))
  write.csv(df,paste0(fnm,"_epsilon",i-1,".csv"))
}
