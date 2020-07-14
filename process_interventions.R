process_interventions <- function(fnm,run_nm){
  load(fnm)
  
  calc_perc_reduction <- function(y){
    perc_reduction_y <- 100*apply(y[,2:ncol(y)],2,function(x){(y[,1]-x)/y[,1]})
    perc_reduction_y[y[,1]==0,] <- 0
    return(perc_reduction_y)
  }
  
  # Function to calculate percentage of outbreaks "averted": number of 
  # simulations with no new infections for each intervention strategy out of 
  # number of simulations with no interventions in which there was an outbreak
  calc_perc_outbreaks_averted <- function(y){
    perc_outbreaks_averted <- 100*colSums(y[,2:ncol(y)]<3)/sum(y[,1]>=3)
    return(perc_outbreaks_averted)
  }
  
  perc_reduction_infections <- calc_perc_reduction(total_infections)
  # total_infections_sort <- apply(total_infections,2,sort)
  # perc_reduction_infections <- calc_perc_reduction(total_infections_sort)
  q_reduction_infections <- t(apply(perc_reduction_infections,2,function(x){quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)}))
  # rownames(q_reduction_infections) <- sapply(1:nrow(q_reduction_infections),function(x) paste0("Strategy ",x))
  write.table(q_reduction_infections,paste0("perc_reduction_infections",run_nm,".csv"),sep = ",",col.names = NA)
  
  perc_outbreaks_averted <- calc_perc_outbreaks_averted(total_infections)
  write.table(perc_outbreaks_averted,paste0("perc_outbreaks_averted",run_nm,".csv"),sep = ",",col.names = NA)
  
  perc_reduction_cases <- calc_perc_reduction(total_cases)
  # total_cases_sort <- apply(total_cases,2,sort)
  # perc_reduction_cases <- calc_perc_reduction(total_cases_sort)
  q_reduction_cases <- t(apply(perc_reduction_cases,2,function(x){quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)}))
  # rownames(q_reduction_cases) <- sapply(1:nrow(q_reduction_cases),function(x) paste0("Strategy ",x))
  write.table(q_reduction_cases,paste0("perc_reduction_cases",run_nm,".csv"),sep = ",",col.names = NA)
}
