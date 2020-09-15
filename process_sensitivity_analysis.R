process_sensitivity_analysis <- function(fnm1,run_nm,dir){
  load(fnm1)
  range_by_intvntn <- t(apply(prob_outbreak_averted,2,range))
  colnames(range_by_intvntn) <- c("min","max")
  write.csv(range_by_intvntn,paste0(dir,"range_prob_outbreak_averted",run_nm,".csv"))
  
  probs <- 0.5 #c(0.5,0.025,0.975)
  median_reduction_infections <- matrix(nrow = nrow(pars),ncol = ncol(prob_outbreak_averted))
  median_reduction_cases <- matrix(nrow = nrow(pars),ncol = ncol(prob_outbreak_averted))
  for (i in 1:nrow(pars)){
    load(paste0(dir,"intvntn_sim_output",run_nm,"_SA",i,".RData"))
    list[,median_reduction_infections[i,]] <- calc_perc_reduction(total_infections,probs)
    # median_reduction_infections[i,] <- tmp1[,1]
    list[,median_reduction_cases[i,]] <- calc_perc_reduction(total_cases,probs)
    # median_reduction_cases[i,] <- tmp2[,1]
  }
  range_reduction_infections <- t(apply(median_reduction_infections,2,range))
  colnames(range_reduction_infections) <- c("min","max")
  write.csv(range_reduction_infections,paste0(dir,"range_perc_reduction_infections",run_nm,".csv"))
  range_reduction_cases <- t(apply(median_reduction_cases,2,range))
  colnames(range_reduction_cases) <- c("min","max")
  write.csv(range_reduction_cases,paste0(dir,"range_perc_reduction_cases",run_nm,".csv"))
}