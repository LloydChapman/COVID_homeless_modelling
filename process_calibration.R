process_calibration <- function(fnm1,fnm2,run_nm,end_date,date_first_case,N_pop){
  res <- read.csv(fnm1,stringsAsFactors = F)
  res[,3] <- as.numeric(date_first_case - (end_date - res[,3]) + 1)
  
  q_par <- data.frame(apply(res,2,function(x) quantile(x,probs = c(0.5,0.025,0.975))))
  # q_par[,3] <- as.Date(q_par[,3],origin = "1970-01-01")

  q_cum_inc <- calc_cum_infection_inc(fnm2,N_pop)
    
  write.csv(cbind(q_par,cum_inc=q_cum_inc),paste0("par_ests",run_nm,".csv"),row.names = F)
  
}
