process_calibration <- function(fnm1,fnm2,run_nm,end_date,date_first_case,N_pop){
  res <- read.csv(fnm1,stringsAsFactors = F)
  # Calculate time before first case identified that infected individuals first entered shelter
  res[,3] <- as.numeric(date_first_case - (end_date - res[,3]) + 1)
  
  # Calculate median and IQR of fitted parameters
  q_par <- data.frame(apply(res,2,function(x) quantile(x,probs = c(0.5,0.025,0.975))))

  # Calculate estimated cumulative infection incidence
  q_cum_inc <- calc_cum_infection_inc(fnm2,N_pop)
  
  # Save results 
  write.csv(cbind(q_par,cum_inc=q_cum_inc),paste0("par_ests",run_nm,".csv"),row.names = F)
  
}
