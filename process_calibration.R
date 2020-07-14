process_calibration <- function(fnm,run_nm,end_date){
  res <- read.csv(fnm,stringsAsFactors = F)
  res[,3] <- as.numeric(end_date - res[,3])
  
  q_par <- data.frame(apply(res,2,function(x) quantile(x,probs = c(0.5,0.025,0.975))))
  q_par[,3] <- as.Date(q_par[,3],origin = "1970-01-01")
  
  write.csv(q_par,paste0("par_ests",run_nm,".csv"),row.names = F)
}
