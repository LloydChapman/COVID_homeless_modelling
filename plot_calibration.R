plot_calibration <- function(fn1,fn2,str,D_S,D_T,testing_dates,tol1,tol2){
  # Load results
  res <- read.csv(fn1,stringsAsFactors = F)
  beta_pstr <- res$V1
  load(fn2)
  
  # Plot output
  pdf(paste0("beta_posterior",str,".pdf"),width = 5, height = 4)
  hist(beta_pstr,breaks = 20,freq = F,xlab = "beta",main = "Posterior distribution of beta")
  dev.off()
  
  # pdf("num_PCRpos_vs_beta",str,".pdf",width = 5,height = 5)
  # plot(beta_pstr,PCRpos_mat,xlab = "beta",ylab = "No. PCR positive")
  # dev.off()
  
  pdf(paste0("onsets_vs_time",str,".pdf"),width = 6,height = 5)
  new_infections_melt <- melt(new_infections_mat)
  ggplot(new_infections_melt,aes(x=Var2,y=value,col=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. of onsets") + theme(legend.position = "none")
  dev.off()
  
  # Reconstruct the number of PCR positive individuals over time from individual-level simulation data frame
  getEventTime <- function(state,i){apply(state,c(1,3),function(x) which(x %in% i)[1])}
  
  tE <- getEventTime(state_arr,2)
  tI_m_p <- getEventTime(state_arr,3)
  tI_s_p <- getEventTime(state_arr,4)
  tI_m_sx <- getEventTime(state_arr,5)
  tI_s_sx <- getEventTime(state_arr,6)
  tR <- getEventTime(state_arr,7)
  
  StartPresx <- pmin(tI_m_p,tI_s_p,na.rm = T)
  Onset <- pmin(tI_m_sx,tI_s_sx,na.rm = T)
  EndPCRpos <- matrix(nrow = nrow(Onset), ncol = ncol(Onset))
  for (i in 1:N){
    EndPCRpos[,i] <- Onset[,i] + sim_pop_list[[i]]$DaysPCRpos
  }
  
  num_infctns <- apply(tE,2,function(y) sapply(1:T_sim,function(x) sum(y==x,na.rm = T)))
  num_infctns_melt <- melt(num_infctns)
  pdf(paste0("new_infections",str,".pdf"),width = 5,height = 4)
  ggplot() + geom_line(aes(x=Var1,y=value,group=Var2,col=as.factor(Var2)),num_infctns_melt) + xlab("Day") + ylab("No. new infections") + theme(legend.position = "none") #,group=Var2,col=as.factor(Var2)
  dev.off()
  
  detectable_viral_load_mat <- matrix(nrow = N, ncol = T_sim)
  for (i in 1:T_sim){
    detectable_viral_load_mat[,i] <- colSums(StartPresx<=i & EndPCRpos>i,na.rm = T)
  }
  
  detectable_viral_load_melt <- melt(detectable_viral_load_mat)
  pdf(paste0("num_detectabe_viral_load_vs_time",str,".pdf"),width = 5, height = 4)
  ggplot() + geom_line(aes(x=Var2,y=value,col=as.factor(Var1)),detectable_viral_load_melt) + geom_line(aes(x=day,y=PCRpos),data.frame(day=c(12,12),PCRpos=range(PCRpos_mat))) + geom_point(aes(x=12,y=D_T),col="black",size=2) + xlab("Day") + ylab("No. with detectable viral load") + theme(legend.position = "none")
  dev.off()
  
  testing_dates <- c(as.Date(c("4/4/2020","4/5/2020","4/6/2020","4/7/2020"),format="%m/%d/%Y"),as.Date(c("4/8/2020","4/9/2020"),format="%m/%d/%Y"))
  D <- c(D_S,D_T)
  PCRpos_mat1 <- data.frame(cbind(PCRpos_sx_testing_mat,PCRpos_mat))
  names(PCRpos_mat1) <- testing_dates
  PCRpos_mat1$sim <- row.names(PCRpos_mat1)
  PCRpos_melt1 <- melt(PCRpos_mat1,id.vars="sim")
  PCRpos_melt1$variable <- as.Date(PCRpos_melt1$variable)
  pdf(paste0("calibration",str,".pdf"),width = 6, height = 5)
  ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_melt1,col="#999999") + geom_point(aes(x=date,y=PCRpos),data.frame(date=testing_dates,PCRpos=D),col="black",size=2) + geom_line(aes(x=date,y=CI),data.frame(date=rep(testing_dates[5],2),CI=c(D_T[1]-tol1,D_T[1]+tol1)),size=1) + geom_line(aes(x=date,y=CI,size=1),data.frame(date=rep(testing_dates[6],2),CI=c(D_T[2]-tol2,D_T[2]+tol2)),size=1) + xlab("Date") + ylab("No. PCR positive") + theme(legend.position = "none") + theme_classic()
  dev.off()
  
  num_PCRpos_staff <- rep(0,length(sim_pop_list))
  for (i in 1:length(sim_pop_list)){
    num_PCRpos_staff[i] <- sum(sim_pop_list[[i]]$HxPCR[sim_pop_list[[i]]$Resident==0],na.rm = T) 
  }
  pdf(paste0("num_PCRpos_staff",str,".pdf"),width = 5, height = 4)
  hist(num_PCRpos_staff,breaks=seq(min(num_PCRpos_staff)-0.5,max(num_PCRpos_staff)+0.5),xlab="No. PCR positive staff",main="Distribution of no. of PCR positive staff")
  dev.off()
}
