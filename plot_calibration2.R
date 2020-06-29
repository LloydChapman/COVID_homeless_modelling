plot_calibration <- function(fnm1,fnm2,str,lm.low,lm.upp,end_date,D_S,D_T,D_C=NULL,testing_dates,sx_testing_dates,tol){
  # Load results
  res <- read.csv(fnm1,stringsAsFactors = F)
  N <- dim(res)[1]
  npar <- dim(res)[2]
  load(fnm2)
  
  # Plot posteriors
  pdf(paste0("R0_posterior",str,".pdf"),width = 5, height = 4)
  hist(res[,1],breaks = seq(lm.low[1],lm.upp[1],by = 0.5),freq = F,xlab = "R0",ylab = "Density",main = "")
  dev.off()
  pdf(paste0("E0_posterior",str,".pdf"),width = 5, height = 4)
  hist(res[,2],breaks = (lm.low[2]-0.5):(lm.upp[2]+0.5),freq = F,xlab = "E0",ylab = "Density",main = "")
  dev.off()
  pdf(paste0("T_posterior",str,".pdf"),width = 5, height = 4)
  hist(res[,3],breaks = (lm.low[3]-0.5):(lm.upp[3]+0.5),freq = F,xlab = "T",ylab = "Density",main = "")
  dev.off()
  
  # Plot pairwise parameter correlation
  p <- vector("list",npar^2)
  binwdths <- c(0.5,1,1)
  # for (i in 1:npar){
  #   # print(i)
  #   for (j in 1:npar){
  #     # print(j)
  #     # print(npar*(i-1)+j)
  #     if (i==j){
  #       print(i)
  #       print(j)
  #       p[[npar*(i-1)+j]] <- ggplot(res[,j],aes(x=res[,j])) + geom_histogram() + theme_classic()
  #     } else {
  #       p[[npar*(i-1)+j]] <- ggplot(res[,c(i,j)],aes(x=res[,j],y=res[,i])) + geom_bin2d(binwidth=c(binwdths[j],binwdths[i])) + theme_classic()
  #     }
  #   }
  # }
  p1 <- ggplot(res,aes(x=R0)) + geom_histogram(aes(y=..density..),binwidth = binwdths[1]) + theme_classic()
  p2 <- ggplot(res,aes(x=E0,y=R0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[2],binwdths[1]),show.legend = F) + theme_classic()
  p3 <- ggplot(res,aes(x=T,y=R0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[3],binwdths[1]),show.legend = F) + theme_classic()
  p4 <- ggplot(res,aes(x=R0,y=E0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[1],binwdths[2]),show.legend = F) + theme_classic()
  p5 <- ggplot(res,aes(x=E0)) + geom_histogram(aes(y=..density..),binwidth = binwdths[2]) + theme_classic()
  p6 <- ggplot(res,aes(x=T,y=E0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[3],binwdths[2]),show.legend = F) + theme_classic()
  p7 <- ggplot(res,aes(x=R0,y=T)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[1],binwdths[3]),show.legend = F) + theme_classic()
  p8 <- ggplot(res,aes(x=E0,y=T)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[2],binwdths[3]),show.legend = F) + theme_classic()
  p9 <- ggplot(res,aes(x=T)) + geom_histogram(aes(y=..density..),binwidth = binwdths[3]) + theme_classic()
  pdf(paste0("parameter_correlation",str,".pdf"),width = 6,height = 6)
  print(grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3,ncol=3))
  dev.off()
  
  # ggpairs(data.frame(R0=res[,1],E0=as.factor(res[,2]),T=as.factor(res[,3])),lower=list(combo="box_no_facet"),diag = list(continuous="barDiag",discrete="barDiag"))
  
  # Daily numbers of infections among residents
  infections_mat <- matrix(0,nrow = N,ncol = lm.upp[3])
  for (i in 1:N){
    infections_mat[i,lm.upp[3]+((-res[i,3]+1):0)] <- infections_list[[i]]
  }
  colnames(infections_mat) <- as.character(seq.Date(end_date-lm.upp[3]+1,end_date,by=1))
  infections_melt <- melt(infections_mat)
  infections_melt$Var2 <- as.Date(infections_melt$Var2)
  pdf(paste0("new_infections",str,".pdf"),width = 5,height = 4)
  print(ggplot() + geom_line(aes(x=Var2,y=value,group=Var1,col=as.factor(Var1)),infections_melt) + xlab("Date") + ylab("No. new infections") + theme(legend.position = "none")) #,group=Var2,col=as.factor(Var2)
  dev.off()
  
  # Daily numbers of new symptomatic cases among residents
  cases_mat1 <- matrix(0,nrow = N,ncol = lm.upp[3])
  for (i in 1:N){
    cases_mat1[i,lm.upp[3]+((-res[i,3]+1):0)] <- cases_list[[i]]
  }
  colnames(cases_mat1) <- as.character(seq.Date(end_date-lm.upp[3]+1,end_date,by=1))
  cases_melt <- melt(cases_mat1)
  cases_melt$Var2 <- as.Date(cases_melt$Var2)
  pdf(paste0("new_cases",str,".pdf"),width = 5,height = 4)
  print(ggplot() + geom_line(aes(x=Var2,y=value,group=Var1,col=as.factor(Var1)),cases_melt) + xlab("Date") + ylab("No. new symptomatic cases") + theme(legend.position = "none")) #,group=Var2,col=as.factor(Var2)
  dev.off()
  
  # pdf(paste0("onsets_vs_time",str,".pdf"),width = 6,height = 5)
  # new_infections_melt <- melt(new_infections_mat)
  # ggplot(new_infections_melt,aes(x=Var2,y=value,col=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. of onsets") + theme(legend.position = "none")
  # dev.off()
  # 
  # # Reconstruct the number of PCR positive individuals over time from individual-level simulation data frame
  # getEventTime <- function(state,i){apply(state,c(1,3),function(x) which(x %in% i)[1])}
  # 
  # tE <- getEventTime(state_arr,2)
  # tI_m_p <- getEventTime(state_arr,3)
  # tI_s_p <- getEventTime(state_arr,4)
  # tI_m_sx <- getEventTime(state_arr,5)
  # tI_s_sx <- getEventTime(state_arr,6)
  # tR <- getEventTime(state_arr,7)
  # 
  # StartPresx <- pmin(tI_m_p,tI_s_p,na.rm = T)
  # Onset <- pmin(tI_m_sx,tI_s_sx,na.rm = T)
  # EndPCRpos <- matrix(nrow = nrow(Onset), ncol = ncol(Onset))
  # for (i in 1:N){
  #   EndPCRpos[,i] <- Onset[,i] + sim_pop_list[[i]]$DaysPCRpos
  # }
  # 
  # num_infctns <- apply(tE,2,function(y) sapply(1:T_sim,function(x) sum(y==x,na.rm = T)))
  # num_infctns_melt <- melt(num_infctns)
  # pdf(paste0("new_infections",str,".pdf"),width = 5,height = 4)
  # ggplot() + geom_line(aes(x=Var1,y=value,group=Var2,col=as.factor(Var2)),num_infctns_melt) + xlab("Day") + ylab("No. new infections") + theme(legend.position = "none") #,group=Var2,col=as.factor(Var2)
  # dev.off()
  # 
  # detectable_viral_load_mat <- matrix(nrow = N, ncol = T_sim)
  # for (i in 1:T_sim){
  #   detectable_viral_load_mat[,i] <- colSums(StartPresx<=i & EndPCRpos>i,na.rm = T)
  # }
  # 
  # detectable_viral_load_melt <- melt(detectable_viral_load_mat)
  # pdf(paste0("num_detectabe_viral_load_vs_time",str,".pdf"),width = 5, height = 4)
  # ggplot() + geom_line(aes(x=Var2,y=value,col=as.factor(Var1)),detectable_viral_load_melt) + geom_line(aes(x=day,y=PCRpos),data.frame(day=c(12,12),PCRpos=range(PCRpos_mat))) + geom_point(aes(x=12,y=D_T),col="black",size=2) + xlab("Day") + ylab("No. with detectable viral load") + theme(legend.position = "none")
  # dev.off()
  
  PCRpos_df <- data.frame(PCRpos_mat)
  names(PCRpos_df) <- testing_dates
  PCRpos_df$sim <- row.names(PCRpos_df)
  PCRpos_melt <- melt(PCRpos_df,id.vars="sim")
  PCRpos_melt$variable <- as.Date(PCRpos_melt$variable)
  
  # if (length(PCRpos_sx_testing_mat)==0){
  #   pdf(paste0("calibration",str,".pdf"),width = 6, height = 5)
  #   print(ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_melt,col="#999999") + 
  #           geom_point(aes(x=date,y=PCRpos),data.frame(date=testing_dates,PCRpos=D_T),col="black",size=2) + 
  #           geom_line(aes(x=date,y=CI,group=as.factor(date)),data.frame(date=rep(testing_dates,2),CI=pmax(c(D_T-tol,D_T+tol),0)),size=1) + 
  #           xlab("Date") + ylab("No. PCR positive") + theme(legend.position = "none") + theme_classic())
  #   dev.off()    
  # } else {
  #   PCRpos_sx_testing_df <- data.frame(PCRpos_sx_testing_mat)
  #   names(PCRpos_sx_testing_df) <- sx_testing_dates
  #   PCRpos_sx_testing_df$sim <- row.names(PCRpos_sx_testing_df)
  #   PCRpos_sx_testing_melt <- melt(PCRpos_sx_testing_df,id.vars="sim")
  #   PCRpos_sx_testing_melt$variable <- as.Date(PCRpos_sx_testing_melt$variable,origin="1970-01-01")
  #   
  #   pdf(paste0("calibration",str,".pdf"),width = 6, height = 5)
  #   print(ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_melt,col="#999999") + 
  #           geom_point(aes(x=date,y=PCRpos),data.frame(date=testing_dates,PCRpos=D_T),col="black",size=2) + 
  #           geom_line(aes(x=date,y=CI,group=as.factor(date)),data.frame(date=rep(testing_dates,2),CI=pmax(c(D_T-tol,D_T+tol),0)),size=1) + 
  #           geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_sx_testing_melt,col="blue",linetype="longdash") + 
  #           geom_point(aes(x=date,y=PCRpos),data.frame(date=sx_testing_dates,PCRpos=D_S),col="black",size=2,pch=15) + 
  #           xlab("Date") + ylab("No. PCR positive") + theme(legend.position = "none") + theme_classic())
  #   dev.off()    
  # }
  
  # Calibration plot with numbers of PCR positives
  p <- ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_melt,col="#999999") +
          geom_point(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_melt,col="#999999",pch=19) +
          geom_point(aes(x=date,y=PCRpos),data.frame(date=testing_dates,PCRpos=D_T),col="black",size=2) + 
          geom_line(aes(x=date,y=CI,group=as.factor(date)),data.frame(date=rep(testing_dates,2),CI=pmax(c(D_T-tol,D_T+tol),0)),size=1) + 
          xlab("Date") + ylab("No. PCR positive") + theme(legend.position = "none") + theme_classic()
  if (length(sx_testing_dates)!=0){
    PCRpos_sx_testing_df <- data.frame(PCRpos_sx_testing_mat)
    names(PCRpos_sx_testing_df) <- sx_testing_dates
    PCRpos_sx_testing_df$sim <- row.names(PCRpos_sx_testing_df)
    PCRpos_sx_testing_melt <- melt(PCRpos_sx_testing_df,id.vars="sim")
    PCRpos_sx_testing_melt$variable <- as.Date(PCRpos_sx_testing_melt$variable,origin="1970-01-01")
    p <- p + geom_line(aes(x=variable,y=value,group=as.factor(sim)),PCRpos_sx_testing_melt,col="blue",linetype="longdash") + 
            geom_point(aes(x=date,y=PCRpos),data.frame(date=sx_testing_dates,PCRpos=D_S),col="black",size=2,pch=15)
  }

  pdf(paste0("calibration",str,".pdf"),width = 6, height = 5)
  print(p)
  dev.off()
  
  # num_PCRpos_staff <- rep(0,length(sim_pop_list))
  # for (i in 1:length(sim_pop_list)){
  #   num_PCRpos_staff[i] <- sum(sim_pop_list[[i]]$HxPCR[sim_pop_list[[i]]$Resident==0],na.rm = T) 
  # }
  # pdf(paste0("num_PCRpos_staff",str,".pdf"),width = 5, height = 4)
  # hist(num_PCRpos_staff,breaks=seq(min(num_PCRpos_staff)-0.5,max(num_PCRpos_staff)+0.5),xlab="No. PCR positive staff",main="Distribution of no. of PCR positive staff")
  # dev.off()
}
