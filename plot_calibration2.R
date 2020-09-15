plot_calibration <- function(fnm1,fnm2,run_nm,lm.low,lm.upp,end_date,D_S,D_T,D_C=NULL,testing_dates,sx_testing_dates,tol,ttl){
  # Load results
  res <- read.csv(fnm1,stringsAsFactors = F)
  N <- dim(res)[1]
  npar <- dim(res)[2]
  load(fnm2)
  
  # Plot posteriors
  pdf(paste0("R0_posterior",run_nm,".pdf"),width = 5, height = 4)
  hist(res[,1],breaks = seq(lm.low[1],lm.upp[1],by = 0.5),freq = F,xlab = "R0",ylab = "Density",main = "")
  dev.off()
  pdf(paste0("E0_posterior",run_nm,".pdf"),width = 5, height = 4)
  hist(res[,2],breaks = (lm.low[2]-0.5):(lm.upp[2]+0.5),freq = F,xlab = "E0",ylab = "Density",main = "")
  dev.off()
  pdf(paste0("T_posterior",run_nm,".pdf"),width = 5, height = 4)
  hist(res[,3],breaks = (lm.low[3]-0.5):(lm.upp[3]+0.5),freq = F,xlab = "T",ylab = "Density",main = "")
  dev.off()
  pdf(paste0("beta_posterior",run_nm,".pdf"),width = 5, height = 4)
  hist(res[,4],freq = F,xlab = "beta",ylab = "Density",main = "")
  dev.off()

  # Plot pairwise parameter correlation
  p <- vector("list",npar^2)
  # binwdths <- c((max(res[,1])-min(res[,1]))/20,1,1)
  binwdths <- c(0.5,1,1)
  p1 <- ggplot(res,aes(x=R0)) + geom_histogram(aes(y=..density..),binwidth = binwdths[1]) + theme_classic()
  p2 <- ggplot(res,aes(x=E0,y=R0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[2],binwdths[1]),show.legend = F) + theme_classic()
  p3 <- ggplot(res,aes(x=T,y=R0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[3],binwdths[1]),show.legend = F) + theme_classic()
  p4 <- ggplot(res,aes(x=R0,y=E0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[1],binwdths[2]),show.legend = F) + theme_classic()
  p5 <- ggplot(res,aes(x=E0)) + geom_histogram(aes(y=..density..),binwidth = binwdths[2]) + theme_classic()
  p6 <- ggplot(res,aes(x=T,y=E0)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[3],binwdths[2]),show.legend = F) + theme_classic()
  p7 <- ggplot(res,aes(x=R0,y=T)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[1],binwdths[3]),show.legend = F) + theme_classic()
  p8 <- ggplot(res,aes(x=E0,y=T)) + stat_bin2d(aes(fill=..density..),binwidth=c(binwdths[2],binwdths[3]),show.legend = F) + theme_classic()
  p9 <- ggplot(res,aes(x=T)) + geom_histogram(aes(y=..density..),binwidth = binwdths[3]) + theme_classic()
  pdf(paste0("parameter_correlation",run_nm,".pdf"),width = 6,height = 6)
  print(grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,nrow=3,ncol=3))
  dev.off()

  # Daily numbers of new infections
  infections_mat <- matrix(0,nrow = N,ncol = lm.upp[3])
  for (i in 1:N){
    infections_mat[i,lm.upp[3]+((-res[i,3]+1):0)] <- infections_list[[i]]
  }
  colnames(infections_mat) <- as.character(seq.Date(end_date-lm.upp[3]+1,end_date,by=1))
  infections_bnds <- rbind(apply(infections_mat,2,mean),apply(infections_mat,2,range))
  infections_bnds_melt <- melt(infections_bnds)
  infections_bnds_wide <- dcast(infections_bnds_melt,Var2 ~ Var1)
  names(infections_bnds_wide)[names(infections_bnds_wide)=="1"] <- "mean"
  names(infections_bnds_wide)[names(infections_bnds_wide)=="2"] <- "LB"
  names(infections_bnds_wide)[names(infections_bnds_wide)=="3"] <- "UB"
  infections_bnds_wide$Var2 <- as.Date(infections_bnds_wide$Var2)
  pdf(paste0("new_infections",run_nm,".pdf"),width = 4.5,height = 3.75)
  print(ggplot(infections_bnds_wide) + geom_line(aes(x=Var2,y=mean),size=1) + geom_ribbon(aes(x=Var2,ymin=LB,ymax=UB),alpha=0.2) + xlab("Date") + ylab("No. new infections") + theme(legend.position = "none") + theme_classic())
  dev.off()
  
  # Daily numbers of new symptomatic cases
  cases_mat1 <- matrix(0,nrow = N,ncol = lm.upp[3])
  for (i in 1:N){
    cases_mat1[i,lm.upp[3]+((-res[i,3]+1):0)] <- cases_list[[i]]
  }
  colnames(cases_mat1) <- as.character(seq.Date(end_date-lm.upp[3]+1,end_date,by=1))
  cases_melt <- melt(cases_mat1)
  cases_melt$Var2 <- as.Date(cases_melt$Var2)
  pdf(paste0("new_cases",run_nm,".pdf"),width = 5,height = 4)
  print(ggplot() + geom_line(aes(x=Var2,y=value,group=Var1,col=as.factor(Var1)),cases_melt) + xlab("Date") + ylab("No. new symptomatic cases") + theme(legend.position = "none")) #,group=Var2,col=as.factor(Var2)
  dev.off()

  # Numbers of PCR positives from testing
  PCRpos_df <- data.frame(PCRpos_mat)
  names(PCRpos_df) <- testing_dates
  PCRpos_df$sim <- row.names(PCRpos_df)
  PCRpos_melt <- melt(PCRpos_df,id.vars="sim")
  # If only one testing date, subtract a little from the date so simulation output is visible
  if (ncol(PCRpos_mat)==1){
    PCRpos_melt$variable <- as.Date(PCRpos_melt$variable)-0.05  
  } else {
    PCRpos_melt$variable <- as.Date(PCRpos_melt$variable)
  }

  # Calibration plot with numbers of PCR positives
  p <- ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim),col="#bfbfbf"),PCRpos_melt) +
          geom_point(aes(x=variable,y=value,col="#bfbfbf"),PCRpos_melt) +
          geom_point(aes(x=date,y=PCRpos,col="black"),data.frame(date=testing_dates,PCRpos=D_T),size=3) +
          geom_line(aes(x=date,y=CI,group=as.factor(date)),data.frame(date=rep(testing_dates,2),CI=pmax(c(D_T-tol,D_T+tol),0)),size=1) +
          xlab("Date") + xlim(min(PCRpos_melt$variable)-1,max(PCRpos_melt$variable)+1) + ylab("No. PCR positive") + ggtitle(ttl) + theme_classic()
  if (length(sx_testing_dates)!=0){
    PCRpos_sx_testing_df <- data.frame(PCRpos_sx_testing_mat)
    names(PCRpos_sx_testing_df) <- sx_testing_dates
    PCRpos_sx_testing_df$sim <- row.names(PCRpos_sx_testing_df)
    PCRpos_sx_testing_melt <- melt(PCRpos_sx_testing_df,id.vars="sim")
    PCRpos_sx_testing_melt$variable <- as.Date(PCRpos_sx_testing_melt$variable,origin="1970-01-01")
    p <- p + geom_line(aes(x=variable,y=value,group=as.factor(sim),col="#9090ff"),PCRpos_sx_testing_melt,linetype="longdash") +
      geom_point(aes(x=variable,y=value,col="#9090ff"),PCRpos_sx_testing_melt) +
      geom_point(aes(x=date,y=PCRpos,col="#0000FF"),data.frame(date=sx_testing_dates,PCRpos=D_S),size=3) +
      # scale_color_identity(name="Testing",guide="legend",labels=c("early symptomatic cases (simulated)","random (simulated)","random (observed)","early symptomatic cases (observed)"))
      scale_color_manual(name="Type of testing",labels=c("#0000FF"="early symptomatic cases (observed)","#9090ff"="early symptomatic cases (simulated)","#bfbfbf"="random (simulated)","black"="random (observed)"),values=c("#0000FF","#9090ff","#bfbfbf","black"))
    pdf(paste0("calibration",run_nm,".pdf"),width = 7.3, height = 3.75)
  } else {
    p <- p + scale_color_identity()
    pdf(paste0("calibration",run_nm,".pdf"),width = 4.5, height = 3.75)
  }
  print(p)
  dev.off()

  # Calibration plot for numbers of symptomatic cases 
  if (!is.null(D_C)){
    cases_df <- data.frame(cases_mat)
    names(cases_df) <- case_dates
    cases_df$sim <- row.names(cases_df)
    cases_melt <- melt(cases_df,id.vars="sim")
    cases_melt$variable <- as.Date(cases_melt$variable)
    pdf(paste0("cases_calibration",run_nm,".pdf"),width = 4.5, height = 3.75)
    print(ggplot() + geom_line(aes(x=variable,y=value,group=as.factor(sim)),cases_melt,col="#bfbfbf") +
      geom_point(aes(x=variable,y=value,group=as.factor(sim)),cases_melt,col="#bfbfbf",pch=19) +
      geom_point(aes(x=date,y=cases),data.frame(date=case_dates,cases=D_C),col="black",size=3) +
      xlab("Date") + ylab("No. symptom onsets") + theme(legend.position = "none") + theme_classic())
    dev.off()
  }

}
