plot_interventions <- function(fnm,run_nm,ttls){
  load(fnm)
  
  infections_melt <- melt(infections)
  # for (i in 1:6){
  #   print(ggplot(infections_melt[infections_melt$Var3==i,],aes(x=Var2,y=value,group=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. new infections") + theme(legend.position = "none"))
  # }
  pdf(paste0("new_infections_interventions",run_nm,".pdf"),width = 4.5,height = 3.75)
  print(ggplot(infections_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3))) + geom_line() + xlab("Day") + ylab("No. new infections") + scale_color_discrete(name = "Strategy", labels = c("None",as.character(1:(dim(infections)[3]-1))))) # + scale_alpha(guide = "none") #,alpha=1/Var3
  dev.off()
  
  total_infections_melt <- melt(total_infections)
  for (i in 2:ncol(total_infections)){
    pdf(paste0("outbreak_size_distn_strategy_",i-1,run_nm,".pdf"),width = 4.5,height = 3.75)
    print(ggplot() + geom_histogram(aes(x=value,y=..density..,fill=as.factor(Var2)),data = total_infections_melt[total_infections_melt$Var2==1,],binwidth = 10,alpha=0.8)
          + geom_histogram(aes(x=value,y=..density..,fill=as.factor(Var2)),data = total_infections_melt[total_infections_melt$Var2==i,],binwidth = 10,alpha=0.8)
          + xlab("Total no. infections") + ylab("Density")
          + scale_fill_discrete(name = "Strategy", labels = c("None",as.character(i-1)))
          + ggtitle(paste0(i-1,") ",ttls[i]))
          + theme_classic())
    dev.off()
  }
  
  cases_melt <- melt(cases)
  pdf(paste0("new_cases_interventions",run_nm,".pdf"),width = 4.5,height = 3.75)
  print(ggplot(cases_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3))) + geom_line() + xlab("Day") + ylab("No. new cases") + scale_color_discrete(name = "Strategy", labels = c("None",as.character(1:(dim(cases)[3]-1))))) # + scale_alpha(guide = "none") #,alpha=1/Var3
  dev.off()
}

plot_interventions_by_strategy <- function(nsims,T_sim,run_nms,epsilons,interventions,R0s,R0lbls,ttls){
  infections_arr <- array(dim = c(nsims,T_sim,length(run_nms)))
  clrs <- gg_color_hue(length(run_nms))
  lbls <- paste(round(R0s,1),R0lbls)
  for (i in 1:length(epsilons)){
    for (j in 1:length(interventions)){
      for (k in 1:length(run_nms)){
        load(paste0("intvntn_sim_output",run_nms[k],"_epsilon",i-1,".RData"))
        infections_arr[,,k] <- infections[,,j]  
      }
      infections_bnds <- abind(apply(infections_arr,c(2,3),mean),apply(infections_arr,c(2,3),range),along = 1)
      infections_bnds_melt <- melt(infections_bnds)
      infections_bnds_wide <- dcast(infections_bnds_melt,Var2 + Var3 ~ Var1)
      names(infections_bnds_wide)[names(infections_bnds_wide)=="1"] <- "mean"
      names(infections_bnds_wide)[names(infections_bnds_wide)=="2"] <- "LB"
      names(infections_bnds_wide)[names(infections_bnds_wide)=="3"] <- "UB"
      infections_bnds_wide$iR0 <- factor(infections_bnds_wide$Var3,levels=unique(infections_bnds_wide$Var3)[order(unique(infections_bnds_wide$Var3),decreasing = T)])
      pdf(paste0("new_infections_interventions",j-1,"_epsilons",i-1,".pdf"),width = 4.5,height = 3.75)
      # print(ggplot(infections_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3),alpha=1/Var3)) + geom_line() + xlab("Day") + ylab("No. new infections") + scale_color_discrete(name = "R0",labels = paste(round(R0s,1),c("","(Seattle)","(Boston)","(SF)"))) + theme_classic()) # + scale_alpha(guide = "none") #,alpha=1/Var3
      print(ggplot(infections_bnds_wide,aes(x=Var2)) + geom_line(aes(x=Var2,y=mean,color=iR0),size=1) + geom_ribbon(aes(ymin=LB,ymax=UB,fill=iR0),alpha=0.2) + xlab("Day") + ylab("No. new infections") + ylim(0,45) + scale_color_manual(name = "R0",labels = lbls[length(lbls):1],values = clrs[length(clrs):1]) + scale_fill_manual(guide = F,values = clrs[length(clrs):1]) + ggtitle(ttls[j]) + theme_classic()) # + scale_alpha(guide = "none") #,alpha=1/Var3
      dev.off()
    }    
  }
}

plot_epsilon_SA <- function(fnm,run_nms,homeless_RR,R0lbls){
  load(fnm)
  
  intvntn_nms <- c("Symptom screening","Routine PCR testing","Universal mask wearing","Relocation of high-risk individuals","Routine PCR testing of staff only","Combination strategy")
  colnames(prob_outbreak_averted) <- as.character(1:dim(prob_outbreak_averted)[2])
  prob_outbreak_averted_melt <- melt(prob_outbreak_averted)
  
  prob_outbreak_averted_melt$R0 <- R0s[prob_outbreak_averted_melt$Var3]
  prob_outbreak_averted_melt$epsilon <- epsilons1[prob_outbreak_averted_melt$Var2]
  prob_outbreak_averted_melt$mean_bckgrnd_inc <- prob_outbreak_averted_melt$epsilon*1e6/homeless_RR
  
  for (i in 1:length(R0s)){
    pdf(paste0("prob_outbreak_averted_vs_mean_bckgrnd_inc",run_nms[i],".pdf"),width = 6,height = 3.75)
    print(ggplot(prob_outbreak_averted_melt[prob_outbreak_averted_melt$Var3==i,],aes(x=mean_bckgrnd_inc,y=value,group=as.factor(Var1),color=as.factor(Var1))) + geom_point() + xlab("Community incidence (infections/1,000,000/day)") + ylab("Probability of averting outbreak") + ylim(0,1) + scale_color_discrete(name = "Strategy",labels = intvntn_nms) + ggtitle(paste("R0 =",round(R0s[i],1),R0lbls[i])) + theme_classic())
    dev.off()
  }
  
} 

plot_PCR_testing_freq_SA <- function(fnm,run_nm,R0lbls){
  load(fnm)
  
  colnames(prob_outbreak_averted1) <- as.character(1:ncol(prob_outbreak_averted1))
  
  prob_outbreak_averted_melt <- melt(prob_outbreak_averted1)
  
  prob_outbreak_averted_melt$R0 <- R0s[prob_outbreak_averted_melt$Var2]
  # prob_outbreak_averted_melt$testing_freq <- testing_freqs[prob_outbreak_averted_melt$Var1]
  prob_outbreak_averted_melt$days_btw_tests <- days_btw_tests[prob_outbreak_averted_melt$Var1]
  
  pdf(paste0("prob_outbreak_averted_vs_PCR_testing_freq",run_nm,".pdf"), width = 6, height = 4)
  # print(ggplot(prob_outbreak_averted_melt,aes(x=testing_freq,y=value,group=as.factor(R0),color=as.factor(round(R0,1)))) + geom_point() + xlab("PCR testing frequency (testing events/week)") + ylab("Probability of averting outbreak") + scale_color_discrete(name = "R0"))
  print(ggplot(prob_outbreak_averted_melt,aes(x=days_btw_tests,y=value,group=as.factor(R0),color=as.factor(round(R0,1)))) + geom_point() + xlab("Frequency of testing (number of days between tests)") + ylab("Probability of averting outbreak") + ylim(0,1) + scale_color_discrete(name = "R0",labels = paste(round(R0s,1),R0lbls)) + theme_classic())
  dev.off()
  
}