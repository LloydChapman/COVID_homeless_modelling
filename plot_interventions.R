plot_interventions <- function(fnm,run_nm){
  load(fnm)
  
  infections_melt <- melt(infections)
  # for (i in 1:6){
  #   print(ggplot(infections_melt[infections_melt$Var3==i,],aes(x=Var2,y=value,group=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. new infections") + theme(legend.position = "none"))
  # }
  pdf(paste0("new_infections_interventions",run_nm,".pdf"),width = 6,height = 5)
  print(ggplot(infections_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3))) + geom_line() + xlab("Day") + ylab("No. new infections") + scale_color_discrete(name = "Strategy", labels = c("None","1","2","3","4","5"))) # + scale_alpha(guide = "none") #,alpha=1/Var3
  dev.off()
  
  total_infections_melt <- melt(total_infections)
  for (i in 2:ncol(total_infections)){
    pdf(paste0("outbreak_size_distn_strategy_",i-1,run_nm,".pdf"),width = 6,height = 5)
    print(ggplot() + geom_histogram(aes(x=value,y=..density..,fill=as.factor(Var2)),data = total_infections_melt[total_infections_melt$Var2==1,],binwidth = 10,alpha=0.8)
          + geom_histogram(aes(x=value,y=..density..,fill=as.factor(Var2)),data = total_infections_melt[total_infections_melt$Var2==i,],binwidth = 10,alpha=0.8)
          + xlab("Total no. infections") + ylab("Density")
          + scale_fill_discrete(name = "Strategy", labels = c("None",as.character(i-1))))
    dev.off()
  }
  
  cases_melt <- melt(cases)
  pdf(paste0("new_cases_interventions",run_nm,".pdf"),width = 6,height = 5)
  print(ggplot(cases_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3))) + geom_line() + xlab("Day") + ylab("No. new cases") + scale_color_discrete(name = "Strategy", labels = c("None","1","2","3","4","5"))) # + scale_alpha(guide = "none") #,alpha=1/Var3
  dev.off()
}
