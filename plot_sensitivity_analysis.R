plot_sensitivity_analysis <- function(fnm1,fnm2,base_pars,lbls,ttl,run_nm,dir){
  load(fnm1)
  prob_outbreak_averted_base <- t(read.csv(fnm2,stringsAsFactors = F,row.names = 1))

  perc_diff_pars <- t(apply(pars,1,function(x){100*(x - base_pars)/base_pars}))
  perc_diff <- t(apply(prob_outbreak_averted,1,function(x){100*(x - prob_outbreak_averted_base)/prob_outbreak_averted_base}))
  
  ord <- order(order(lbls))
  lbls <- lbls[order(lbls)] # order label strings alphabetically so colors match labels in plots
  clrs <- gg_color_hue(ncol(pars))
  clrs <- clrs[order(clrs)] # order color strings alphabetically so colors match labels in plots
  
  for (i in 1:ncol(prob_outbreak_averted)){
    df <- data.frame(cbind(perc_diff_pars,prob=prob_outbreak_averted[,i]))
    pdf(paste0(dir,"spider_diag_interventions",i,run_nm,"_SA.pdf"),width = 6,height = 4)
    print(ggplot(df) + geom_line(aes(x=h,y=prob,color=clrs[ord[1]]),data = aggregate(prob ~ h,df,mean),size=1) #+ geom_line(aes(x=h,y=prob,group=as.factor(h)),color=clrs[2],alpha=0.2) #+ scale_color_manual(labels = "h",values=clrs[1])
              + geom_line(aes(x=alpha,y=prob,color=clrs[ord[2]]),data = aggregate(prob ~ alpha,df,mean),size=1) #+ geom_line(aes(x=alpha,y=prob,group=as.factor(alpha)),color=clrs[1],alpha=0.2) 
              + geom_line(aes(x=sens_sx,y=prob,color=clrs[ord[3]]),data = aggregate(prob ~ sens_sx,df,mean),size=1) #+ geom_line(aes(x=sens_sx,y=prob,group=as.factor(sens_sx)),color=clrs[7],alpha=0.2)  
              + geom_line(aes(x=spec_sx,y=prob,color=clrs[ord[4]]),data = aggregate(prob ~ spec_sx,df,mean),size=1) #+ geom_line(aes(x=spec_sx,y=prob,group=as.factor(spec_sx)),color=clrs[9],alpha=0.2) 
              + geom_line(aes(x=sx_pos_PCR_test_compliance,y=prob,color=clrs[ord[5]]),data = aggregate(prob ~ sx_pos_PCR_test_compliance,df,mean),size=1) #+ geom_line(aes(x=sx_pos_PCR_test_compliance,y=prob,group=as.factor(sx_pos_PCR_test_compliance)),color=clrs[10],alpha=0.2)
              + geom_line(aes(x=sens,y=prob,color=clrs[ord[6]]),data = aggregate(prob ~ sens,df,mean),size=1) #+ geom_line(aes(x=sens,y=prob,group=as.factor(sens)),color=clrs[6],alpha=0.2) 
              + geom_line(aes(x=spec,y=prob,color=clrs[ord[7]]),data = aggregate(prob ~ spec,df,mean),size=1) #+ geom_line(aes(x=spec,y=prob,group=as.factor(spec)),color=clrs[8],alpha=0.2) 
              + geom_line(aes(x=routine_PCR_test_compliance,y=prob,color=clrs[ord[8]]),data = aggregate(prob ~ routine_PCR_test_compliance,df,mean),size=1) #+ geom_line(aes(x=routine_PCR_test_compliance,y=prob,group=as.factor(routine_PCR_test_compliance)),color=clrs[5],alpha=0.2) 
              + geom_line(aes(x=mask_eff,y=prob,color=clrs[ord[9]]),data = aggregate(prob ~ mask_eff,df,mean),size=1) #+ geom_line(aes(x=mask_eff,y=prob,group=as.factor(mask_eff)),color=clrs[4],alpha=0.2) 
              + geom_line(aes(x=mask_compliance,y=prob,color=clrs[ord[10]]),data = aggregate(prob ~ mask_compliance,df,mean),size=1) #+ geom_line(aes(x=mask_compliance,y=prob,group=as.factor(mask_compliance)),color=clrs[3],alpha=0.2) 
              + xlab("Percentage difference in parameter") + ylab("Probability of averting outbreak")
              + scale_color_manual(name = "Parameter",labels = lbls,values=clrs)
              + ggtitle(ttl)
              + theme_classic())
    dev.off()
  }
}
