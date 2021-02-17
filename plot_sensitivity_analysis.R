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
    
    agg_df <- vector("list",ncol(perc_diff_pars))
    for (j in 1:ncol(perc_diff_pars)){
      agg_df[[j]] <- aggregate(df$prob,by=list(df[,j]),FUN=mean)
      names(agg_df[[j]]) <- c("value","prob")
      agg_df[[j]]$par <- names(df)[j]
    }
    write.csv(do.call(rbind,agg_df),paste0(dir,"mean_prob_outbreak_averted_interventions",i,run_nm,"_SA.csv"),row.names = F)
    
    pdf(paste0(dir,"spider_diag_interventions",i,run_nm,"_SA.pdf"),width = 6,height = 4)
    print(ggplot(df) + geom_line(aes(x=value,y=prob,color=clrs[ord[1]]),data = agg_df[[1]],size=1) #+ geom_line(aes(x=h,y=prob,group=as.factor(h)),color=clrs[ord[1]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[2]]),data = agg_df[[2]],size=1) #+ geom_line(aes(x=alpha,y=prob,group=as.factor(alpha)),color=clrs[ord[2]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[3]]),data = agg_df[[3]],size=1) #+ geom_line(aes(x=sens_sx,y=prob,group=as.factor(sens_sx)),color=clrs[ord[3]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[4]]),data = agg_df[[4]],size=1) #+ geom_line(aes(x=spec_sx,y=prob,group=as.factor(spec_sx)),color=clrs[ord[4]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[5]]),data = agg_df[[5]],size=1) #+ geom_line(aes(x=sx_pos_PCR_test_compliance,y=prob,group=as.factor(sx_pos_PCR_test_compliance)),color=clrs[ord[5]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[6]]),data = agg_df[[6]],size=1) #+ geom_line(aes(x=sens,y=prob,group=as.factor(sens)),color=clrs[ord[6]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[7]]),data = agg_df[[7]],size=1) #+ geom_line(aes(x=spec,y=prob,group=as.factor(spec)),color=clrs[ord[7]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[8]]),data = agg_df[[8]],size=1) #+ geom_line(aes(x=routine_PCR_test_compliance,y=prob,group=as.factor(routine_PCR_test_compliance)),color=clrs[ord[8]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[9]]),data = agg_df[[9]],size=1) #+ geom_line(aes(x=mask_eff,y=prob,group=as.factor(mask_eff_susc),color=clrs[ord[9]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[10]]),data = agg_df[[10]],size=1) #+ geom_line(aes(x=mask_eff,y=prob,group=as.factor(mask_eff_inf)),color=clrs[ord[10]],alpha=0.2)
              + geom_line(aes(x=value,y=prob,color=clrs[ord[11]]),data = agg_df[[11]],size=1) #+ geom_line(aes(x=mask_compliance,y=prob,group=as.factor(mask_compliance)),color=clrs[ord[11]],alpha=0.2)
              + xlab("Percentage difference in parameter") + ylab("Probability of averting outbreak")
              + scale_color_manual(name = "Parameter",labels = lbls,values=clrs)
              + ggtitle(ttl)
              + theme_classic())
    dev.off()
  }
}
