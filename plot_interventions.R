# plot_interventions <- function()

infections_melt <- melt(infections)
for (i in 1:6){
  print(ggplot(infections_melt[infections_melt$Var3==i,],aes(x=Var2,y=value,group=as.factor(Var1))) + geom_line() + xlab("Day") + ylab("No. new infections") + theme(legend.position = "none"))
}

pdf("new_infections_interventions_1.pdf",width = 6,height = 5)
ggplot(infections_melt,aes(x=Var2,y=value,group=interaction(as.factor(Var1),as.factor(Var3)),col=as.factor(Var3))) + geom_line() + xlab("Day") + ylab("No. new infections") + scale_color_discrete(name = "Strategy", labels = c("None","1","2","3","4","5")) # + scale_alpha(guide = "none") #,alpha=1/Var3
dev.off()