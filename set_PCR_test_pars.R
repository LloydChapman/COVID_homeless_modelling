# Set PCR test parameters

# Sensitivity as a function of days since start of infectiousness
FNR <- read.csv("digitised_sens_graph.csv",header = F, stringsAsFactors = F)
x <- round(FNR[(round(mu_E)+1):nrow(FNR),1],0)-round(mu_E) # [ ] - Think about whether it should be mu_E or 3 here?
y <- 1 - FNR[(round(mu_E)+1):nrow(FNR),2]
fit <- lm(y ~ bs(x,knots = c(3,4,6,7,10,15,16)-round(mu_E)))
fit_extrap <- approxExtrap(x,y,(x[length(x)]+1):max_days_PCR_pos)

# # Plot PCR sensitivity curve
# pdf("PCR_sens_vs_time_since_presx_infctsnss.pdf",width = 5,height = 4)
# xx <- seq(0,max_days_PCR_pos)
# plot(x,y,pch=19,xlim = c(0,max_days_PCR_pos), ylim = c(0,max(y)), xlab = "Days since start of presymptomatic infectiousness", ylab = "PCR sensitivity")
# lines(xx, sens(xx,fit,fit_extrap,max_days_PCR_pos), col="red")
# dev.off()

# PCR test specificity
spec <- c(1,1,NA,NA,NA,NA,NA) # specificities for states 1 to 7