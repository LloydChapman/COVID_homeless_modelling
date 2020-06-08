library(actuar)
library(splines)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(gsubfn)

setwd("~/Dropbox/Homeless/Code")
source("COVID_homeless_interventions.R")
source("COVID_homeless_functions.R")

# Number of simulations per intervention strategy
nsims <- 100

## Set intervention strategies
# 1 - Active symptom screening
# 2 - PCR testing once upon entry
# 3 - Routine PCR testing
# 4 - Passive symptom screening
# 5 - Universal masking
# 6 - Removal of high-risk individuals
interventions <- list(NULL,1,c(1,2),c(3,4),c(1,5),c(1,6),c(1,2,3,5,6))

# Matrix to store total number of infections from each simulation
new_infections <- matrix(NA,length(interventions),nsims)

# Load CCMS data from MSC South outbreak
CCMS_data <- read.csv("../Data/CCMS_data.csv",stringsAsFactors = F)
names(CCMS_data)[1] <- "Date"
CCMS_data$Date <- as.Date(CCMS_data$Date,format = "%d-%b")
# Remove empty rows from after Apr 10
CCMS_data <- CCMS_data[CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(CCMS_data))

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
SF_case_data <- aggregate(Case.Count ~ Date,SF_data,sum)

# Load posterior distribution for beta from calibration
beta_pstr <- read.csv("results_ABC_SMC_MNN_gen_10_9.csv",stringsAsFactors = F)
beta <- median(beta_pstr$V1)

# Run simulations for each intervention strategy
for (i in 1:length(interventions)){
  for (j in 1:nsims){
    res <- COVID_homeless_intervention_model(beta,CCMS_data,SF_case_data,interventions[[i]])
    new_infections[i,j] <- sum(res$new_infections)
  }  
}

print(rowMeans(new_infections))
# Save total number of new infections per simulation
write.table(cbind(new_infections,rowMeans(new_infections)),"new_onsets_interventions.csv",sep = ",",col.names = F, row.names = F)

# Quality checks
# Daily PCR testing with no max number of weekly tests # [ ] - changes hard coded in interventions code, need to UPDATE
new_infections_PCR <- numeric(nsims)
for (j in 1:nsims){
  res3 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,c(3,4))
  new_infections_PCR[j] <- sum(res3$new_infections)
}

barplot(res3$new_infections)
View(res3$sim_pop)
View(res3$state)
View(res3$presence)
tE3 <- getEventTime1(res3$state,2)
tI_m_p3 <- getEventTime1(res3$state,3)
tI_s_p3 <- getEventTime1(res3$state,4)
View(cbind(tE3,tI_m_p3,tI_s_p3,res3$sim_pop$DayRemoved))
mean(new_infections_PCR)

# write.table(new_infections_PCR,"new_onsets_daily_PCR_no_max_per_week.csv",sep = ",",col.names = F, row.names = F)

## beta=0.01, epsilon=0, no migration in/out only hospitalisation of some clinical cases 
# No interventions
set.seed(1)
res4 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,NULL)
View(res4$sim_pop)
View(res4$state)
View(res4$presence)
plot(colSums(res4$presence))
barplot(res4$new_infections)
sum(res4$new_infections) # 302
# Passive sx screening and daily PCR testing with 100% compliance and sensitivity=1
set.seed(1)
res5 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,c(3,4))
View(res5$sim_pop)
View(res5$state)
View(res5$presence)
plot(colSums(res5$presence))
barplot(res5$new_infections)
sum(res5$new_infections) # 295

## beta=0.005
# No interventions
set.seed(1)
res6 <- COVID_homeless_intervention_model(0.005,CCMS_data,SF_case_data,NULL)
View(res6$sim_pop)
View(res6$state)
View(res6$presence)
plot(colSums(res6$presence))
barplot(res6$new_infections)
sum(res6$new_infections) # 302
# Passive sx screening and daily PCR testing with 100% compliance and sensitivity=1
set.seed(1)
res7 <- COVID_homeless_intervention_model(0.005,CCMS_data,SF_case_data,c(3,4))
View(res7$sim_pop)
View(res7$state)
View(res7$presence)
plot(colSums(res7$presence))
barplot(res7$new_infections)
sum(res7$new_infections) # 260

## 100% effective masking
set.seed(1)
res8 <- COVID_homeless_intervention_model(0.01,CCMS_data,SF_case_data,5)
View(res8$sim_pop)
View(res8$state)
View(res8$presence)
plot(colSums(res8$presence))
barplot(res8$new_infections)
sum(res8$new_infections) # 2 - only the initial cases


