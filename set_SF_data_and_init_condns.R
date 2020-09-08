# Load CCMS data from MSC South outbreak
CCMS_data <- read.csv("../Data/CCMS_data.csv",stringsAsFactors = F)
names(CCMS_data)[1] <- "Date"
CCMS_data$Date <- as.Date(CCMS_data$Date,format = "%d-%b")
# Remove empty rows from after Apr 10
CCMS_data <- CCMS_data[CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(CCMS_data))
# Replace NA's with 0's
CCMS_data[is.na(CCMS_data)] <- 0

# start_date <- min(CCMS_data$Date)
start_date <- as.Date("3/28/2020",format = "%m/%d/%Y")
end_date <- max(CCMS_data$Date)
mass_testing_start_date <- as.Date("4/8/2020", format = "%m/%d/%Y")
mass_testing_end_date <- as.Date("4/9/2020", format = "%m/%d/%Y")
mass_testing_start_day <- as.integer(mass_testing_start_date-start_date+1)
mass_testing_end_day <- as.integer(mass_testing_end_date-start_date+1)

# Date first case identified
date_first_case <- as.Date("4/5/2020",format = "%m/%d/%Y")

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("../Data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
SF_case_data <- aggregate(Case.Count ~ Date,SF_data,sum)

# Load DPH data
linelist <- read.csv("../Data/LineList2.csv",stringsAsFactors = F)
PCR_data <- read.csv("../Data/LabData2.csv",stringsAsFactors = F)

datecols <- grep("Date",names(linelist))
for (i in 1:length(datecols)){
  linelist[,datecols[i]] <- as.Date(linelist[,datecols[i]],format = "%m/%d/%y")
}
linelist$Age[linelist$Age=="unknown"]<-""
linelist$Age <- as.integer(linelist$Age)

# Aggregate PCR testing data and symptomatic case data
idx <- (linelist$TestingDate>=start_date & linelist$TestingDate<=end_date & (linelist$Symptomatic=="" | linelist$TestingDate>=mass_testing_start_date))
agg_PCR_data <- aggregate(ID ~ TestingDate, linelist[idx,], length)
names(agg_PCR_data)[2] <- "Tests"
sum(agg_PCR_data$Tests)
tmp <- aggregate(ID ~ TestingDate, linelist[linelist$Result=="positive" & idx,], length)
names(tmp)[2] <- "PositiveTests"
agg_PCR_data <- merge(agg_PCR_data,tmp,all.x = T)
agg_PCR_data$PositiveTests[is.na(agg_PCR_data$PositiveTests)] <- 0

idx <- (linelist$TestingDate>=start_date & linelist$Symptomatic %in% c("Y","Y_IQ") & linelist$TestingDate<mass_testing_start_date)
agg_PCR_sx_testing_data <- aggregate(ID ~ TestingDate, linelist[idx,], length)
names(agg_PCR_sx_testing_data)[2] <- "Tests"
tmp <- aggregate(ID ~ TestingDate, linelist[linelist$Result=="positive" & idx,], length)
names(tmp)[2] <- "PositiveTests"
agg_PCR_sx_testing_data <- merge(agg_PCR_sx_testing_data,tmp,all.x = T)
agg_PCR_sx_testing_data$PositiveTests[is.na(agg_PCR_sx_testing_data$PositiveTests)] <- 0

# Count symptomatic individuals (exclude those with negative PCR results who reported having had symptoms at I & Q sites)
tmp <- aggregate(ID ~ SymptomOnsetDate, linelist[!(linelist$Symptomatic=="Y_IQ" & linelist$Result=="negative"),], length)
names(tmp)[2] <- "Cases"
case_data <- data.frame(date = seq.Date(start_date,end_date,1),Cases=rep(0,as.integer(end_date-start_date+1)))
case_data$Cases[case_data$date %in% tmp$SymptomOnsetDate] <- tmp$Cases[tmp$SymptomOnsetDate %in% case_data$date]

# Set number of residents and staff in shelter and duration of simulation
N_res <- 350 # Set such that ~255 unique individuals are present in shelter at some point between 3/29 and 4/10  # CCMS_data$Total_Census[1]
N_staff <- 65
N_pop <- N_res + N_staff
# T_sim <- as.integer(end_date - start_date + 1) #13 # testing 5 days after 1st cases confirmed on 4/5 + 2 days presymptomatic infectiousness (assumed) + 6 days of symptoms of 1st case  # 30 # Days; Total of 1 month

# Set weights for presence of residents and staff in shelter
w <- rep(1,N_pop) #c(rep(1,N_res),rep(1/2,N_staff))

# Set natural history parameters
source("set_nat_hist_pars.R")

# Set background transmission rate
reporting_delay <- 7 # days from start of infectiousness = 2 days presymptomatic infectious + 5 days from symptom onset to reporting
trnsmssn_window <- 21 # days
underreporting <- 10 # under-reporting ratio for confirmed cases vs infections
homeless_RR <- 2 # relative-risk of infection for homeless individuals
mean_daily_cases <- mean(SF_case_data$Case.Count[SF_case_data$Date>=end_date-trnsmssn_window+reporting_delay & SF_case_data$Date<=end_date+reporting_delay]) # mean of confirmed cases for period of interest
mean_daily_inc <- mean_daily_cases/881549 # population estimate from US Census Bureau [https://www.census.gov/quickfacts/sanfranciscocountycalifornia]
epsilon <- mean_daily_inc*underreporting*homeless_RR # adjusted transmission rate outside shelter

# Flag for whether to count hospitalisations and deaths
hospitalisation <- F # false as data not available for MSC South outbreak

# Set PCR test parameters
source("set_PCR_test_pars.R")

# PCR testing frequency
# testing_dates <- as.Date(c("4/8/2020","4/9/2020"),format="%m/%d/%Y") # testing dates
# testing_days <- as.integer(testing_dates - start_date + 1) # days on which testing occurred
# N_tested <- rep((144-10)/2,2) # number tested during mass testing = total tested from press release - number tested during first 4 days = 144 - 10, # number tested on each testing day #c(72,72) # assume all testing happened during mass testing and was equally split between 4/8 and 4/9
# sx_testing_dates <- as.Date(c("4/4/2020","4/5/2020","4/6/2020","4/7/2020"),format="%m/%d/%Y")
# sx_testing_days <- as.integer(sx_testing_dates - start_date + 1)
testing_dates <- agg_PCR_data$TestingDate # dates of PCR testing of "random" individuals # [ ] - N.B. in practice these are likely to be close contacts/suspected cases
N_tested <- agg_PCR_data$Tests[agg_PCR_data$TestingDate %in% testing_dates]  # number tested on each testing day
sx_testing_dates <- agg_PCR_sx_testing_data$TestingDate # dates of PCR testing of symptomatic individuals
N_sx_tested <- agg_PCR_sx_testing_data$Tests[agg_PCR_sx_testing_data$TestingDate %in% sx_testing_dates] # number of symptomatic individuals tested on each symptomatic testing day
case_dates <- case_data$date # dates with numbers of symptomatic cases recorded

# Initialise variables
Number <- 1:N_pop
Resident <- rep(0,N_pop)
Resident[1:N_res] <- 1
res_present0 <- sample(1:N_res,CCMS_data$Total_Census[1])
res_absent0 <- setdiff(1:N_res,res_present0)
staff_present0 <- (N_res+1):N_pop
Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
v <- round(N_res*CCMS_data[1,2:ncol(CCMS_data)]/CCMS_data[1,2],0)
Hi_Risk_60_only_present0 <- sample(res_present0,CCMS_data$Hi_Risk_60_only[1])
Hi_Risk_60_only_absent0 <- sample(res_absent0,v$Hi_Risk_60_only - CCMS_data$Hi_Risk_60_only[1])
Risk[c(Hi_Risk_60_only_present0,Hi_Risk_60_only_absent0)] <- 2
Hi_Risk_Dx_only_present0 <- sample(setdiff(res_present0,Hi_Risk_60_only_present0),CCMS_data$Hi_Risk_Dx_only[1])
Hi_Risk_Dx_only_absent0 <- sample(setdiff(res_absent0,Hi_Risk_60_only_absent0),v$Hi_Risk_Dx_only - CCMS_data$Hi_Risk_Dx_only)
Risk[c(Hi_Risk_Dx_only_present0,Hi_Risk_Dx_only_absent0)] <- 3
Hi_Risk_Both_Age_Dx_present0 <- sample(setdiff(res_present0,c(Hi_Risk_60_only_present0,Hi_Risk_Dx_only_present0)),CCMS_data$Hi_Risk_Both_Age_Dx[1])
Hi_Risk_Both_Age_Dx_absent0 <- sample(setdiff(res_absent0,c(Hi_Risk_60_only_absent0,Hi_Risk_Dx_only_absent0)),v$Hi_Risk_Both_Age_Dx-CCMS_data$Hi_Risk_Both_Age_Dx[1])
Risk[c(Hi_Risk_Both_Age_Dx_present0,Hi_Risk_Both_Age_Dx_absent0)] <- 4
Age <- rep(NA,N_pop)
Age[Risk %in% c(1,3)] <- sample(x=seq(min(linelist$Age,na.rm = T),59), size=sum(Risk %in% c(1,3)), replace=TRUE) # [ ] - UPDATE to true age distribution?
Age[Risk %in% c(2,4)] <- sample(x=seq(60,max(linelist$Age,na.rm = T)), size=sum(Risk %in% c(2,4)), replace=TRUE) # [ ] - UPDATE to true age distribution?

# Observed data
D_S <- agg_PCR_sx_testing_data$PositiveTests
D_T <- agg_PCR_data$PositiveTests # number of positives in residents and staff in random testing
D_C <- case_data$Cases # set to NULL if no data available on symptom onsets

#  Lower and upper boundaries for priors for R0, E0 and T
lm.low <- c(1,1,14)
lm.upp <- c(8,5,30)
