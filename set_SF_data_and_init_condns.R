# Load CCMS data from MSC South outbreak
CCMS_data <- read.csv("data/CCMS_data.csv",stringsAsFactors = F)
names(CCMS_data)[1] <- "Date"
CCMS_data$Date <- as.Date(paste0(CCMS_data$Date,"-2020"),format = "%d-%b-%Y")
# Remove empty rows from after Apr 10
CCMS_data <- CCMS_data[CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(CCMS_data))
# Replace NA's with 0's
CCMS_data[is.na(CCMS_data)] <- 0

# start_date <- min(CCMS_data$Date)
start_date <- as.Date("3/28/2020",format = "%m/%d/%Y")
end_date <- max(CCMS_data$Date)

# Date first case identified
date_first_case <- as.Date("4/5/2020",format = "%m/%d/%Y")

# Load SF data from SF DPH [https://data.sfgov.org/COVID-19/COVID-19-Cases-Summarized-by-Date-Transmission-and/tvq9-ec9w]
SF_data <- read.csv("data/COVID-19_Cases_Summarized_by_Date__Transmission_and_Case_Disposition.csv",stringsAsFactors = F)
# Remove deaths
SF_data <- SF_data[!(SF_data$Case.Disposition=="Death"),]
names(SF_data)[names(SF_data)=="Specimen.Collection.Date"] <- "Date"
names(SF_data)[names(SF_data)=="Case.Count"] <- "Cases"
SF_data$Date <- as.Date(SF_data$Date) 
SF_case_data <- aggregate(Cases ~ Date,SF_data,sum)

# Load aggregate PCR testing data and symptomatic case data from SFDPH
agg_PCR_data <- read.csv("data/SF_shelter_PCR_data.csv",colClasses = c("Date","integer","integer"))
agg_PCR_sx_testing_data <- read.csv("data/SF_shelter_PCR_data_sx.csv",colClasses = c("Date","integer","integer"))
case_data <- read.csv("data/SF_shelter_case_data.csv",colClasses = c("Date","integer"))

# Set number of residents and staff in shelter and duration of simulation
N_res <- 350 # Set such that ~255 unique individuals are present in shelter at some point between 3/29 and 4/10  # CCMS_data$Total_Census[1]
N_staff <- 65
N_pop <- N_res + N_staff

# Set weights for presence of residents and staff in shelter
w <- rep(1,N_pop)

# Set natural history parameters
source("set_nat_hist_pars.R")

# Set background transmission rate
reporting_delay <- 7 # days from start of infectiousness = 2 days presymptomatic infectious + 5 days from symptom onset to reporting
trnsmssn_window <- 21 # days
underreporting <- 10 # under-reporting ratio for confirmed cases vs infections
homeless_RR <- 2 # relative-risk of infection for homeless individuals
epsilon <- calc_epsilon(SF_case_data,end_date-trnsmssn_window+reporting_delay,end_date+reporting_delay,881549,underreporting,homeless_RR) # population estimate from US Census Bureau [https://www.census.gov/quickfacts/sanfranciscocountycalifornia]

# Flag for whether to count hospitalisations and deaths
hospitalisation <- F # false as data not available for MSC South outbreak

# Set PCR test parameters
sens <- sensitivity("constant",max_days_PCR_pos,const_sens = 0.75) # sensitivity as a function of days since start of infectiousness
spec <- c(1,1,NA,NA,NA,NA,NA) # specificities for states 1 to 7

# PCR testing frequency
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
Age[Risk %in% c(1,3)] <- sample(x=seq(20,59), size=sum(Risk %in% c(1,3)), replace=TRUE) # [ ] - UPDATE to true age distribution?
Age[Risk %in% c(2,4)] <- sample(x=seq(60,77), size=sum(Risk %in% c(2,4)), replace=TRUE) # [ ] - UPDATE to true age distribution?

# Observed data
D_S <- agg_PCR_sx_testing_data$PositiveTests
D_T <- agg_PCR_data$PositiveTests # number of positives in residents and staff in random testing
D_C <- case_data$Cases # set to NULL if no data available on symptom onsets

#  Lower and upper boundaries for priors for R0, E0 and T
lm.low <- c(1,1,14)
lm.upp <- c(8,5,30)
