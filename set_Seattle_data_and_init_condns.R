# Seattle data from MMWR paper and report

# Load CCMS data if available
CCMS_data <- NULL

# SF CCMS data for proportions in different risk groups
SF_CCMS_data <- read.csv("../Data/CCMS_data.csv",stringsAsFactors = F)
names(SF_CCMS_data)[1] <- "Date"
SF_CCMS_data$Date <- as.Date(SF_CCMS_data$Date,format = "%d-%b")
# Remove empty rows from after Apr 10
SF_CCMS_data <- SF_CCMS_data[SF_CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(SF_CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(SF_CCMS_data))
# Replace NA's with 0's
SF_CCMS_data[is.na(SF_CCMS_data)] <- 0

# Start and end date of period for which data is available
start_date <- as.Date("4/1/2020",format = "%m/%d/%Y")
end_date <- as.Date("4/8/2020",format = "%m/%d/%Y")

# Load Seattle case data obtained from https://www.kingcounty.gov/depts/health/covid-19/data/daily-summary.aspx
Seattle_case_data <- read.csv("../Data/Seattle_cases.csv",stringsAsFactors = F)
Seattle_case_data$Date <- as.Date(Seattle_case_data$Date,format = "%m/%d/%y")

# Set number of residents and staff in shelter and duration of simulation
N_res <- 245
N_staff <- 27
N_pop <- N_res + N_staff

# Set weights for presence of residents and staff in shelter
w <- c(rep(1,N_res),rep(1/2,N_staff))

# Set natural history parameters
source("set_nat_hist_pars.R")

# Set background transmission rate
reporting_delay <- 7 # days from start of infectiousness = 2 days presymptomatic infectious + 5 days from symptom onset to reporting
trnsmssn_window <- 21 # days
underreporting <- 10 # under-reporting ratio for confirmed cases vs infections
homeless_RR <- 2 # relative-risk of infection for homeless individuals
mean_daily_cases <- mean(Seattle_case_data$Cases[Seattle_case_data$Date>=end_date-trnsmssn_window+reporting_delay & Seattle_case_data$Date<=end_date+reporting_delay]) # mean of confirmed cases for period of interest
mean_daily_inc <- mean_daily_cases/753675 # population estimate from US Census Bureau [https://www.census.gov/quickfacts/fact/table/seattlecitywashington/PST045219]
epsilon <- mean_daily_inc*underreporting*homeless_RR # adjusted transmission rate outside shelter

# Flag for whether to count hospitalisations and deaths
hospitalisation <- F # false as data not available

# Set PCR test parameters
source("set_PCR_test_pars.R")

# PCR testing frequency
testing_dates <- as.Date(c("4/1/2020","4/8/2020"),format="%m/%d/%Y") # testing dates
N_tested <- c(154,103) # numbers tested during mass testing
sx_testing_dates <- as.Date(integer(), origin = "1970-01-01")
N_sx_tested <- integer() # number of symptomatic individuals tested on each symptomatic testing day

# Initialise variables
Number <- 1:N_pop
Resident <- rep(0,N_pop)
Resident[1:N_res] <- 1
res_present0 <- 1:N_res # assume all residents present for duration of simulation for now # [ ] - UPDATE this with info from MMWR paper # sample(1:N_res,SF_CCMS_data$Total_Census[1])
res_absent0 <- setdiff(1:N_res,res_present0)
staff_present0 <- (N_res+1):N_pop
Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
v <- round(N_res*SF_CCMS_data[1,2:ncol(SF_CCMS_data)]/SF_CCMS_data[1,2],0)
Hi_Risk_60_only_present0 <- sample(res_present0,round(length(res_present0)/4))
Hi_Risk_60_only_absent0 <- sample(res_absent0,round(length(res_absent0)/4))
Risk[c(Hi_Risk_60_only_present0,Hi_Risk_60_only_absent0)] <- 2
Hi_Risk_Dx_only_present0 <- sample(setdiff(res_present0,Hi_Risk_60_only_present0),round(length(res_present0)/4))
Hi_Risk_Dx_only_absent0 <- sample(setdiff(res_absent0,Hi_Risk_60_only_absent0),round(length(res_absent0)/4))
Risk[c(Hi_Risk_Dx_only_present0,Hi_Risk_Dx_only_absent0)] <- 3
Hi_Risk_Both_Age_Dx_present0 <- sample(setdiff(res_present0,c(Hi_Risk_60_only_present0,Hi_Risk_Dx_only_present0)),length(res_present0)-3*round(length(res_present0)/4))
Hi_Risk_Both_Age_Dx_absent0 <- sample(setdiff(res_absent0,c(Hi_Risk_60_only_absent0,Hi_Risk_Dx_only_absent0)),length(res_absent0)-3*round(length(res_absent0)/4))
Risk[c(Hi_Risk_Both_Age_Dx_present0,Hi_Risk_Both_Age_Dx_absent0)] <- 4
Age <- rep(NA,N_pop)
Age[Risk %in% c(1,3)] <- sample(x=seq(50,59), size=sum(Risk %in% c(1,3)), replace=TRUE) # [ ] - UPDATE to true age distribution?
Age[Risk %in% c(2,4)] <- sample(x=seq(60,77), size=sum(Risk %in% c(2,4)), replace=TRUE) # [ ] - UPDATE to true age distribution?

# Observed data
D_S <- NULL
D_T <- c(15,16) # number of positives during mass testing
D_C <- NULL # set to NULL if no data available on symptom onsets

#  Lower and upper boundaries for priors for R0, E0 and T
lm.low <- c(1,1,14)
lm.upp <- c(8,5,30)
