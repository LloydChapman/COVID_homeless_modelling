# Boston data from Baggett JAMA 2020 and Mosites MMWR 2020

# Load CCMS data if available
CCMS_data <- NULL

# SF CCMS data for proportions in different risk groups
SF_CCMS_data <- read.csv("data/CCMS_data.csv",stringsAsFactors = F)
names(SF_CCMS_data)[1] <- "Date"
SF_CCMS_data$Date <- as.Date(paste0(SF_CCMS_data$Date,"-2020"),format = "%d-%b-%Y")
# Remove empty rows from after Apr 10
SF_CCMS_data <- SF_CCMS_data[SF_CCMS_data$Date<=as.Date("04/10/2020",format="%m/%d/%Y"),]
# Reformat names
names(SF_CCMS_data) <- gsub("\\.\\.\\.|\\.\\.|\\.","_",names(SF_CCMS_data))
# Replace NA's with 0's
SF_CCMS_data[is.na(SF_CCMS_data)] <- 0

# Start and end date of period for which data is available
start_date <- as.Date("4/2/2020",format = "%m/%d/%Y")
end_date <- as.Date("4/3/2020",format = "%m/%d/%Y")

# Date first case identified
date_first_case <- as.Date("3/26/2020",format = "%m/%d/%Y")

# Load Boston case data obtained from https://dashboard.cityofboston.gov/t/Guest_Access_Enabled/views/COVID-19/Dashboard1?:showAppBanner=false&:display_count=n&:showVizHome=n&:origin=viz_share_link&:isGuestRedirectFromVizportal=y&:embed=y
Boston_data <- read.delim("data/Boston_Cases_over_time_data.csv",sep = "\t",stringsAsFactors = F, fileEncoding = "UTF-16")
Boston_case_data <- data.frame(Date = as.Date(Boston_data$Day.of.Timestamp,format = "%B %d, %Y"))
Boston_case_data$Cases <- c(1,diff(Boston_data$Cases..Total.Positive,lag = 1))

# Set number of residents and staff in shelter and duration of simulation
N_res <- 408
N_staff <- 50
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
epsilon <- calc_epsilon(Boston_case_data,end_date-trnsmssn_window+reporting_delay,end_date+reporting_delay,692600,underreporting,homeless_RR) # population estimate from US Census Bureau [https://www.census.gov/quickfacts/fact/table/bostoncitymassachusetts,US/PST045219]

# Flag for whether to count hospitalisations and deaths
hospitalisation <- F # false as data not available

# Set PCR test parameters
sens <- sensitivity("constant",max_days_PCR_pos,const_sens = 0.75) # sensitivity as a function of days since start of infectiousness
spec <- c(1,1,NA,NA,NA,NA,NA) # specificities for states 1 to 7

# PCR testing frequency
testing_dates <- as.Date("4/3/2020",format="%m/%d/%Y") # testing dates
N_tested <- 408 + 50 # number of residents and staff tested during mass testing
sx_testing_dates <- as.Date(integer(), origin = "1970-01-01")
N_sx_tested <- integer() # number of symptomatic individuals tested on each symptomatic testing day

# Initialise variables
Number <- 1:N_pop
Resident <- rep(0,N_pop)
Resident[1:N_res] <- 1
res_present0 <- 1:N_res # assume all residents present for duration of simulation
res_absent0 <- setdiff(1:N_res,res_present0)
staff_present0 <- (N_res+1):N_pop
Present <- ((1:N_pop) %in% c(res_present0,staff_present0))
Risk <- rep(1,N_pop) # N.B. assumes all staff are low risk and under-60
N_under60 <- 44 + 119 + floor(245/2)
N_60plus <- ceiling(245/2)
Hi_Risk_60_only <- sample(res_present0,round(N_60plus*SF_CCMS_data$Hi_Risk_60_only[1]/(SF_CCMS_data$Hi_Risk_60_only[1]+SF_CCMS_data$Hi_Risk_Both_Age_Dx[1])))
Risk[Hi_Risk_60_only] <- 2
Hi_Risk_Dx_only <- sample(setdiff(res_present0,Hi_Risk_60_only),round(N_under60*SF_CCMS_data$Hi_Risk_Dx_only[1]/(SF_CCMS_data$Hi_Risk_Dx_only[1]+SF_CCMS_data$Total_Low_Risk[1])))
Risk[Hi_Risk_Dx_only] <- 3
Hi_Risk_Both_Age_Dx <- sample(setdiff(res_present0,c(Hi_Risk_60_only,Hi_Risk_Dx_only)),round(N_60plus*SF_CCMS_data$Hi_Risk_Both_Age_Dx[1]/(SF_CCMS_data$Hi_Risk_60_only[1]+SF_CCMS_data$Hi_Risk_Both_Age_Dx[1])))
Risk[Hi_Risk_Both_Age_Dx] <- 4
Age <- rep(NA,N_pop)
Age[Risk %in% c(1,3)] <- sample(x=seq(18,59), size=sum(Risk %in% c(1,3)), replace=TRUE, prob = c(rep(44/N_under60,17),rep(119/N_under60,15),rep((N_under60-44-119)/N_under60,10)))
Age[Risk %in% c(2,4)] <- sample(x=seq(60,79), size=sum(Risk %in% c(2,4)), replace=TRUE)

# Observed data
D_S <- NULL
D_T <- 147 + 15 # number of positives in residents and staff during mass testing
D_C <- NULL # set to NULL if no data available on symptom onsets

#  Lower and upper boundaries for priors for R0, E0 and T
lm.low <- c(1,1,14)
lm.upp <- c(8,5,30)
