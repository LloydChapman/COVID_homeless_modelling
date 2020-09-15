rm(list=ls())

# Load CCMS data from MSC South outbreak
CCMS_data <- read.csv("data/CCMS_data.csv",stringsAsFactors = F)
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
idx <- (linelist$TestingDate>=start_date & linelist$TestingDate<=end_date & (linelist$Symptomatic=="" | linelist$TestingDate>=mass_testing_start_date) & !is.na(linelist$TestingDate))
agg_PCR_data <- aggregate(ID ~ TestingDate, linelist[idx,], length)
names(agg_PCR_data)[2] <- "Tests"
sum(agg_PCR_data$Tests)
tmp <- aggregate(ID ~ TestingDate, linelist[linelist$Result=="positive" & idx,], length)
names(tmp)[2] <- "PositiveTests"
agg_PCR_data <- merge(agg_PCR_data,tmp,all.x = T)
agg_PCR_data$PositiveTests[is.na(agg_PCR_data$PositiveTests)] <- 0

# Save PCR data from "random" testing
write.csv(agg_PCR_data,"data/SF_shelter_PCR_data.csv",row.names = F)

idx <- (linelist$TestingDate>=start_date & linelist$Symptomatic %in% c("Y","Y_IQ") & linelist$TestingDate<mass_testing_start_date & !is.na(linelist$TestingDate))
agg_PCR_sx_testing_data <- aggregate(ID ~ TestingDate, linelist[idx,], length)
names(agg_PCR_sx_testing_data)[2] <- "Tests"
tmp <- aggregate(ID ~ TestingDate, linelist[linelist$Result=="positive" & idx,], length)
names(tmp)[2] <- "PositiveTests"
agg_PCR_sx_testing_data <- merge(agg_PCR_sx_testing_data,tmp,all.x = T)
agg_PCR_sx_testing_data$PositiveTests[is.na(agg_PCR_sx_testing_data$PositiveTests)] <- 0

# Save PCR data from testing of early symptomatics
write.csv(agg_PCR_sx_testing_data,"data/SF_shelter_PCR_data_sx.csv",row.names = F)

# Count symptomatic individuals (exclude those with negative PCR results who reported having had symptoms at I & Q sites)
tmp <- aggregate(ID ~ SymptomOnsetDate, linelist[!(linelist$Symptomatic=="Y_IQ" & linelist$Result=="negative"),], length)
names(tmp)[2] <- "Cases"
case_data <- data.frame(date = seq.Date(start_date,end_date,1),Cases=rep(0,as.integer(end_date-start_date+1)))
case_data$Cases[case_data$date %in% tmp$SymptomOnsetDate] <- tmp$Cases[tmp$SymptomOnsetDate %in% case_data$date]

# Save symptomatic case data
write.csv(case_data,"data/SF_shelter_case_data.csv",row.names = F)
