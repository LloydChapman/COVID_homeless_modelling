rm(list=ls())

library(mvnfast)
library(actuar)
library(splines)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(gsubfn)
library(gridExtra)

setwd("~/Dropbox/Homeless/Code") 
source("ABC_SMC_MNN.R")
source("ABC_SMC_MNN_functions.R")
source("set_data_and_init_condns.R")
source("COVID_homeless_calibration2.R")
source("COVID_homeless_functions.R")
source("plot_calibration2.R")
source("process_calibration.R")

# Number of particles
N <- 1000

# Run ABC for SF
ABC_SMC_MNN("SF",N,"_SF_22")

# Run ABC for Seattle
ABC_SMC_MNN("Seattle_A",N,"_Seattle_A_12")
ABC_SMC_MNN("Seattle_B",N,"_Seattle_B_12")
ABC_SMC_MNN("Seattle_C",N,"_Seattle_C_12")

# Run ABC for Boston
ABC_SMC_MNN("Boston",N,"_Boston_11")
