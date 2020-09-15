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
source("ABC_SMC.R")
source("ABC_SMC_functions.R")
source("set_data_and_init_condns.R")
source("COVID_homeless_calibration.R")
source("COVID_homeless_functions.R")
source("plot_calibration.R")
source("process_calibration.R")

# Number of particles
N <- 1000

# Run ABC for SF
ABC_SMC("SF",N,"_SF_22")

# Run ABC for Seattle
ABC_SMC("Seattle_A",N,"_Seattle_A_12")
ABC_SMC("Seattle_B",N,"_Seattle_B_12")
ABC_SMC("Seattle_C",N,"_Seattle_C_12")

# Run ABC for Boston
ABC_SMC("Boston",N,"_Boston_11")
