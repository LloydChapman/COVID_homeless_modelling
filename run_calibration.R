# install.packages(c("tmvtnorm","actuar"))
# library(tmvtnorm)
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

# Run ABC for SF
ABC_SMC_MNN("SF","_15")

# Run ABC for Seattle
ABC_SMC_MNN("Seattle_A","_Seattle_A_5")
ABC_SMC_MNN("Seattle_B","_Seattle_B_5")
ABC_SMC_MNN("Seattle_C","_Seattle_C_5")

# # Run ABC for Boston
ABC_SMC_MNN("Boston","_Boston_4")
