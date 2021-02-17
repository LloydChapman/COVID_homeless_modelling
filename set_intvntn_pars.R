# Set PCR test parameters
sens <- sensitivity("constant",max_days_PCR_pos,const_sens = 0.75) # sensitivity as a function of days since start of infectiousness
spec <- c(1,1,NA,NA,NA,NA,NA) # specificities for states 1 to 7

# PCR testing frequency
testing_freq <- 2 # testing events per week
testing_days <- floor(seq(1,T_sim,by = 7/testing_freq)) # testing twice per week on 1st and 4th day

# Set intervention parameters
max_PCR_tests_per_week <- 2 # maximum number of PCR tests per week
min_days_btw_tests <- 3 # minimum number of days between PCR tests

# PCR testing once upon entry
entry_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing on entry

# Routine PCR testing
routine_PCR_test_compliance <- 0.8 # 80% compliance with routine PCR testing
sx_pos_PCR_test_compliance <- 0.8 # 80% compliance with PCR testing among those who screen symptom positive

# Mask wearing
mask_compliance <- 0.6 # 60% compliance with universal masking
mask_eff_susc <- 0.4 # 40% reduction in infection rate for susceptible individuals wearing masks
mask_eff_inf <- 0.3 # 30% reduction in transmission from infected individuals wearing masks

# Symptom screening sensitivity and specificity
sens_sx <- c(NA,NA,NA,NA,NA,0.4,NA) # sensitivities for states 1 to 7
spec_sx <- c(0.9,0.9,0.9,0.9,0.9,NA,0.9) # specificities for states 1 to 7