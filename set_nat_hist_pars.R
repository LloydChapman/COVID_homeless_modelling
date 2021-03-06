# Set natural history parameters
mu_E <- 3 # mean of negative-binomially-distributed latent period
r_E <- 4 # shape parameter of negative-binomially-distributed latent period
p_E <- r_E/(r_E+mu_E-1) # 'success' probability of negative-binomially-distributed latent period

# Vector of probabilities of developing severe symptoms by age group and co-morbidity status (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
p_s <- c((0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30,(0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30) # Values from Davies medRxiv 2020 - no co-morbidity dependence as no information on this in paper

h <- 1 # relative infectiousness of asymptomatic individuals (individuals with mild symptoms) from Li Science 2020

mu_p <- 2.3 # mean of negative binomial presymptomatic period
r_p <- 4 # shape parameter of negative-binomially-distributed presymptomatic period
p_p <- r_p/(r_p+mu_p-1) # 'success' probability of negative-binomially-distributed presymptomatic period

mu_sx <- 8 # mean of negative-binomially-distributed duration of symptoms
r_sx <- 4 # shape parameter of negative-binomially-distributed duration of symptoms
p_sx <- r_sx/(r_sx+mu_sx-1) # 'success' probability of negative-binomially-distributed duration of symptoms

alpha <- 2 # relative infectiousness of presymptomatic infection

# Vector of probabilities of hospitalisation by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
p_h <- c((0.021+0.025+0.035+0.077)/4,(0.159+0.262+0.446)/3,(0.044+0.054+0.075+0.165)/4,(0.340+0.561+0.954)/3) # Values from Tuite CMAJ 2020

p_ICU <- 0.261 # probability of ICU admission from Wang JAMA 2020
# Vector of probabilities of death for hospitalised patients admitted to ICU by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
RR_death_comorb <- 16/21/(16/31) # correction to RR of death from Tuite CMAJ 2020
p_d <- c((20*0.17+10*0.23+10*0.31)/40,(0.41+0.55+0.60)/3,RR_death_comorb*(20*0.17+10*0.23+10*0.31)/40,RR_death_comorb*(0.41+0.55+0.60)/3) # Tuite CMAJ 2020

# Construct discrete distribution for duration of PCR positivity
mean_days_PCR_pos <- 20 # mean duration of PCR positivity following symptom onset based on viral shedding dynamics papers
min_days_PCR_pos <- 5 # He Nat Med 2020, Xiao Jrnl Clin Vir 2020
max_days_PCR_pos <- 37 # Zhou Lancet 2020
discrnorm <- discretize(pnorm(x,mean_days_PCR_pos,5), from = min_days_PCR_pos-0.5, to = max_days_PCR_pos+0.5, step = 1, method = "rounding")
discrnorm <- discrnorm/sum(discrnorm) # normalise
# pdf("detectable_viral_load_duration_distn.pdf",width = 5.5,height = 4)
# plot(min_days_PCR_pos:max_days_PCR_pos,discrnorm,xlab = "Duration of detectable viral load from start of infectiousness (days)", ylab = "Probability",pch=19) #, main = "Distribution of duration of detectable viral load")
# dev.off()