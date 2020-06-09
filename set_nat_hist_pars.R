# Set natural history parameters
mu_E <- 3 #5 # mean of negative-binomially-distributed latent period - REDUCED TO 3 DAYS TO ACCOUNT FOR INCLUSION OF PRE-SYMPTOMATIC INFECTIOUS STATE WITH MEAN OF 2.3 DAYS 
r_E <- 4 # shape parameter of negative-binomially-distributed latent period
p_E <- r_E/(r_E+mu_E-1) # 'success' probability of negative-binomially-distributed latent period

# Vector of probabilities of developing severe symptoms by age group and co-morbidity status (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
p_s <- c((0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30,(0.37+0.42+0.51+0.59)/4,(10*0.72+20*0.76)/30) # Values from Davies medRxiv 2020 - no co-morbidity dependence as no information on this in paper

h <- 0.55 # relative infectiousness of asymptomatic individuals (individuals with mild symptoms) from Li Science 2020

mu_p <- 2.3 # mean of negative binomial presymptomatic period
r_p <- 4 # shape parameter of negative-binomially-distributed presymptomatic period
p_p <- r_p/(r_p+mu_p-1) # 'success' probability of negative-binomially-distributed presymptomatic period

mu_sx <- 8 # mean of negative-binomially-distributed duration of symptoms
r_sx <- 4 # shape parameter of negative-binomially-distributed duration of symptoms
p_sx <- r_sx/(r_sx+mu_sx-1) # 'success' probability of negative-binomially-distributed duration of symptoms

alpha <- 0.44/(1-0.44)*mu_sx/mu_p # relative infectiousness of presymptomatic infection (He 2020 Nat Med)

# Vector of probabilities of hospitalisation by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
p_h <- c((0.021+0.025+0.035+0.077)/4,(0.159+0.262+0.446)/3,(0.044+0.054+0.075+0.165)/4,(0.340+0.561+0.954)/3) # Values from Tuite CMAJ 2020

p_ICU <- 0.261 # probability of ICU admission from Wang JAMA 2020
# Vector of probabilities of death for hospitalised patients admitted to ICU by age-group and co-morbidity status for severe cases (rows: 1- under-60 & no co-morbidities, 2- 60+ & no co-morbidities, 3- under-60 & co-morbidities, 4- 60+ & co-morbidities)
p_d <- c((20*0.17+10*0.23+10*0.31)/40,(0.41+0.55+0.60)/3,(20*0.45+10*0.60+10*0.81)/40,(1+1+1)/3)

# Construct discrete distribution for duration of PCR positivity
mean_days_PCR_pos <- 20 # mean duration of PCR positivity following symptom onset based on viral shedding dynamics papers
min_days_PCR_pos <- 5 # He Nat Med 2020, Xiao Jrnl Clin Vir 2020
max_days_PCR_pos <- 37 # Zhou Lancet 2020 #25
discrnorm <- discretize(pnorm(x,mean_days_PCR_pos,5), from = min_days_PCR_pos-0.5, to = max_days_PCR_pos+0.5, step = 1, method = "rounding")
discrnorm <- discrnorm/sum(discrnorm) # normalise

epsilon <- (62/2959)/mean_days_PCR_pos * 2 #0.001 # mean_daily_inc/0.14 # transmission rate outside shelter assuming 1/0.14=7.1x as many infections as confirmed cases from Li Science 2020 - [ ] consider making this time-dependent