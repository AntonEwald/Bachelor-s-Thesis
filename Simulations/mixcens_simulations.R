#--------------------------------------Code for Simulation-----------------------------------------------------
#@author: Anton Holm
# This is a simulation used for my Bachelor's Thesis in Mathematical Statistics at Stockholm University
# Spring, 2020
#--------------------------------------Code for Simulation-----------------------------------------------------

library(lmec)
library(tidyverse)
load("yearly_variances.Rda")

set.seed(19931031)

kovariat <- rep(0:10-mean(0:10), each = 12) #Creates covariates
nr_of_covariates <- length(unique(kovariat)) # Number of unique covariates
N = length(kovariat) #Denotes the number of datapoints per simulation to N
LOQ_fraction2 = c(0.8) #The proportion of data to be censored when intercept and slope equals zero
simulations = length(LOQ_fraction2) # How many simulations with different input (in this case 2, one per LOQ fraction)
sd = c(0.05, 1.4) #Lowest and highest std for NI based on locations to be used when simulation error term
slope = c(log(1.01), log(1.1)) # Value of the slope in log-scale for a 1% and 5% yearly increase in original scale.
mean = 0 #Error term are log-normal with mean 0, var = sd^2
i=1 #Index to use for looping over repetition in a certain simulation
repititions = 100 #Number of reps of each simulation

true_coefs_2 <- matrix(nrow=repititions*simulations, ncol=6) #Creates matrices to save value of variables for each method
tobit_coefs_2 <- matrix(nrow=repititions*simulations, ncol=6)
museum_coefs_2 <- matrix(nrow=repititions*simulations, ncol=6)


start_time <- Sys.time()
  while(i<101){ #Repeating a simulation the number of times decided by the variable repititions above
    epsilon_errors <- rlnorm(N, mean, sd[2]) #Get the error terms for each response variable
    year1 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[1,1]) #Simulates between-years errorterms
    year2 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[2,1])
    year3 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[3,1])
    year4 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[4,1])
    year5 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[5,1])
    year6 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[6,1])
    year7 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[7,1])
    year8 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[8,1])
    year9 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[9,1])
    year10 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[10,1])
    year11 <- rlnorm(N/nr_of_covariates, mean, yearly_variances[11,1])
    random_effects <- c(year1, year2, year3, year4, year5, year6, year7, year8, year9, year10, year11)
    n <- length(epsilon_errors) #Number of error terms noted as n
    LOQ <- sort(epsilon_errors*random_effects)[0.8*n] #Decides the LOQ for the specific repitition
    predictors <- exp(kovariat*slope[2]) * epsilon_errors*random_effects #Intercept = 0, Slope = selected value of slope vector,  this builds our model to simulate from
    dataset_censored <- ifelse(predictors<LOQ, -LOQ, predictors) #Censores values under LOQ
    
    #Tar fram tobitmodellen
    y <- log(abs(dataset_censored))
    cens <- dataset_censored < 0
    x <- (kovariat)
    Z <- matrix(rep_len(1,N), ncol=1)
    cluster <- as.numeric(factor(kovariat))
    fitted <- mixcens(y, cens, x, init = c(0,0,1,1))
    tobit_coefs_2[i,1] <- fitted$beta[1]
    tobit_coefs_2[i,2] <- fitted$beta[2]
    tobit_coefs_2[i,3] <- sqrt(fitted$varFix[1,1])
    tobit_coefs_2[i,4] <- sqrt(fitted$varFix[2,2])
    tobit_coefs_2[i,5] <- i
    tobit_coefs_2[i,6] <- LOQ_fraction
    
    #Den verkliga ocensurerade modellen
    true_model <- lm(log(predictors)~kovariat)
    model_summary <- summary(true_model)
    true_coefs_2[i,1] <- true_model$coefficients[1]
    true_coefs_2[i,2] <- true_model$coefficients[2]
    true_coefs_2[i,3] <- model_summary$coefficients[1,2]
    true_coefs_2[i,4] <- model_summary$coefficients[2,2]
    true_coefs_2[i,5] <- i
    true_coefs_2[i,6] <- LOQ_fraction
    
    #Modellen enligt museet
    museets_obs <- ifelse(dataset_censored<0, log(abs(LOQ)/sqrt(2)), log(abs(dataset_censored)))
    museum_model <- lm(museets_obs~kovariat)
    museum_summary <- summary(museum_model)
    museum_coefs_2[i,1] <- museum_model$coefficients[1]
    museum_coefs_2[i,2] <- museum_model$coefficients[2]
    museum_coefs_2[i,3] <- museum_summary$coefficients[1,2]
    museum_coefs_2[i,4] <- museum_summary$coefficients[2,2]
    museum_coefs_2[i,5] <- i
    museum_coefs_2[i,6] <- LOQ_fraction
    
    i = i+1
  }
end_time <- Sys.time()
time_of_loop <- end_time - start_time



### Saves dataframe

#Named as (Slope, errorterms, randomeffects)
#creates a df of all coefficents
coefs <- cbind(true_coefs_2, tobit_coefs_2, museum_coefs_2) %>%
  as.data.frame() %>%
  select(V1, V2,V3, V4, V7, V8, V9, V10, V13, V14, V15, V16, V17, V18) %>%
  rename('True Intercept' = V1) %>%
  rename('True Beta' = V2) %>%
  rename('True Intercept sd' = V3) %>%
  rename('True Beta sd' = V4) %>%
  rename('Tobit Intercept' = V7) %>%
  rename('Tobit Beta' = V8) %>%
  rename('Tobit Intercept sd' = V9) %>%
  rename('Tobit Beta sd' = V10) %>%
  rename('Museum Intercept' = V13) %>%
  rename('Museum Beta' = V14) %>%
  rename('Museum Intercept sd' = V15) %>%
  rename('Museum Beta sd' = V16) %>%
  rename('simulation' = V17) %>%
  rename('Limit fraction' = V18) %>%
  mutate('Std' = sd[2]) %>%
  mutate('Slope' = slope[2])


intercepts <- coefs %>%
  dplyr::select(`True Intercept`, `Tobit Intercept`, `Museum Intercept`, simulation, `Limit fraction`, Std, Slope) %>%
  gather(key="Method", value="Intercept", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage"), sep=" ") %>%
  dplyr::select(simulation, Method, Intercept,`Limit fraction`, Std, Slope )


beta <- coefs %>%
  dplyr::select(`True Beta`, `Tobit Beta`, `Museum Beta`, simulation, `Limit fraction`, Std, Slope ) %>%
  gather(key="Method", value="Beta", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage"), sep=" ") %>%
  dplyr::select(simulation, Method, Beta, `Limit fraction`, Std, Slope)

intercepts_sd <- coefs %>%
  dplyr::select(`True Intercept sd`, `Tobit Intercept sd`, `Museum Intercept sd`, simulation, `Limit fraction`, Std, Slope ) %>%
  gather(key="Method", value="Intercept sd", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage", "Garb2"), sep=" ") %>%
  dplyr::select(simulation, Method, `Intercept sd`,`Limit fraction`, Std, Slope )

beta_sd <- coefs %>%
  dplyr::select(`True Beta sd`, `Tobit Beta sd`, `Museum Beta sd`, simulation, `Limit fraction`, Std, Slope ) %>%
  gather(key="Method", value="Beta sd", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage", "Garb2"), sep=" ") %>%
  dplyr::select(simulation, Method, `Beta sd`,`Limit fraction`, Std, Slope )

intercept_beta <- inner_join(intercepts, beta, by = c('simulation','Method', 'Limit fraction', 'Std', 'Slope'))
s0.05_r0.05 <- inner_join(intercept_beta, intercepts_sd, by = c('simulation','Method', 'Limit fraction', 'Std', 'Slope')) %>% 
  mutate('Random Effects' = 'Low') #NEED TO CHANGE

#Named as (Slope, errorterms, randomeffects)
df_1.1_0.5_Low_NEG <- inner_join(s0.05_r0.05, beta_sd, by = c('simulation','Method', 'Limit fraction', 'Std', 'Slope')) %>% 
  select(simulation, `Limit fraction`, Slope, Std, Method, Intercept, Beta, `Intercept sd`, `Beta sd`, `Random Effects`)
save(df_1.1_0.5_Low_NEG, file="slope1.1_res0.5_random_Low_NEG.Rda")