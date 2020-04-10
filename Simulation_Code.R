#--------------------------------------Code for Simulation-----------------------------------------------------
#@author: Anton Holm
# This is a simulation used for my Bachelor's Thesis in Mathematical Statistics at Stockholm University
# Spring, 2020
#--------------------------------------Code for Simulation-----------------------------------------------------

library(lmec)
library(tidyverse)

set.seed(19931031)

kovariat <- rep(0:10-mean(0:10), each = 12) #Creates covariates
N = length(kovariat) #Denotes the number of datapoints per simulation to N
LOQ_fraction2 = c(0.3, 0.6) #The proportion of data to be censored when intercept and slope equals zero
simulations = length(LOQ_fraction2) # How many simulations with different input (in this case 2, one per LOQ fraction)
sd = c(0.05, 1.4) #Lowest and highest std for NI based on locations to be used when simulation error term
slope = c(log(1.01), log(1.05)) # Value of the slope in log-scale for a 1% and 5% yearly increase in original scale.
mean = 0 #Error term are log-normal with mean 0, var = sd^2
i=1 #Index to use for looping over repetition in a certain simulation
k=1 #Index for different simulations (i.e different limits)
repititions = 100 #Number of reps of each simulation

true_coefs_2 <- matrix(nrow=repititions*simulations, ncol=4) #Creates matrices to save value of variables for each method
tobit_coefs_2 <- matrix(nrow=repititions*simulations, ncol=4)
museum_coefs_2 <- matrix(nrow=repititions*simulations, ncol=4)


start_time <- Sys.time()
while(k<simulations+1){ #Simulations for different LOQ
  reps = k*repititions #Increases max value of i in order to get more repetitions for the next simulation
  i = (k-1)*repititions + 1 #Force i to start at a suitable value when changing to another simulation
  LOQ_fraction = LOQ_fraction2[k] #Pick the fraction of censored values (k:th value from above vector)
  while(i<reps+1){ #Repeating a simulation the number of times decided by the variable repititions above
    epsilon_errors <- rlnorm(N, mean, sd[1]) #Get the error terms for each response variable
    n <- length(epsilon_errors) #Number of error terms noted as n
    LOQ <- sort(epsilon_errors)[LOQ_fraction*n] #Decides the LOQ for the specific repitition
    predictors <- kovariat*slope[2] + epsilon_errors #Intercept = 0, Slope = first value of slope vector,  this builds our model to simulate from
    dataset_censored <- ifelse(predictors<LOQ, -LOQ, predictors) #Censores values under LOQ
    
    #Tar fram tobitmodellen
    yL <- log(abs(dataset_censored))
    cens <- dataset_censored < 0
    X <- cbind(rep_len(1,N), kovariat)
    Z <- matrix(rep_len(1,N), ncol=1)
    cluster <- as.numeric(factor(kovariat))
    fitted <- lmec(yL, cens, X, Z, cluster, maxstep = 20, method = "ML")
    tobit_coefs_2[i,1] <- fitted$beta[1]
    tobit_coefs_2[i,2] <- fitted$beta[2]
    tobit_coefs_2[i,3] <- i
    tobit_coefs_2[i,4] <- LOQ_fraction
    
    #Den verkliga ocensurerade modellen
    true_model <- lm(log(predictors)~kovariat)
    true_coefs_2[i,1] <- true_model$coefficients[1]
    true_coefs_2[i,2] <- true_model$coefficients[2]
    true_coefs_2[i,3] <- i
    true_coefs_2[i,4] <- LOQ_fraction
    
    #Modellen enligt museet
    museets_obs <- ifelse(dataset_censored<0, abs(LOQ)/sqrt(2), dataset_censored)
    museum_model <- lm(museets_obs~kovariat)
    museum_coefs_2[i,1] <- museum_model$coefficients[1]
    museum_coefs_2[i,2] <- museum_model$coefficients[2]
    museum_coefs_2[i,3] <- i
    museum_coefs_2[i,4] <- LOQ_fraction
    
    i = i+1
  }
  k = k+1
}
end_time <- Sys.time()
time_of_loop <- end_time - start_time

#creates a df of all coefficents
coefs_slope0.05_res0.05 <- cbind(true_coefs_2, tobit_coefs_2, museum_coefs_2) %>%
  as.data.frame() %>%
  select(V1, V2, V5, V6, V9, V10, V11, V12) %>%
  rename('True Intercept' = V1) %>%
  rename('True Beta' = V2) %>%
  rename('Tobit Intercept' = V5) %>%
  rename('Tobit Beta' = V6) %>%
  rename('Museum Intercept' = V9) %>%
  rename('Museum Beta' = V10) %>%
  rename('simulation' = V11) %>%
  rename('Limit fraction' = V12) %>%
  mutate('Std' = sd[1]) %>%
  mutate('Slope' = slope[1])


intercepts_slope0.05_res0.05 <- coefs_slope0.05_res0.05 %>%
  dplyr::select(`True Intercept`, `Tobit Intercept`, `Museum Intercept`, simulation, `Limit fraction`, Std, Slope ) %>%
  gather(key="Method", value="Intercept", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage"), sep=" ") %>%
  dplyr::select(simulation, Method, Intercept,`Limit fraction`, Std, Slope )


beta_slope0.05_res0.05 <- coefs_slope0.05_res0.05 %>%
  dplyr::select(`True Beta`, `Tobit Beta`, `Museum Beta`, simulation, `Limit fraction`, Std, Slope ) %>%
  gather(key="Method", value="Beta", -c(simulation, `Limit fraction`, Std, Slope )) %>%
  separate(Method, c("Method", "Garbage"), sep=" ") %>%
  dplyr::select(simulation, Method, Beta, `Limit fraction`, Std, Slope)

df_slope0.05_res0.05 <- inner_join(intercepts_slope0.05_res0.05, beta_slope0.05_res0.05, by = c('simulation','Method', 'Limit fraction', 'Std', 'Slope'))


#_____________________________ Dataframes saved below!!!!_________________________________
df_slope0.01_res0.05 #The Dataframe of repetitions for simulations over each LOQ at slope = 1% increase and standard deviation = 0.05