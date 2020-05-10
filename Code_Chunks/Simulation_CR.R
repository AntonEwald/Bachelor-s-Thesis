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
mean = 0 #Error term are log-normal with mean 0, var = sd^2
i=1 #Index to use for looping over repetition in a certain simulation
repititions = 100 #Number of reps of each simulation


slope = log(1.05)
lmec_CR_1.05_0.05 <- matrix(nrow=100, ncol=4) #Creates matrices to save value of variables for each method
subs_CR_1.05_0.05 <- matrix(nrow=100, ncol=4)


start_time <- Sys.time()

  while(i<repititions+1){ #Repeating a simulation the number of times decided by the variable repititions above
    epsilon_errors <- rlnorm(N, mean, 0.05) #Get the error terms for each response variable
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
    predictors <- exp(kovariat*slope) * epsilon_errors*random_effects #Intercept = 0, Slope = selected value of slope vector,  this builds our model to simulate from
    dataset_censored <- ifelse(predictors<LOQ, -LOQ, predictors) #Censores values under LOQ

    #Tar fram tobitmodellen
    yL <- log(abs(dataset_censored))
    cens <- dataset_censored < 0
    X <- cbind(rep_len(1,N), kovariat)
    Z <- matrix(rep_len(1,N), ncol=1)
    cluster <- as.numeric(factor(kovariat))
    fitted <- lmec(yL, cens, X, Z, cluster, maxstep = 20, method = "ML")
    lmec_CR_1.05_0.05[i,1] <- fitted$beta[2]
    lmec_CR_1.05_0.05[i,2] <- sqrt(fitted$varFix[2,2])
    lmec_CR_1.05_0.05[i,3] <- i
    lmec_CR_1.05_0.05[i,4] <- log(slope)

        #Modellen enligt museet
    museets_obs <- ifelse(dataset_censored<0, log(abs(LOQ)/sqrt(2)), log(abs(dataset_censored)))
    museum_model <- lm(museets_obs~kovariat)
    museum_summary <- summary(museum_model)
    subs_CR_1.05_0.05[i,1] <- museum_model$coefficients[2]
    subs_CR_1.05_0.05[i,2] <- museum_summary$coefficients[2,2]
    subs_CR_1.05_0.05[i,3] <- i
    subs_CR_1.05_0.05[i,4] <- log(slope)
    
    
    i = i+1
    print(i)
  }
end_time <- Sys.time()
time_of_loop <- end_time - start_time


df_80_1.05_0.05 <- full_join(as.data.frame(lmec_CR_1.05_0.05),as.data.frame(subs_CR_1.05_0.05)) %>% 
  rename(Beta_Hat = V1, Standard_Error = V2, Simulation = V3, Slope = V4) %>% 
  mutate(Noise = 0.05) %>% mutate(Slope = exp(Slope))

save(df_80_1.05_0.05, file = "df_80_1.05_0.05.Rda")
