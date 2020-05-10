#--------------------------------------Code for Simulation-----------------------------------------------------
#@author: Anton Holm
# This is a simulation used for my Bachelor's Thesis in Mathematical Statistics at Stockholm University
# Spring, 2020
#--------------------------------------Code for Simulation-----------------------------------------------------

library(lmec)
library(tidyverse)
load("yearly_variances.Rda")

mixcens <- function(y, cens, x, init = c(0, 0, 1, 1)){
  library(mnormt)
  l_i <- function(data, alpha, beta, sigma2b, sigma2e){
    n <- nrow(data)
    n_c <- sum(data$cens)
    Sigma11 <- diag(sigma2e, n_c) + matrix(sigma2b, nrow = n_c, ncol = n_c)
    Sigma22 <- diag(sigma2e, n - n_c) + matrix(sigma2b, nrow = n - n_c, ncol = n - n_c)
    Sigma21 <- matrix(sigma2b, nrow = n - n_c, ncol = n_c)
    Sigma12 <- matrix(sigma2b, nrow = n_c, ncol = n - n_c)
    mu <- alpha + beta * data$x
    mu1 <- mu[1:n_c]
    y1 <- data$y[1:n_c]
    mu2 <- mu[(n_c + 1):n]
    y2 <- data$y[(n_c + 1):n]
    ll <- 0
    if ((n_c > 0) & (n_c < n)){
      ll <- ll + log(pmnorm(y1, 
                            mean = mu1 + as.numeric(Sigma12 %*% solve(Sigma22) %*% (y2 - mu2)),
                            varcov = Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma21
      ))
    }
    if (n_c < n){
      ll <- ll + log(dmnorm(y2, mean = mu2, varcov = Sigma22))
    }
    if (n_c == n){
      ll <- log(pmnorm(y1, mean = mu1, varcov = Sigma11))
    }
    ll
  }
  l <- function(data, pars){
    data <- arrange(data, x, -cens)
    alpha <- pars[1]
    beta <- pars[2]
    sigma2b <- pars[4]
    sigma2e <- pars[3]
    loglik <- data %>% group_by(x) %>% 
      nest() %>% 
      mutate(ll = map2_dbl(x, data, ~l_i(data = mutate(.y, x = .x), alpha, beta, sigma2b, sigma2e))) %>% 
      pull(ll) %>% sum()
    - loglik
  }
  data <- tibble(y = y, x = x, cens = cens)
  fitted <- optim(par = init, fn = function(x) l(data = data, pars = c(x[1], x[2], exp(x[3]), exp(x[4]))), hessian = TRUE)
  list(beta = fitted$par[1:2], varFix = solve(fitted$hessian)[1:2, 1:2])
}

set.seed(19931031)

kovariat <- rep(0:10-mean(0:10), each = 12) #Creates covariates
nr_of_covariates <- length(unique(kovariat)) # Number of unique covariates
N = length(kovariat) #Denotes the number of datapoints per simulation to N
mean = 0 #Error term are log-normal with mean 0, var = sd^2
i=0 #Index to use for looping over repetition in a certain simulation
k=1
mixcens_1.01 <- matrix(nrow=100, ncol=3) #Creates matrices to save value of variables for each method
lmec_1.01 <- matrix(nrow=100, ncol=3)
subs_1.01 <- matrix(nrow=100, ncol=3)

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

start_time <- Sys.time()
  while(i<1){ #Repeating a simulation the number of times decided by the variable repititions above
    LOQ <- max(sort(epsilon_errors*random_effects)[i*n], 0, na.rm = TRUE) #Decides the LOQ for the specific repitition
    predictors <- exp(kovariat*log(1.01)) * epsilon_errors*random_effects #Intercept = 0, Slope = selected value of slope vector,  this builds our model to simulate from
    dataset_censored <- ifelse(predictors<LOQ, -LOQ, predictors) #Censores values under LOQ
    
    #Tar fram mixcens
    y <- log(abs(dataset_censored))
    cens <- dataset_censored < 0
    x <- (kovariat)
    fit_mix <- try(mixcens(y, cens, x, init = c(0,0,1,1)), silent = TRUE)
    mixcens_1.01[k,1] <- try(fit_mix$beta[2], silent = TRUE)
    mixcens_1.01[k,2] <- try(sqrt(fit_mix$varFix[2,2]), silent = TRUE)
    mixcens_1.01[k,3] <- i
    
    #Tar fram tobitmodellen
    X <- cbind(rep_len(1,N), kovariat)
    Z <- matrix(rep_len(1,N), ncol=1)
    cluster <- as.numeric(factor(kovariat))
    fitted <- lmec(y, cens, X, Z, cluster, maxstep = 20, method = "ML")
    lmec_1.01[k,1] <- fitted$beta[2]
    lmec_1.01[k,2] <- sqrt(fitted$varFix[2,2])
    lmec_1.01[k,3] <- i
    

    #Modellen enligt museet
    museets_obs <- ifelse(dataset_censored<0, log(abs(LOQ)/sqrt(2)), log(abs(dataset_censored)))
    museum_model <- lm(museets_obs~kovariat)
    museum_summary <- summary(museum_model)
    subs_1.01[k,1] <- museum_model$coefficients[2]
    subs_1.01[k,2] <- museum_summary$coefficients[2,2]
    subs_1.01[k,3] <- i
    
    i = i+0.01
    k = k+1
    print(i)
  }
end_time <- Sys.time()
time_of_loop <- end_time - start_time

df_cens_1.01 <- as.data.frame(mixcens_1.01) %>% mutate(lmec_estim = lmec_1.01[,1]) %>% mutate(lmec_sd = lmec_1.01[,2]) %>% 
  mutate(subs_estim = subs_1.01[,1]) %>% mutate(subs_sd = subs_1.01[,2]) %>% mutate(slope = log(1.01)) %>% rename(mixcens_estim = V1, mixcens_sd = V2, censored = V3)

save(df_cens_1.01, file = "df_cens_1.01.Rda")

mix_df_1.01 <- mixcens_1.01 %>% as.data.frame() %>% mutate(V1 = as.numeric(as.character(V1))) %>% mutate(V2 = as.numeric(as.character(V2))) %>% mutate(V3 = as.numeric(as.character(V3)))
lmec_df_1.01 <- lmec_1.01 %>% as.data.frame()
subs_df_1.01 <- subs_1.01 %>% as.data.frame()

cens_df_1.01_low <- full_join(mix_df_1.01, lmec_df_1.01) %>% full_join(subs_df_1.01) %>% 
  rename(Beta_Hat = V1, Sd = V2, Censoring = V3)

save(cens_df_1.01_low, file = "../cens_df_1.01_low.Rda")
