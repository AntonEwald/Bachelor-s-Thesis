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