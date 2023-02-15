rm(list = ls())

load("hw3_data.RData")
suppressWarnings()




# EXERCISE 3 --------------------------------------------------------------

# Load packages.
packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)


# Input.
k <- 5 # Dimensions
n0 <- 80  # Sample size 0
n1 <- 90  # Sample size 1
alpha <- .05 # significance level


# Reproducibility.
set.seed(324)


# Take random sample from a multi-normal distribution.
Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}


# Generate k random means.
generate_mean <- function(k){
  
  return(sample(1:10, k))
}


# Generate random Covariance matrix (k*k).
generate_sigma <- function(k){
  sym_mat <- matrix(1, ncol = k, nrow = k, byrow = TRUE)
  diag_mat <- diag(sample(0:20), k)
  sigma <- diag_mat + sym_mat
  
  return(sigma)
}


# Sigmoid Function.
sigmoid <- function(x, theta){
  n_features <- length(theta)
  n <- exp(theta[1] + sum(x[-length(x)]*theta[2:n_features]))
  
  return(n / (1+n))
}



mu <- generate_mean(k)
sigma <- generate_sigma(k)

x <- Take_sample_normal(n = n0, mu=mu, sigma = sigma, label = 0)
z <- Take_sample_normal(n = n1, mu=mu, sigma = sigma, label = 1)
u <- rbind(x, z) # Actual data




# Friedman Procedure ------------------------------------------------------

P <- 1000
?predict
Friedman_procedure <- function(P){
  
  kolm_t <- rep(NA, P)

  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
    z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  mann_t <<- mann_t
  kolm_t <<- kolm_t
  

  
  hist(kolm_t, main = "
       Kolmogorov-Smirnov statistic \n distribution under H_0",
       col = "skyblue", border = "white", breaks= 30)
  p_kolm <<- quantile(kolm_t , 1 - alpha)
  abline(v = p_kolm , col = "red" , lty = 3 , lwd = 2)
  box()
  

}

Friedman_procedure(P)
par(mfrow = c(1,1))

# Simulation to get info about ALPHA --------------------------------------

alpha_info <- function(P, percentile_kolm){
  prop_rej <- rep(NA , P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu=mu, sigma = sigma, label = 0)  
    z <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u <- rbind(x, z) # Combine the data
    true_coef <- glm(label ~ ., data = u)$coefficients
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = true_coef)
    z_scores <- apply(z, MARGIN = 1, sigmoid, theta = true_coef)
    true_kolm <- ks.test(x_scores ,z_scores, alternative = "two.sided")$statistic
    # Perform the Test
    prop_rej[i] <- true_kolm < percentile_kolm

  }


  
  data = as.data.frame(prop_rej)
  colnames(data) <- c("Kolmogorov-Smirnov")
  
  return(data)
}


data <- alpha_info(P = P, percentile_kolm = p_kolm)


barplot(prop.table(table(data$`Kolmogorov-Smirnov`)), col = c("red" , "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using KS statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))

#barplot(prop.table(table(data$Mann)), col = c("red", "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using Mann Statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))



?ks.test()

# Simulation to get info about POWER --------------------------------------

# Distance functions.
beta_q <- function(sigma1, sigma2){
  trace_1 <- tr(sigma)
  trace_2 <- tr(sigma2)
  trace_3 <- tr(sqrtm(sqrtm(sigma1)$B %*% sigma2 %*% sqrtm(sigma1)$B)$B)
  
  return(trace_1 + trace_2 - 2*trace_3)
}

Wasserstein_distance <- function(mu1, mu2, sigma1, sigma2){
  norma <- norm(mu1-mu2, type = "2")**2
  beta_quadro <- beta_q(sigma1, sigma2)
  
  return(norma + beta_quadro)
}


n0 <- 100
n1 <- 100

power_info <- function(P, k, percentile_ks, percentile_mann, increment){
  mu2 <- mu + increment 
  sigma2 <- sigma
  distance <- Wasserstein_distance(mu, mu2, sigma, sigma2)
  prop_ks <- rep(NA, P)
  prop_mann <- rep(NA, P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    true_coef <- glm(label ~ ., data = u)$coefficients
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = true_coef)
    z_scores <- apply(z, MARGIN = 1, sigmoid, theta = true_coef)
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    true_mann <- wilcox.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_ks[i] <- true_kolm > percentile_ks
    prop_mann[i] <- true_mann > percentile_mann
  }
  data <- c(distance, mean(prop_ks), mean(prop_mann))
  names(data) <- c('Distance', 'K_S', 'Mann')
  
  return(data)
}

k <- 5
mu <- 1:k
sigma <- generate_sigma(length(mu))
power_info(100, k = length(mu), percentile_ks = p_kolm, percentile_mann  = p_mann, increment = 10)





# Plot
results <- matrix(NA, ncol = 3, nrow = 0)

k <- 10
P <- 1000
M <- seq(0,1, .05)
mu <- generate_mean(k)
sigma <- generate_sigma(length(mu))

for(l in M){
  a <- power_info(P, k = length(mu), percentile_ks = p_kolm,
                  percentile_mann = p_mann, increment = l)
  results <- rbind(results, a)
  
}

results <- as.data.frame(results)

par(mfrow = c(1,1))
plot(results$Distance, results$K_S, type = 'l', main='Relation between distance and power', lwd=2, col='skyblue', xlab='Distance', ylab='Power')
points(results$Distance, results$Mann, type = 'l', main='Relation between distance and power', lwd=2, col='lightgreen')
legend('bottomright', legend = c('K-S', 'Mann'), col= c('skyblue', 'lightgreen'), lwd=3)

