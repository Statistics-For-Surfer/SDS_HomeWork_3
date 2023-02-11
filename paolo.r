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
set.seed(123)


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

Friedman_procedure <- function(P){
  
  kolm_t <- rep(NA, P)
  mann_t <- rep(NA, P)
  
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
    z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    mann_t[i] <- wilcox.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  par(mfrow = c(1,2))
  
  hist(kolm_t, main = "
       Kolmogorov-Smirnov statistic \n distribution under H_0",
       col = "skyblue", border = "white", breaks= 20)
  p_kolm <- quantile(kolm_t , 1 - alpha)
  abline(v = p_kolm , col = "red" , lty = 3 , lwd = 2)
  box()
  
  
  hist(mann_t, main = "
       Mann statistic distribution \n under H_0",
       col = "lightgreen", border = "white", breaks= 30)
  p_mann <- quantile(mann_t , 1 - alpha)
  abline(v = p_mann  , col = "red" , lty = 3 , lwd = 2)
  box()

}

Friedman_procedure(P)



# Simulation to get info about ALPHA --------------------------------------

alpha_info <- function(P, percentile_kolm, percentile_mann){
  prop_rej <- matrix(NA , P , 2)
  prop_rej_kolm <- rep(NA , P)
  prop_rej_mann <- rep(NA , P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu=mu, sigma = sigma, label = 0)  
    z <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u <- rbind(x, z) # Combine the data
    true_coef <- glm(label ~ ., data = u)$coefficients
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = true_coef)
    z_scores <- apply(z, MARGIN = 1, sigmoid, theta = true_coef)
    true_kolm <- ks.test(x_scores ,z_scores, alternative = "two.sided")$statistic
    true_mann <- wilcox.test(x_scores, z_scores, alternative = "two.sided")$statistic
    # Perform the Test
    prop_rej[i,1] <- true_kolm < percentile_kolm
    prop_rej[i,2] <- true_mann < percentile_mann
  }


  
  data = as.data.frame(prop_rej)
  colnames(data) <- c("Kolmogorov-Smirnov", "Mann")
  
  return(data)
}


data <- alpha_info(P = P, percentile_kolm = p_kolm, percentile_mann = p_mann)


barplot(prop.table(table(data$`Kolmogorov-Smirnov`)), col = c("red" , "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using KS statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))

barplot(prop.table(table(data$Mann)), col = c("red", "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using Mann Statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))





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



power_info <- function(P, k, percentile_ks, percentile_mann){
  mu2 <- mu + 1   # insert by hand
  sigma2 <- sigma
  distance <- Wasserstein_distance(mu, mu2, sigma, sigma2)
  prop_ks <- rep(NA, P)
  prop_mann <- rep(NA, P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n0, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    true_coef <- glm(label ~ ., data = u)$coefficients
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = true_coef)
    z_scores <- apply(z, MARGIN = 1, sigmoid, theta = true_coef)
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    true_mann <- wilcox.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_ks[i] <- true_kolm < percentile_ks
    prop_mann[i] <- true_mann < percentile_mann
  }
  data <- c(distance, mean(prop_ks), mean(prop_ks))
  names(data) <- c('Distance', 'K-S', 'Mann')
  
  return(data)
}


mu <- c(1, 2, 3, 4, 5)
sigma <- generate_sigma(length(mu))
power_info(1000, k = length(mu), percentile_ks = p_kolm, percentile_mann  = p_mann)





