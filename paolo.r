rm(list = ls())
load("hw3_data.RData")
suppressWarnings()


####
k <-  5 # Dimensions
n0 <- 80  # Sample size 0
n1 <- 90  # Sample size 1
alpha <- .05 # significance level
#### Generate Data from a multivariate normal
library(MASS)






##### Esercise 3
set.seed(123) # Reproducibility

Take_sample_normal <- function(n, mu,sigma,label){
  x <- mvrnorm(n,mu,sigma)
  x <- as.data.frame(cbind(x,label))
  colnames(x[length(colnames(x))]) <-"label" 
  return(x)
}


generate_mean <- function(k){
  return(sample(1:7, k))
}

generate_sigma <- function(k){
  sym_mat <- matrix(1, ncol = k, nrow = k , byrow = TRUE)
  diag_mat <- diag(sample(0:20),k)
  # Sum the diagonal and symmetric matrices
  sigma <- diag_mat + sym_mat
  return(sigma)
  
}


sigma <- generate_sigma(k)
x <- Take_sample_normal(n = n0,mu=mu,sigma = sigma,label = 0)
z <- Take_sample_normal(n = n1,mu=mu,sigma = sigma,label = 1)
u <- rbind(x,z) # Actual data


######
# Sigmoid Function
sigmoid <- function(x, theta){
  n_features <- length(theta)
  n <- exp(theta[1] + sum(x*theta[2:n_features]))
  return(n / (1+n))}######



#### Friedmann Procedure

P <- 10000 # simulation size
kolm_t <- rep(NA , P)
mann_t <- rep(NA, P)

for(i in 1:P){
  z_p <- Take_sample_normal(n1,mu,sigma,label =1) # Under H_0
  u_p <- rbind(x,z_p)
  glm_coef <- glm(label ~. , data = u_p)$coefficients
  x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = glm_coef)
  z_scores <- apply(z_p ,MARGIN = 1 , sigmoid ,theta = glm_coef) 
  kolm_t[i] <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
  mann_t[i] <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
}

hist(kolm_t, main = "
     Kolmogorov-Smirnov statistic distribution \n under H_0",
     col = "skyblue" , border = "white")
box()
p_kolm <- quantile(kolm_t , 1 - alpha)
abline(v = p_kolm , col = "red" , lty = 3 , lwd = 2)

p_mann <- quantile(mann_t , 1 - alpha)
hist(mann_t, main = "
     Mann statistic distribution \n under H_0",
     col = "lightgreen" , border = "white")
abline(v = p_mann  , col = "red" , lty = 3 , lwd = 2)

box()
########


# Simulation to get info about alpha

info_alpha <- function(P,percentile_kolm,percentile_mann){
  prop_rej <- matrix(NA , P , 2)
  prop_rej_kolm <- rep(NA , P)
  prop_rej_mann <- rep(NA , P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0,mu=mu,sigma = sigma,label = 0)  
    z <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u <- rbind(x,z) # Combine the data
    true_coef <- glm(label ~. , data = u)$coefficients
    x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = true_coef)
    z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = true_coef)
    true_kolm <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
    true_mann <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
    # Perform the Test
    prop_rej[i,1] <- true_kolm < percentile_kolm
    prop_rej[i,2] <- true_mann < percentile_mann
  }
  
  data = as.data.frame(prop_rej)
  colnames(data) <- c("Kolmogorov-Smirnov","Mann")
  
  return(data)
}



data <- info_alpha(P = P , percentile_kolm  = p_kolm , percentile_mann  = p_mann)
barplot(table(data$`Kolmogorov-Smirnov`), col = c("red" , "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using KS statistic", names.arg = c("Reject" , "Accept"))

barplot(table(data$Mann), col = c("red", "blue"),  main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using Mann Statistic",names.arg = c("Reject" , "Accept"))

# Simulation to get info about power [TODO]

# Since we are dealing with a family of multivariate normal distributions we can get the formula dor the distance analitycally
library(pracma)
library(psych)



beta_q <- function(sigma1,sigma2){
  traccia1 <- tr(sigma)
  traccia2 <- tr(sigma2)
  traccia3 <- tr(sqrtm(sqrtm(sigma1) %*% sigma2 %*% sqrtm(sigma1)))
  return(traccia1 + traccia2 -2*traccia3)
}
Wasserstein_distance <- function(mu1,mu2,sigma1,sigma2){
  norma <- norm(mu1-mu2 , type = "2")**2
  beta_quadro <- beta_q(sigma1,sigma2)
  return(norma + beta_quadro)
  
}

### GET THE INFO
mu

power_info <- function(P,k, percentile_ks , percentile_mann){
  mu2 <- c(1,5)
  sigma2 <- sigma
  distance <- Wasserstein_distance(mu,mu2,sigma,sigma2)
  prop_ks <- rep(NA , P)
  prop_mann <- rep(NA , P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu , sigma = sigma,label = 0)
    z <- Take_sample_normal(n = n0, mu=  mu2, sigma = sigma2,label = 1)
    u <- rbind(x,z) # Combine the data
    true_coef <- glm(label ~. , data = u)$coefficients
    x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = true_coef)
    z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = true_coef)
    true_kolm <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
    true_mann <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
    prop_ks[i] <- true_kolm < percentile_ks
    prop_mann[i] <- true_mann < percentile_mann
  }
  data <- c(distance, mean(prop_ks),mean(prop_ks))
  return(data)
}




mu <- c(1,2)
sigma <- generate_sigma(2)
power_info(100,k = length(mu), percentile_ks = p_kolm ,percentile_mann  = p_mann)




mu
sigma
p_kolm
p_mann








P = 3
aa <- info_power(P , p_kolm , p_mann ,  theta_space)
aa

#################






abline(v = quantile(kolm_t, 1 - alpha), col = "red", lwd = 3)
abline(v= true_kolm , col = "blue")

hist(mann_t)
abline(v = quantile(mann_t, 1 - alpha), col = "red")
abline(v= true_mann , col = "blue")
########
# Simulation to get the power of the test




########



