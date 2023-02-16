

# EXERCISE 3 --------------------------------------------------------------

rm(list = ls())

load("hw3_data.RData")
suppressWarnings()




# Load packages.
packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)






# Reproducibility.
set.seed(234)


# Take random sample from a multi-normal distribution.
Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}


Friedman_procedure <- function(P, xdata){
  kolm_t <- rep(NA, P)
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label = 1) # Under H_0
    u_p <- rbind(xdata, z_p)
    glm_model <-glm(label ~ ., data = u_p)
    x_scores <- predict(glm_model , xdata[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  return(kolm_t)
}


alpha_info <- function(M, P){
  prop_rej <- rep(NA, M)
  p_values <- rep(NA , M)
  
  for(i in 1:M){
    x_p <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)  
    z_p <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1) 
    u_p<- rbind(x_p, z_p) 
    
    kk <- Friedman_procedure(P, xdata  = x_p, zdata = z_p) 
    glm_model <-glm(label ~ ., data = u_p)
    
    x_scores <- predict(glm_model , x_p[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_rej[i] <- true_kolm < quantile(kk, 1 - alpha)
  }
  
  return(prop_rej)
}



# Input.
k <- 5 
P <- 100
M <- 100
n0 <- 80  
n1 <- 90  
alpha <- .05 


mu <- rep(0, k)
sigma <- diag(1, k)






# alpha_Friedmann <- rep(NA, length(n0s))
# for(i in 1:length(n0s)){
#   n0 <- n0s[i]
#   alpha_Friedmann[i] <- 1-mean(alpha_info(P))
# }



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



power_info <- function(M, P, n0, n1, mu, mu2, sigma, sigma2){
  
  prop_ks <- rep(NA, M)
  dist <- Wasserstein_distance(mu, mu2, sigma, sigma2)
  
  for(i in 1:M){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu= mu2, sigma = sigma2, label = 1)
    u <- rbind(x, z) 
    
    p_model <- glm(label ~ ., data = u)
    x_scores <- predict(p_model, x)
    z_scores <- predict(p_model, z) 
    
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    kk <- Friedman_procedure(P, x)
    prop_ks[i] <- true_kolm > quantile(kk, 1-alpha)
  }
  data <- c(dist, mean(prop_ks))
  
  return(data)
}




k <- 100
P <- 100
M <- 10
n0 <- 80  
n1 <- 90  
alpha <- .05 

ns <- c(20, 40, 60, 80, 100)
ks <- c(5, 20, 50, 100)
colors <- c('skyblue', 'lightgreen', 'orchid', 'darkorange')


mu <- rep(0, k)
sigma <- diag(1, k)


grid <- seq(0, .4, .1)



distance_k <- matrix(NA, 0, 2)

for(i in grid){
  print(i)
  power_k <- power_info(M, P, n0, n1, mu, mu+i, sigma, sigma)
  distance_k <- rbind(distance_k, power_k)
}

plot(distance_k, xlim=c(0,8), main = "Relation between distance and power\n for different k's", xlab='k', ylab='power', las=1, type='l', lwd=3, col=colors[k])

legend('bottomright', col=colors, legend=ks, lwd=3)

Wasserstein_distance(mu, mu+.3, sigma, sigma)
