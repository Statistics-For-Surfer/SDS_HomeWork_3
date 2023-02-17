rm(list = ls())
#load packages
packs <- c('MASS', 'pracma', 'psych','caTools','e1071',"dgof")
lapply(packs, require, character.only = TRUE)
###########

k <- 5 # Dimensions
n0 <- 70  # Sample size 0
n1 <- 70  # Sample size 1
alpha <- .05 # significance level
set.seed(324)# reproducibility
P <- 1000
acc_rej_col <- c("#d3305d" , "#ABCDEF")
###########

# Take random sample from a multi-normal distribution.
Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}



mu <- rep(0,k)
sigma <- diag(1 ,k)

############



Friedman_procedure <- function(P,xdata,zdata){
  kolm_t <- rep(NA, P)
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(xdata, z_p)
    glm_model <-glm(label ~ ., data = u_p)
    x_scores <- predict(glm_model , xdata[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  return(kolm_t)
  
}

x_data <- Take_sample_normal(n0,mu,sigma,0)

Friedman_procedure_2 <- function(P,x_data){
  x_fri <- x_data
  kolm_t <- rep(NA, P)
  labels <- c(rep(0,n1),rep(1,n1))
  for(i in 1:P){
    idx <- sample(x = 1:(n0+n1), n0+n1)
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    x_fri$label <- labels[idx[1:n0]]
    z_p$label <- length(labels[idx[(n0+1):(n0+n1)]])
    u_p <- as.data.frame(rbind(x_fri,z_p))
    glm_f <- glm(label~., data = u_p)
    scores <- predict(glm_f ,u_p[,1:k])
    kolm_t[i] <- ks.test(scores[u_p$label == 0] , scores[u_p$label == 1])$statistic

  }
  return(kolm_t)
  
  
}
############

P <- 100
alpha_info <- function(P){
  prop_rej <- rep(NA, P)
  p_values <- rep(NA , P)
  
  for(i in 1:P){
    x_p <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)  
    z_p <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u_p<- rbind(x_p, z_p) # Combine the data
    
    kk <- Friedman_procedure_2(P = P, x_data  = x_p ) # As mavi said
    glm_model <-glm(label ~ ., data = u_p)
    
    x_scores <- predict(glm_model , x_p[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_rej[i] <- true_kolm < quantile(kk , 1 -alpha)
    
  }
  

  return(prop_rej)
}


data <- alpha_info(P)



t <- proportions(table(data))
barplot(t, col = acc_rej_col , main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using KS statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))

n0s <- c(20,50,70,100)
k <- 10
P <- 150
mu <- rep(1,k)
sigma <- diag(1,k)
alpha_Friedmann <- rep(NA ,length(n0s) )
for(i in 1:length(n0s)){
  n0 <- n0s[i]
  alpha_Friedmann[i] <- 1-mean(alpha_info(P))

}

plot(n0s,alpha_Friedmann,type = "o",ylim = c(0,.4), col = "gold",lwd = 3)

ks <- c(5,20,70,100)
n0 <- 70
alpha_Friedmann_2 <- rep(NA ,length(ks))
P = 150
for(i in 1:length(ks)){
  k <- ks[i]
  mu <- rep(1,k)
  sigma <- diag(1,k)
  alpha_Friedmann_2[i] <- 1-mean(alpha_info(P))
}

plot(ks,alpha_Friedmann_2,type = "o",ylim = c(0,.4), col = "gold",lwd = 3)
#############
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


#############

power_info <- function(P, k, increment){
  
  prop_ks <- rep(NA, P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    p_model <- glm(label ~ ., data = u)
    x_scores <- predict(p_model, x)
    z_scores <- predict(p_model, z) 
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    kk <- Friedman_procedure(P,x,z)
    prop_ks[i] <- true_kolm > quantile(kk,1-alpha)
  }
  data <-  mean(prop_ks)
  
  
  return(data)
}

?glm
##########
P <- 100
k <- 5
mu <- 1:k
sigma <- diag(1,k)
max_in <- 1.5
ll <- seq(from = 0 , to = max_in , by = 0.3 )
powers <- c()
distance <- c()
for (i in (ll)){
  mu2 <- mu + i 
  print(mu2)
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers,power_info(100, k = length(mu),increment = i))
}
plot(distance,powers,col = 'darkblue' , type = "o" , lwd = 2)

Wasserstein_distance(mu, mu +1.5 , sigma, sigma)
####

###


###

#######
#######
