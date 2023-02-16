rm(list=ls())

load("hw3_data.RData")




# Scale dataset -----------------------------------------------------------
scale_datasets_list <- function(ls){
  scaled_list <- list()
  names <- names(ls)
  for(i in 1:length(ls)){
    
    #### Scaling by column.
    patient <- data.frame(apply(ls[[i]], 2, scale))
    colnames(patient) <- colnames(ls[[i]])
    scaled_list[[names[i]]] <- patient
  }
  return(scaled_list)
}

td_scale <- scale_datasets_list(td_data)
asd_scale <- scale_datasets_list(asd_data)



# Frequency of number observation  ----------------------------------------

# Count how many patients in the ASD group have time series of a given length
sort( table( sapply(td_scale, function(x) nrow(x)) ), decreasing = T )
sort( table( sapply(asd_scale, function(x) nrow(x)) ), decreasing = T )




# Plot of time series of each ROI for each patient -------------------------

library(manipulate)

manipulate(plot(asd_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "skyblue", lwd = 2, 
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(asd_data)), region = slider(1,116))


manipulate(plot(td_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "darkred" , lwd = 2 ,  
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(td_data)), region = slider(1,116))



# Take summaries of ROI's patients -----------------------------------------
list_summary <- function(list, fun){
  
  len <- length(list)
  tab <- matrix(NA, nrow = len, ncol = 116)

  for(i in 1:len){
    tab[i,] <- unname(sapply(list[[i]], fun))  }

  return(tab)
}





td_ROI_mean<- list_summary(td_scale, mean)
td_ROI_sd<- list_summary(td_scale, sd)



b <- list_summary(asd_scale, mean)

a <- cbind(a, rep(0, nrow(a)))
b <- cbind(b, rep(0, nrow(b)))







library(corrplot)
corrplot(cor(asd_data[[1]]), order="hclust")













# EXERCISE 3 --------------------------------------------------------------

rm(list = ls())

load("hw3_data.RData")
suppressWarnings()




# Load packages.
packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)


# Input.
k <- 5 # Dimensions
n0 <- 80  # Sample size 0
n1 <- 90  # Sample size 1
alpha <- .05 # significance level


# Reproducibility.
set.seed(234)


# Take random sample from a multi-normal distribution.
Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}


# Generate k random means.
generate_mean <- function(k){
  
  return(sample(1:10, k, replace = T))
}


# Generate random Covariance matrix (k*k).
generate_sigma <- function(k){
  sym_mat <- matrix(1, ncol = k, nrow = k, byrow = TRUE)
  diag_mat <- diag(rep(1.5, k))
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
  
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
    z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  kolm_t <<- kolm_t
  
  hist(kolm_t, main = "
       Kolmogorov-Smirnov statistic \n distribution under H_0",
       col = "skyblue", border = "white", breaks= 30)
  p_kolm <<- quantile(kolm_t , 1 - alpha)
  abline(v = p_kolm , col = "red" , lty = 3 , lwd = 2)
  box()
}

Friedman_procedure(P)


# Simulation to get info about ALPHA --------------------------------------

alpha_info <- function(P, percentile_kolm){
  prop_rej <- rep(NA, P)
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

  return(prop_rej)
}


data <- alpha_info(P = P, percentile_kolm = p_kolm)

barplot(prop.table(table(data)), col = c("red" , "blue"), main = "Proportion of times we accept-reject \n the null hypothesis  when is actually true \n using KS statistic", names.arg = c("Reject" , "Accept"), ylim = c(0,1))






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



power_info <- function(P, k, percentile_ks, increment){
  mu2 <- mu + increment 
  sigma2 <- sigma
  distance <- Wasserstein_distance(mu, mu2, sigma, sigma2)
  prop_ks <- rep(NA, P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    true_coef <- glm(label ~ ., data = u)$coefficients
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = true_coef)
    z_scores <- apply(z, MARGIN = 1, sigmoid, theta = true_coef)
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_ks[i] <- true_kolm > percentile_ks
  }
  data <- c(distance, mean(prop_ks))
  names(data) <- c('Distance', 'K_S')
  
  return(data)
}

k <- 5
mu <- 1:k
sigma <- generate_sigma(length(mu))
power_info(1000, k = length(mu), percentile_ks = p_kolm, increment = 10)





# Plot
n0 <- 200
n1 <- 200

P <- 1000
k <- 50
M <- seq(0, .5, .05)
mu <- generate_mean(k)
sigma <- generate_sigma(length(mu))

Wasserstein_distance(mu, mu+.5, sigma, sigma)

results <- matrix(NA, ncol = 2, nrow = 0)

for(l in M){
  a <- power_info(P, k = length(mu), percentile_ks = p_kolm, increment = l)
  results <- rbind(results, a)
}

results <- as.data.frame(results)

plot(results$Distance, results$K_S, type = 'l', main=paste0("Relation between Distance and Test's Power \n with k = ", k), lwd=3, col='skyblue', xlab='Distance', ylab="Test's Power", las=1, xlim=c(0,10))
abline(h=1, lty =3)


k20_n200 <- results







colors <- c('skyblue', 'lightgreen', 'orchid', 'darkorange')
par(mfrow=c(2,2))

### k = 2
plot(0:10, seq(0, 1, .1), type='n', ylim=c(0,1), main="Relation between Distance and Test's Power\n k = 2", xlab='Distance', ylab="Test's Power")
points(k2_n20, type='l', col=colors[1], lwd=2)
points(k2_n50, type='l', col=colors[2], lwd=2)
points(k2_n100, type='l', col=colors[3], lwd=2)
points(k2_n200, type='l', col=colors[4], lwd=2)
abline(h=1, lty =3)

legend('bottomright', legend = c(20, 50, 100, 200), col = colors, lwd=3, title='n0 & n1', bty = "n")




### k = 5
plot(0:10, seq(0, 1, .1), type='n', ylim=c(0,1), main="Relation between Distance and Test's Power\n k = 5", xlab='Distance', ylab="Test's Power")
points(k5_n20, type='l', col=colors[1], lwd=2)
points(k5_n50, type='l', col=colors[2], lwd=2)
points(k5_n100, type='l', col=colors[3], lwd=2)
points(k5_n200, type='l', col=colors[4], lwd=2)
abline(h=1, lty =3)

legend('bottomright', legend = c(20, 50, 100, 200), col = colors, lwd=3, title='n0 & n1', bty = "n")




### k = 10
plot(0:10, seq(0, 1, .1), type='n', ylim=c(0,1), main="Relation between Distance and Test's Power\n k = 10", xlab='Distance', ylab="Test's Power")
points(k10_n20, type='l', col=colors[1], lwd=2)
points(k10_n50, type='l', col=colors[2], lwd=2)
points(k10_n100, type='l', col=colors[3], lwd=2)
points(k10_n200, type='l', col=colors[4], lwd=2)
abline(h=1, lty =3)

legend('bottomright', legend = c(20, 50, 100, 200), col = colors, lwd=3, title='n0 & n1', bty = "n")



### k = 20
plot(0:10, seq(0, 1, .1), type='n', ylim=c(0,1), main="Relation between Distance and Test's Power\n k = 20", xlab='Distance', ylab="Test's Power")
points(k20_n20, type='l', col=colors[1], lwd=2)
points(k20_n50, type='l', col=colors[2], lwd=2)
points(k20_n100, type='l', col=colors[3], lwd=2)
points(k20_n200, type='l', col=colors[4], lwd=2)
abline(h=1, lty =3)

legend('bottomright', legend = c(20, 50, 100, 200), col = colors, lwd=3, title='n0 & n1', bty = "n")










packs <- c('MASS', 'pracma', 'psych','caTools','e1071')
lapply(packs, require, character.only = TRUE)



















k <- 5 # Dimensions
n0 <- 100  # Sample size 0
n1 <- 100  # Sample size 1
alpha <- .05 # significance level
set.seed(123)# reproducibility
P <- 1000
acc_rej_col <- c("#d3305d" , "#ABCDEF")


Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}



mu <- rep(0,k)
sigma <- diag(1 ,k)


Friedman_procedure <- function(P,x,z){
  kolm_t <- rep(NA, P)
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_model <- glm(label ~ ., data = u_p)
    x_scores <- predict(glm_model, x[,1:k])
    z_scores <- predict(glm_model, z_p[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  return(kolm_t)
  
}



P <- 1
alpha_info <- function(M){
  prop_rej <- rep(NA, M)
  p_values <- rep(NA , M)
  
  for(i in 1:M){
    x_p <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)  
    z_p <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u_p<- rbind(x_p, z_p) # Combine the data
    
    kk <- Friedman_procedure(P,x_p,z_p) # As mavi said
    glm_model <- glm(label ~ ., data = u_p)
    
    x_scores <<- predict(glm_model , x_p[,1:k])
    z_scores <<- predict(glm_model,  z_p[,1:k])
    
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_rej[i] <- true_kolm < quantile(kk , 1 - alpha)
    p_values[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$p.value
    
  }
  
  data = as.data.frame(cbind(prop_rej,p_values))
  colnames(data) <- c("KS","p")
  
  return(data)
}

data <- alpha_info(P = P)


hist(data$p)














# NEW TRy -----------------------------------------------------------------
packs <- c('MASS', 'pracma', 'psych','caTools','e1071', 'dplyr')
lapply(packs, require, character.only = TRUE)



k <- 5 # Dimensions
n0 <- 200  # Sample size 0
n1 <- 200  # Sample size 1
alpha <- .05 # significance level
set.seed(324)# reproducibility
P <- 100
acc_rej_col <- c("#d3305d" , "#ABCDEF")

mu <- rep(0,k)
sigma <- diag(1 ,k)


Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}




Fried <- function(P,x){
  
  kolm_t <- rep(NA, P)
  
  for(i in 1:P){
    z_l <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)
    u <- rbind(x, z_l)
    
    model <- glm(label ~ ., data = sample_n(u, ((n1 + n0) * .8)))
    x_scores <- predict(model, x[,1:k])
    z_scores <- predict(model, z_l[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  
  }
  
  return(kolm_t)

}

x <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)
a <- Fried(P, x)




P <- 15
alpha_info <- function(P){
  prop_rej <- rep(NA, P)
  p_values <- rep(NA , P)
  
  for(i in 1:P){
    
    x <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)  
    
    z_p <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u <- rbind(x, z_p) # Combine the data
    
    kk <- Fried(100,x) # As mavi said
    
    label0 <- rep(0,n0)
    label1 <- rep(1,n1)
    labels <- c(label0,label1)
    per_kolm <- rep(NA, P)

    idx <- sample(x = 1:(n0+n1), size = n0+n1)
    u$label <- labels[idx]
    glm_model <- glm(label ~ ., data = u ,maxit= 100)
    scores <- predict(glm_model , u[,1:k])
    p_values[i] <- ks.test(scores[u$label == 0],scores[u$label == 1])$p.value

    
    
    # glm_model <- glm(label ~ ., data = u_p, maxit = 10)
    # 
    # x_scores <- predict(glm_model , x[,1:k])
    # z_scores <- predict(glm_model,  z_p[,1:k])
    # 
    # test <- wilcox.test(x_scores, z_scores, alternative = "two.sided")
    # true_kolm <- test$statistic
    # 
    # prop_rej[i] <- true_kolm < quantile(kk , 1 - alpha)
    # p_values[i] <- test$p.value
    
  }
  
  data = as.data.frame(cbind(prop_rej,p_values))
  colnames(data) <- c("KS","p")
  
  return(data)
}

data <- alpha_info(P = P)

hist(data$p)
data$p
mean(data$p < .05)









Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}




M <- 1000

k <- 20
mu <- rep(0, k)
sigma <- diag(1, k)

n0 <- 40
n1 <- 50

x <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 1)

t <- rep(NA, M)
s <- rep(NA, M)

for(i in 1:M){
  z <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 0)
  
  
  u <- rbind(x, z)
  idx <- sample(x = 1:(n0+n1), size = n0+n1)
  u$label <- u$label[idx]
  
  model <- glm(label ~ ., data = u)
  
  x_scores <- predict(model, x)
  z_scores <- predict(model, z)
  
  
  t[i] <- ks.test(x_scores, z_scores)$p.value
  s[i] <- ks.test(x_scores, z_scores)$statistic

  }


mean(t < .05)
hist(t)

D <- 1.358 * sqrt((n0+n1)/(n0*n1))
mean(s > D)
D





# FINAL -------------------------------------------------------------------

Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}



k <- 5 # Dimensions
n0 <- 70  # Sample size 0
n1 <- 70  # Sample size 1
alpha <- .05 # significance level
set.seed(324)# reproducibility

M <- 100
P <- 100
acc_rej_col <- c("#d3305d" , "#ABCDEF")

mu <- rep(0,k)
sigma <- diag(1 ,k)


Friedman_procedure <- function(P,x,z){
  kolm_t <- rep(NA, P)
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label =1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_model <-glm(label ~ ., data = u_p)
    x_scores <- predict(glm_model , x[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  return(kolm_t)
  
}


alpha_info <- function(M){
  prop_rej <- rep(NA, M)
  
  for(i in 1:M){
    x_p <- Take_sample_normal(n = n0, mu = mu, sigma = sigma, label = 0)  
    z_p <- Take_sample_normal(n = n1, mu = mu, sigma = sigma, label = 1)  # same distributions
    u_p<- rbind(x_p, z_p) # Combine the data
    
    kk <- Friedman_procedure(P,x_p,z_p) # As mavi said
    glm_model <-glm(label ~ ., data = u_p)
    
    x_scores <- predict(glm_model , x_p[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    prop_rej[i] <- true_kolm < quantile(kk , 1 - alpha)
    
  }
  
  return(prop_rej)
}


prop_alpha <- alpha_info(M)

t <- proportions(table(prop_alpha))
barplot(t, col = acc_rej_col , main = "Proportion of times we accept-reject \n H0 when is actually true \n using KS statistic", ylim = c(0,1), names.arg = c("Reject" , "Accept"))




ks <- c(5,20,50,75,100)

alpha_k <- c()
for(k in ks){
  print(k)
  mu <- rep(0,k)
  sigma <- diag(1,k)
  alpha <- 1 - mean(alpha_info(M))
  alpha_k <- c(alpha_k, alpha)
  
}


plot(ks, alpha_k, type='b', col='gold3', pch=16)


