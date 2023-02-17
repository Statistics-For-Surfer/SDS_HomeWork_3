
packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)


Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}


Friedman_procedure <- function(P, x_data , z_data , permut = FALSE){
  x_fri <- x_data       # Copy the data
  z_p <- z_data
  kolm_t <- rep(NA, P)  # Pre-set the Kolmogorov-Smirnov statistic
  labels <- c(rep(0,n1),rep(1,n1))   #  Set of the label
  for(i in 1:P){ 
    if(permut == F){
      z_p <- Take_sample_normal(n1, mu, sigma, label =1)
    } # Under H_0
    idx <- sample(x = 1:(n0+n1), n0+n1)    # Shuffle the label
    
    x_fri$label <- labels[idx[1:n0]]                   # Permuted label
    z_p$label <- labels[idx[(n0+1):(n0+n1)]]           # Permuted label
    u_p <- as.data.frame(rbind(x_fri,z_p))             # Row bind the two data frame
    glm_f <- glm(label~., data = u_p)                  # Train the model
    scores <- predict(glm_f ,u_p[,1:k])                # Compute the scores
    
    kolm_t[i] <- ks.test(scores[u_p$label == 0] , scores[u_p$label == 1])$statistic  # Save the i -th statistics value
    
  }
  return(kolm_t)
  
}


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
  
  dist <- Wasserstein_distance(mu, mu2, sigma, sigma2)
  
  prop_ks <- rep(NA, M)
  for(i in 1:M){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    p_model <- glm(label ~ ., data = u)
    x_scores <- predict(p_model, x)
    z_scores <- predict(p_model, z) 
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    kk <- Friedman_procedure(P,x_data = x, z_data = z)
    prop_ks[i] <- true_kolm > quantile(kk, 1-alpha)
  }
  
  return(c(dist, mean(prop_ks)))
}






k <- 100
mu <- rep(0, k)
sigma <- diag(1, k)

Wasserstein_distance(mu, mu + .3, sigma, sigma)

grid <- c(0, .07, .1, .14, .17, .2, .3)
grid <- seq(0, .3, .05)



# Parameters
P <- 100
M <- 100
n0 <- 80  
n1 <- 90  
alpha <- .05 

ns <- c(20, 40, 60, 80, 100)
ks <- c(5, 20, 50, 100)
colors <- c('cornflowerblue', 'chartreuse3', 'darkorchid2', 'chocolate1')



# k = 5
k <- 5
mu <- rep(0, k)
sigma <- diag(1, k)
print(paste0("k = ", k , " -->"))

grid <- c(0, .2, .3, .4, .5, 1, 2)
distance_k5 <- matrix(NA, 0, 2)

for(i in grid){
  print(i)
  power_k <- power_info(M, P, n0, n1, mu, mu+i, sigma, sigma)
  distance_k5 <- rbind(distance_k5, power_k)
}

# k = 20
k <- 20
mu <- rep(0, k)
sigma <- diag(1, k)
print(paste0("k = ", k , " -->"))

grid <- c(0, .1, .15, .2, .25, .3, .6)
distance_k20 <- matrix(NA, 0, 2)

for(i in grid){
  print(i)
  power_k <- power_info(M, P, n0, n1, mu, mu+i, sigma, sigma)
  distance_k20 <- rbind(distance_k20, power_k)
}

# k = 50
k <- 50
mu <- rep(0, k)
sigma <- diag(1, k)
print(paste0("k = ", k , " -->"))

grid <- c(0, .07, .1, .13, .15, .17, .4)
distance_k50 <- matrix(NA, 0, 2)

for(i in grid){
  print(i)
  power_k <- power_info(M, P, n0, n1, mu, mu+i, sigma, sigma)
  distance_k50 <- rbind(distance_k50, power_k)
}

# k = 100
k <- 100
mu <- rep(0, k)
sigma <- diag(1, k)
print(paste0("k = ", k , " -->"))

grid <- c(0, .07, .1, .14, .17, .2, .3)
distance_k <- matrix(NA, 0, 2)

for(i in grid){
  print(i)
  power_k1 <- power_info(M, P, n0, n1, mu, mu+i, sigma, sigma)
  distance_k100 <- rbind(distance_k100, power_k)
}





# Final plot
plot(distance_k5, xlim=c(0,6),  xlab='distance', ylab='power', las=1, type='l', lwd=3, col=colors[1])
points(distance_k20, type='l', lwd=3, col=colors[2])
points(distance_k50, type='l', lwd=3, col=colors[3])
points(distance_k100, type='l', lwd=3, col=colors[4])

legend('bottomright', legend=ks, col=colors, lwd=3, title='k')

grid()
