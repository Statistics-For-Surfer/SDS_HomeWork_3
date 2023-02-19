# NEW procedure -----------------------------------------------------------


packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)
colors <- c('cornflowerblue', 'chartreuse3', 'darkorchid2', 'chocolate1')




Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)    # Random generate from a multivariate normal 
  x <- as.data.frame(cbind(x, label))   # add label
  colnames(x[length(colnames(x))]) <- "label"   # col label
  
  return(x)
}


Friedman_procedure <- function(P,x_data , z_data , permut = FALSE){
  x_fri <- x_data       # Copy the data
  z_p <- z_data         # copy the z data
  kolm_t <- rep(NA, P)  # Pre-set the Kolmogorov-Smirnov statistic
  labels <- c(rep(0,n0),rep(1,n1))   # Set of the label
  for(i in 1:P){                     # Loop 
    if(permut == F){                 # Check goodnees of fit
      z_p <- Take_sample_normal(n1, mu, sigma, label =1)  # Sample under the null F
    } 
    if(permut == T){                 # Two sample test
      idx <- sample(x = 1:(n0+n1), n0+n1)    # Shuffle the label
      
      x_fri$label <- labels[idx[1:n0]]            # Permuted label
      z_p$label <- labels[idx[(n0+1):(n0+n1)]]}   # Permuted label
    u_p <- as.data.frame(rbind(x_fri,z_p))      # Row bind the two data frame
    glm_f <- glm(label~., data = u_p)           # Train the model
    scores <- predict(glm_f ,u_p[,1:k])         # Compute the scores
    
    kolm_t[i] <- ks.test(scores[u_p$label == 0] , scores[u_p$label == 1])$statistic  
  }
  return(kolm_t)
}


beta_q <- function(sigma1, sigma2){
  trace_1 <- tr(sigma)    # Trace of sigma X
  trace_2 <- tr(sigma2)   # Trace of sigma Y
  trace_3 <- tr(sqrtm(sqrtm(sigma1)$B %*% sigma2 %*% sqrtm(sigma1)$B)$B)
  
  return(trace_1 + trace_2 - 2*trace_3)
}

Wasserstein_distance <- function(mu1, mu2, sigma1, sigma2){
  norma <- norm(mu1-mu2, type = "2")**2   
  beta_quadro <- beta_q(sigma1, sigma2)
  
  return(norma + beta_quadro)
}


power_info <- function(P, k, increment , permut = F){
  
  prop_ks <- rep(NA, P)    # Pre-set proportions of accept/reject
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)  # Take sample of x
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1) # Take sample of z (From a different distribution)
    u <- rbind(x,z)    # Combine the dataset
    p_model <- glm(label ~ ., data = u)   # Perform the obs model
    x_scores <- predict(p_model, x)       # Predict the obs scores
    z_scores <- predict(p_model, z)       # ""
    kolm_obs <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic                    # Obs statistic
    if(permut == FALSE){
      kk <- Friedman_procedure(P,x_data = x, z_data = z)
    }
    if(permut == TRUE){
      kk <- Friedman_procedure(P,x_data = x, z_data = z , permut = T)
    }
    
    prop_ks[i] <- kolm_obs > quantile(kk,1-alpha)
  }
  data <-  mean(prop_ks)
  
  
  return(data)
}




k <- 20
mu <- rep(0, k)
sigma <- diag(1, k)

Wasserstein_distance(mu, mu + .34, sigma, sigma)

grid <- c(0, .07, 1,.13, .15, .18, .23)


i <- .34
mu2 <- mu + i
sigma2 <- sigma
distance <- Wasserstein_distance(mu, mu2, sigma, sigma2)
powers <- power_info(P, k = length(mu),increment = i , permut = T)

a <- c(distance, powers)
b <- p_k20[7, ]
p_k20 <- p_k20[-c(7,8),]
p_k20 <- rbind(p_k20, a,b)


P <- 100
n0 <- 70
n1 <- 60
alpha <- .05




# k = 5
k <- 5
mu <- rep(0,k)
sigma <- diag(1,k)
grid <- c(0, .2, .3, .4, .5, 1, 2)
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_k5 <- cbind(distance, powers)



# k = 20
k <- 20
mu <- rep(0,k)
sigma <- diag(1,k)
grid <- c(0, .1, .15, .2, .25, .3, .6)
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_k20 <- cbind(distance, powers)

# k = 50
k <- 50
mu <- rep(0,k)
sigma <- diag(1,k)
grid <- c(0, .07, .1, .15, .18, .23, .4)
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_k50 <- cbind(distance, powers)




# k = 100
k <- 100
mu <- rep(0,k)
sigma <- diag(1,k)
grid <- c(0, .07, .1, .14, .17, .2, .3, .4)
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_k100 <- cbind(distance, powers)




tipo <- 'b'

plot(p_k5, type=tipo, col=colors[1], lwd=2, ylim=c(0,1), xlim=c(0, 7))
points(p_k20,  type=tipo, col=colors[2], lwd=2)
points(p_k50,  type=tipo, col=colors[3], lwd=2)
points(p_k100,  type=tipo, col=colors[4], lwd=2)



distance_power_k <- list(p_k5, p_k20, p_k50, p_k100)
save(distance_power_k, file='distance_power_k.RData')




load('distance_power_k.RData')







# n0 ----------------------------------------------------------------------

P <- 80
n0s <- c(20,50,70,100) 

n1 <- 70
k <- 20
mu <- rep(0,k)
sigma <- diag(1,k)
grid <- c(0, .1, .15, .2, .25, .3, .35, .4)


n0 <- 20
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_n0_20 <- cbind(distance, powers)




n0 <- 50
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_n0_50 <- cbind(distance, powers)



n0 <- 70
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_n0_70 <- cbind(distance, powers)





n0 <- 100
powers <- c()
distance <- c()
for (i in (grid)){
  print(i)
  mu2 <- mu + i
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers ,power_info(P, k = length(mu),increment = i , permut = T))
}

p_n0_100 <- cbind(distance, powers)


tipo='b'
plot(p_n0_20, type=tipo, col=colors[1])
points(p_n0_50, type=tipo, col=colors[2])
points(p_n0_70, type=tipo, col=colors[3])
points(p_n0_100, type=tipo, col=colors[4])


distance_power_n0 <- list(p_n0_20, p_n0_50, p_n0_70, p_n0_100)
save(distance_power_n0, file = 'distance_power_n0.RData')
