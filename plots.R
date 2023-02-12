
rm(list = ls())





# Load packages.
packs <- c('MASS', 'pracma', 'psych')
lapply(packs, require, character.only = TRUE)

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
  
  return(sample(1:10, k, replace = T))
}


# Generate random Covariance matrix (k*k).
generate_sigma <- function(k){
  sym_mat <- matrix(1, ncol = k, nrow = k, byrow = TRUE)
  diag_mat <- diag(rep(.5, k))
  sigma <- diag_mat + sym_mat
  
  return(sigma)
}


# Sigmoid Function.
sigmoid <- function(x, theta){
  n_features <- length(theta)
  n <- exp(theta[1] + sum(x[-length(x)]*theta[2:n_features]))
  
  return(n / (1+n))
}

# info about alpha.
alpha_info <- function(P, percentile_kolm, mu, sigma, n0, n1){
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
  
  return(1-mean(prop_rej))
}



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






# Alpha on k fixed n0 & n1 ------------------------------------------------

alpha <- .05
n0 <- 70
n1 <- 50
P <- 1000
k <- c(5, 20, 50, 100)

alpha_plot <- function(n0, n1, P, k){
  k_alpha <- matrix(NA, nrow=0, ncol=2)
  
  for(el in k){
  mu <- generate_mean(el)
  sigma <- generate_sigma(length(mu))
  
  x <- Take_sample_normal(n = n0, mu=mu, sigma = sigma, label = 0)
  
  kolm_t <- rep(NA, P)
  
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu, sigma, label = 1) # Under H_0
    u_p <- rbind(x, z_p)
    glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
    x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
    z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  p_kolm <- quantile(kolm_t, 1 - alpha)
  
  alpha_hat <- alpha_info(P, p_kolm, mu, sigma, n0, n1)
  k_alpha <- rbind(k_alpha, c(el, alpha_hat))}
  return(k_alpha)
}

k_alpha_val <- alpha_plot(n0, n1, P, k)

plot(k_alpha_val, xaxt = 'n', xlab='k', ylab=expression(hat(alpha)), type = 'b', main = paste0('alpha estimate for different k\n with n0 = ', n0, ' and n1 = ', n1), pch=18, col='gold3', lwd=2)
axis(1, at=k, labels=k)


 



# Alpha on n1 for fixed n0 & k --------------------------------------------

alpha <- .05
n0 <- 50
n1 <- seq(10, 100, 10)
P <- 1000
k <- 20

alpha_plot <- function(n0, n1, P, k){
  k_alpha <- matrix(NA, nrow=0, ncol=2)
  
  for(el in n1){
    mu <- generate_mean(k)
    sigma <- generate_sigma(length(mu))
    
    x <- Take_sample_normal(n = n0, mu=mu, sigma = sigma, label = 0)
    
    kolm_t <- rep(NA, P)
    
    for(i in 1:P){
      z_p <- Take_sample_normal(n = el, mu=mu, sigma = sigma, label = 1)
      u_p <- rbind(x, z_p)
      glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
      x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
      z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
      kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    }
    p_kolm <- quantile(kolm_t, 1 - alpha)
    
    alpha_hat <- alpha_info(P, p_kolm, mu, sigma, n0, el)
    k_alpha <- rbind(k_alpha, c(el, alpha_hat))}
    return(k_alpha)
}

k_alpha_val <- alpha_plot(n0, n1, P, k)

plot(k_alpha_val, xaxt = 'n', xlab='n1', ylab=expression(hat(alpha)), type = 'b', main = paste0('alpha estimate for different n1\n with n0 = ', n0, ' and k = ', k), pch=18, col='gold3', lwd=2)
axis(1, at=n1, labels=n1)








# Alpha on n0 for fixed n1 & k --------------------------------------------

alpha <- .05
n1 <- 50
n0 <- seq(10, 100, 10)
P <- 1000
k <- 20

alpha_plot <- function(n0, n1, P, k){
  k_alpha <- matrix(NA, nrow=0, ncol=2)
  
  for(el in n0){
    mu <- generate_mean(k)
    sigma <- generate_sigma(length(mu))
    
    x <- Take_sample_normal(n = el, mu=mu, sigma = sigma, label = 0)
    
    kolm_t <- rep(NA, P)
    
    for(i in 1:P){
      z_p <- Take_sample_normal(n = n1, mu=mu, sigma = sigma, label = 1)
      u_p <- rbind(x, z_p)
      glm_coef <- unname(glm(label ~ ., data = u_p)$coefficients)
      x_scores <- apply(x, MARGIN = 1, sigmoid, theta = glm_coef)
      z_scores <- apply(z_p, MARGIN = 1, sigmoid, theta = glm_coef) 
      kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    }
    p_kolm <- quantile(kolm_t, 1 - alpha)
    
    alpha_hat <- alpha_info(P, p_kolm, mu, sigma, el, n1)
    k_alpha <- rbind(k_alpha, c(el, alpha_hat))}
  return(k_alpha)
}

k_alpha_val <- alpha_plot(n0, n1, P, k)

plot(k_alpha_val, xaxt = 'n', xlab='n0', ylab=expression(hat(alpha)), type = 'b', main = paste0('alpha estimate for different n0\n with n1 = ', n1, ' and k = ', k), pch=18, col='gold3', lwd=2)
axis(1, at=n0, labels=n0)














M <- seq(0, 1, .1)
Wasserstein_distance(mu, mu+1, sigma, sigma)

results <- matrix(NA, ncol = 2, nrow = 0)

for(l in M){
  a <- power_info(P, k = length(mu), percentile_ks = p_kolm, increment = l)
  results <- rbind(results, a)
}

results <- as.data.frame(results)

plot(results$Distance, results$K_S, type = 'l', main=paste0("Relation between Distance and Test's Power \n with k = ", k), lwd=3, col='skyblue', xlab='Distance', ylab="Test's Power", las=1, xlim=c(0,10), ylim = c(0,1))


















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



















