rm(list = ls())
load("hw3_data.RData")
suppressWarnings()




##### Esercise 3
set.seed(123) # Reproducibility

# Case 1 (The two distributions are the same)
n0 <- 80  # numero campioni con etichetta 0
n1 <- 90  # sample size of labels 1
k <- 5
# Get a sample from a specific distribution (es: normale)
mu = 10
sigma = 2


take_data_distribution <- function(k,n, label = NULL){
  x <- matrix(NA,n,k)
  for(i in 1:k){
    x[,i] <-rnorm(n, mean = mu , sd = sigma)
  }

  x <- as.data.frame(cbind(x,label))
  colnames(x[length(colnames(x))]) <-"label" 
  return(x)
}




take_data_distribution_gamma <- function(k,n, label = NULL){
  x <- matrix(NA,n,k)
  for(i in 1:k){
    x[,i] <-rgamma(n, shape = 8, rate = 5)
  }
  
  x <- as.data.frame(cbind(x,label))
  colnames(x[length(colnames(x))]) <-"label" 
  return(x)
}




x <- take_data_distribution(k,n0,label = 0)
z <- take_data_distribution(k,n1,label = 1)


## Combine the two dataset

u <- rbind(x,z) # Actual data

#### Implement the Friedman procedure
P <- 1000



sigmoid <- function(x, theta){
  n_features <- length(theta)
  
  n <- exp(theta[1] + sum(x*theta[2:n_features]))
  return(n / (1+n))}

kolm_t <- rep(NA , P)
mann_t <- rep(NA, P)

for(i in 1:P){
  z_p <- take_data_distribution(k,n1,label =1)
  u_p <- rbind(x,z_p)
  glm_coef <- glm(label ~. , data = u_p)$coefficients
  x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = glm_coef)
  z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = glm_coef) #[TODO] See if make sense
  kolm_t[i] <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
  mann_t[i] <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
}

par(mfrow = c(1,2))
hist(kolm_t)
hist(mann_t)
par(mfrow = c(1,1))

########
# Use the actual data
alpha <- .05
# Simulation to get info about alpha
for(i in 1:P){
  x <- take_data_distribution(k,n0,label = 0)
  z <- take_data_distribution(k,n1,label = 1)
  u <- rbind(x,z) # Actual data
  true_coef <- glm(label ~. , data = u)$coefficients
  x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = true_coef)
  z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = true_coef)
  true_kolm <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
  true_mann <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
  abline(v= true_kolm , col = "blue")
  

}

# Simulation to get info about power
for(i in 1:P){
  x <- take_data_distribution(k,n0,label = 0)
  z <- take_data_distribution_gamma(k,n1,label = 1)
  u <- rbind(x,z) # Actual data
  true_coef <- glm(label ~. , data = u)$coefficients
  x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = true_coef)
  z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = true_coef)
  true_kolm <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
  true_mann <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
  abline(v= true_kolm , col = "blue")
  
  
}








abline(v = quantile(kolm_t, 1 - alpha), col = "red", lwd = 3)
abline(v= true_kolm , col = "blue")

hist(mann_t)
abline(v = quantile(mann_t, 1 - alpha), col = "red")
abline(v= true_mann , col = "blue")
########
# Simulation to get the power of the test




########




# Esiste un pacchetto per la normale multivariate
x <- rnorm(n0, mu, sigma)     # [TODO] Multivariate case
y_x <- rep(0, n0)

z <- rnorm(n1, mu ,sigma)
y_z <- rep(1, n1)



u <- c(y_x ,y_z)
length(u)

#####  Fisher's permutations
?t.test
P <-  1000
?sample
for(i in 1:P){
  
  u_P <- sample(u , size = length(u))    # We hypotized that the 2 distributions are the same
  glm(u_P ~ x  + z)   # sotto forma di dataset   F(u)
  ### Compute the score of the 2 distributions (positive and negative) 
  p_negative <- dnorm(n0, 0.8, 0.2 )  # F(x)
  p_positive <- dnorm(n1 , 0.5 , 0.4) # F(z)
  # Perform the Kolmogorov , Chi squared , Mann Test
  
  t_hat <- c(t , t.test(p_negative , p_positive , alternative = "two.sided"))
  
}

# We end uo with a set of values for t_hat (distribution of t_hat under H0)

### Perform the test

t_hat_true <- t.test(x, z , alternative = "two.sided")

if( t_hat_true > qunatile(1-alpha , t_hat ) ){decisio  <- "reject H0"}



####### Friedman's Procedure


n0 <- 80  # numero campioni con etichetta 0
n1 <- 90  # sample size of labels 1
k <- 1
# Get a sample from a specific distribution (es: normale)

mu = 10
sigma = 2


# Esiste un pacchetto per la normale multivariate
x <- rnorm(n0, mu, sigma)     # [TODO] Multivariate case
y_x <- rep(0, n0)

z <- rnorm(n1, mu ,sigma)
y_z <- rep(1, n1)



u <- c(y_x ,y_z)
length(u)


?t.test
P <-  1000
for(i in 1:P){
  z_P <- rnorm(n1 , mu + 2 , sigma/2) # Hypotized a distribution for z_s sample
  glm(u_P ~ x + z_P)
  #### Compute the the score for positive and negative lables
  p_negative <- F(x)
  p_positive <- F(z)  # [TODO] check if z or z_P
  # Perform the test (Kolmogorov, chi-squared , Mann)
  t_hat_P <- c(t_hat_P , t.test(p_negative, p_positive , "two.sided"))
}

# We end uo with a set of values for t_hat (distribution of t_hat under H0)

### Perform the test

t_hat_true <- t.test(x, z , alternative = "two.sided")

if( t_hat_true > qunatile(1-alpha , t_hat ) ){decisio  <- "reject H0"}





data("iris")
i <- iris[iris$Species == "setosa" | iris$Species == "versicolor" ,]
m <- glm(Species ~ Petal.Length , data = i , family = binomial)

m$coefficients[2:length(m$coefficients)]






#curve(sigmoid(x,m$coefficients),from = 2 , to = 3)


X <- i$Petal.Length[i$Species == "setosa"]
X
Z <- i$Petal.Length[i$Species == "versicolor"]


Z_p <- rnorm(50 , mean = 4.2, sd = 1)
data <- data.frame(c(X,Z_p), i$Species)
data



model <- glm(i.Species ~ . , data= data,family = binomial)
model$coefficients

a <- sapply(X,sigmoid,theta = model$coefficients)
b <- sapply(Z,sigmoid,theta = model$coefficients)

t.test(a,b)

scores_meno <- apply(X ,MARGIN = 1 , sigmoid ,theta = model$coefficients)

scores_piu <- apply(Z ,MARGIN = 1, sigmoid ,theta = m$coefficients)


t.test(scores_meno,scores_piu )


