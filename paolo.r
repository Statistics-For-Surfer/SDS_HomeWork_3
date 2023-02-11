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

mu <- c(1,2,3,4,5)
k <- 5

# Create a diagonal matrix with positive elements
diag_mat <- diag(1:k)

# Create a symmetric matrix with non-negative elements
sym_mat <- matrix(c(0, 1, 1, 1, 1,
                    1, 0, 1, 1, 1,
                    1, 1, 0, 1, 1,
                    1, 1, 1, 0, 1,
                    1, 1, 1, 1, 0), nrow = k, ncol = k, byrow = TRUE)

# Sum the diagonal and symmetric matrices
sigma <- diag_mat + sym_mat



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


info_power <- function(P, percentile_kolm, percentile_mann, space_theta){
  n <- length(space_theta)
  l <- 0
  prop_rej <- matrix(NA , ncol = 3 , nrow = (n*(n-1)) * P )
  for (theta in space_theta){
    for(theta_2 in space_theta){
      if (theta_2 == theta){
        next
      }
      l <- l + 1
      for(i in 1:P){
        x <- Take_sample_normal(k= k,n = n0,mu= theta , sigma = sigma,label = 0)
        z <- Take_sample_normal(k= k,n = n0,mu= theta_2 , sigma = sigma,label = 1)
        #[TODO] TAKE the distance between the distribitions
        u <- rbind(x,z) # Combine the data
        true_coef <- glm(label ~. , data = u)$coefficients
        x_scores <- apply(x ,MARGIN = 1 , sigmoid ,theta = true_coef)
        z_scores <- apply(z ,MARGIN = 1 , sigmoid ,theta = true_coef)
        true_kolm <- ks.test(x_scores,z_scores,alternative = "two.sided")$statistic
        true_mann <- wilcox.test(x_scores,z_scores, alternative = "two.sided")$statistic
        # Perform the Test
        prop_rej[l*i,1] <- true_kolm < percentile_kolm
        prop_rej[l*i,2] <- true_mann < percentile_mann
        prop_rej[l*i,3] <-  "d"#distance
        
      } }}
  data = as.data.frame(prop_rej)
  return(data)
}



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



#######################################
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


