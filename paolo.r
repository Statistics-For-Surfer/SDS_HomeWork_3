rm(list = ls())
load("hw3_data.RData")



##### Esercise 3


# Case 1 (The two distributions are the same)
 

n0 <- 80  # numero campioni con etichetta 0
n1 <- 90  # sample size of labels 1
k <- 1
# Get a sample from a specific distribution (es: normale)

mu = 10
sigma = 2

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
  p_negative <- dnorm(n0, 0.8, 0.2 )
  p_positive <- dnorm(n1 , 0.5 , 0.4)
  # Perform the Kolmogorov , Chi squared , Mann Test
  
  t_hat <- c(t , t.test(p_negative , p_positive , alternative = "two.sided"))
  
}

# We end uo with a set of values for t_hat (distribution of t_hat under H0)

### Perform the test

t_hat_true <- t.test(x, z , alternative = "two.sided")

if( t_hat_true > qunatile(1-alpha , t_hat ) ){decisio  <- "reject H0"}











t <- c()
for (i in 1:P){
  z <-rnorm(m , mu , sigma)
  
  
} 


length(u)
?glm

glm(u ~  x + z )