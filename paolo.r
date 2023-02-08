rm(list = ls())
load("hw3_data.RData")



##### Esercise 3


# Get a sample from a specific distribution
n <- 80
m <- 90
mu = 10
sigma = 2


x <- rnorm(n , mu , sigma)
y_x <- rep(1, n)
y_z <- rep(0, m)
u <- c(y_x ,y_z)

P = 1000

t <- c()
for (i in 1:P){
  z <-rnorm(m , mu , sigma)
  
  
} 


length(u)
?glm

glm(u ~  x + z )