rm(list = ls())
load("hw3_data.RData")

# Count how many patients in the ASD group have time series of a given length
sort( table( sapply(asd_data, function(x) nrow(x)) ), decreasing = T )

library(manipulate)

manipulate(plot(asd_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "skyblue", lwd = 2, 
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(asd_data)), region = slider(1,116))


manipulate(plot(td_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "darkred" , lwd = 2 ,  
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(td_data)), region = slider(1,116))

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


