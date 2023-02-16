rm(list = ls())
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

# Scale the dataset
td_scale <- scale_datasets_list(td_data)
asd_scale <- scale_datasets_list(asd_data)

## Len of the dataset
sort( table( sapply(td_scale, function(x) nrow(x)) ), decreasing = T )
sort( table( sapply(asd_scale, function(x) nrow(x)) ), decreasing = T )

## Plot times series

library(manipulate)
manipulate(plot(asd_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "skyblue", lwd = 2, 
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(asd_data)), region = slider(1,116))


manipulate(plot(td_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "darkred" , lwd = 2 ,  
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(td_data)), region = slider(1,116))


####  get the variables
list_summary <- function(list, fun){
  
  len <- length(list)
  tab <- matrix(NA, nrow = len, ncol = 116)
  
  for(i in 1:len){
    tab[i,] <- unname(sapply(list[[i]], fun, na.rm = T))  }
  
  return(tab)
}

td_ROI_mean <- list_summary(td_scale, mean)
td_ROI_median <- list_summary(td_scale, median)
td_ROI_sd <- list_summary(td_scale, sd)

x_data <- cbind(td_ROI_mean,td_ROI_median,td_ROI_sd)
x_data <- cbind(x_data,rep(0,93)) 



asd_ROI_mean <- list_summary(asd_scale, mean)
asd_ROI_median <- list_summary(asd_scale, median)
asd_ROI_sd <- list_summary(asd_scale, sd)


z_data <- cbind(asd_ROI_mean,asd_ROI_median,asd_ROI_sd)
z_data <- as.data.frame(cbind(z_data,0)) 

names(z_data)[ncol(z_data)] <- 'label'



### Split the dataset
train_x_data <- z_data[1:60,]
train_x_data <- as.data.frame(train_x_data)
dim(train_x_data)
typeof(train_x_data)
test_x_data <- z_data[61:93,]


### Normality distribution
#[TODO] test the distributions normality

### Best parameter
mu_x_data <- apply(train_x_data[,1:k] , MARGIN = 2 , mean)
mu_x_data[is.na(mu_x_data)] <- mean(mu_x_data , na.rm = T)
sigma_x_data <- cov(train_x_data[,1:k])
sigma_x_data[is.na(sigma_x_data)] <- 0



### Friedmann procedure
k <- dim(train_x_data)[2] - 1
n0 <-  dim(train_x_data)[1]
n1 <- dim(z_data)[1]

###
Take_sample_normal <- function(n, mu, sigma, label){
  x <- mvrnorm(n, mu, sigma)
  x <- as.data.frame(cbind(x, label))
  colnames(x[length(colnames(x))]) <- "label"
  
  return(x)
}

Friedman_procedure <- function(P,xdata,zdata){
  kolm_t <- rep(NA, P)
  for(i in 1:P){
    z_p <- Take_sample_normal(n1, mu = mu_x_data, sigma = sigma_x_data , label =1) # Under H_0
    u_p <- rbind(xdata, z_p)
    glm_model <-glm(label ~ ., data = u_p)
    x_scores <- predict(glm_model , xdata[,1:k])
    z_scores <- predict(glm_model,  z_p[,1:k])
    kolm_t[i] <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
  }
  
  return(kolm_t)
  
}
###
d <- Friedman_procedure(P = 100 , xdata = train_x_data , zdata = z_data)

###