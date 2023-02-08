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



# Simulation 3
rm(list=ls())

A <- matrix(rnorm(1000, 10, 3), 50, 20)
B <- matrix(rnorm(1000, 12, 3), 50, 20)

tab_A <- data.frame(cbind(A, rep(0, 50)))
tab_B <- data.frame(cbind(B, rep(1, 50)))

coeff_A <- unname(glm(X21~., data = tab_A, family = 'binomial')$coefficients)
coeff_B <- glm(X21~., data = tab_B, family = 'binomial')$coefficients


A <- cbind(rep(1, 50), A)
B <- cbind(rep(1, 50), B)


sig_A <- function(x){
  return(1 / (1 + exp(- x %*% coeff_A)))
}


sig_B <- function(x){
  return(1 / (1 + exp(- x %*% coeff_B)))
}


apply(A, 1, sig_A)
apply(B, 1, sig_B)




