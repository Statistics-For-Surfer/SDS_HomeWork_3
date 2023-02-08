<<<<<<< HEAD
rm(list=ls())

=======
  rm(list = ls())
>>>>>>> 3e2b32eddbf66a85eada3eea8060a3b480e170e3
load("hw3_data.RData")




# Scale dataset -----------------------------------------------------------
scale_datasets_list <- function(ls){
  scaled_list <- list()
  names <- names(ls)
  for(i in 1:length(ls)){
    
    #### Scaling by column.
    patient <- data.frame(apply(ls[[i]], 2, scale))
    colnames(patient) <- colnames(ls[[i]])
    scaled_list[[names[i]]] <-patient
  }
  return(scaled_list)
}

td_scale <- scale_datasets_list(td_data)
asd_scale <- scale_datasets_list(asd_data)



# Frequency of number observation  ----------------------------------------

# Count how many patients in the ASD group have time series of a given length
<<<<<<< HEAD
sort( table( sapply(td_scale, function(x) nrow(x)) ), decreasing = T )
sort( table( sapply(asd_scale, function(x) nrow(x)) ), decreasing = T )




# Plot of time series of each ROI foreach patient -------------------------
=======
  sort( table( sapply(asd_data, function(x) nrow(x)) ), decreasing = T )

>>>>>>> 3e2b32eddbf66a85eada3eea8060a3b480e170e3
library(manipulate)

manipulate(plot(asd_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "skyblue", lwd = 2, 
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(asd_data)), region = slider(1,116))


manipulate(plot(td_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "darkred" , lwd = 2 ,  
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(td_data)), region = slider(1,116))

<<<<<<< HEAD


length(asd_data)

# Take summaries of ROI's patients -----------------------------------------

a <- function(list, fun){
  # if(fun == 'mean') {fun <- mean}
  # else if(fun == 'sd') {fun <- sd}
  # else if(fun == 'median') {fun <- median}
  # fun <- fun
  tab <- matrix(NA, nrow = length(list), ncol = 116)
  
  len <- length(list)
  
  for(i in 1:len){
    print(lapply(list[i], mean))
    # print(unname(unlist(lapply(list[i], fun))))
    tab[i,] <- unname(unlist(lapply(list[i], fun)))
  }
  
  # return(tab)
}

a(td_scale, 'mean')

o <- unname(unlist(lapply(td_scale$caltech_0051478, mean)))