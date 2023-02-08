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
    scaled_list[[names[i]]] <-patient
  }
  return(scaled_list)
}

td_scale <- scale_datasets_list(td_data)
asd_scale <- scale_datasets_list(asd_data)



# Frequency of number observation  ----------------------------------------

# Count how many patients in the ASD group have time series of a given length
sort( table( sapply(td_scale, function(x) nrow(x)) ), decreasing = T )
sort( table( sapply(asd_scale, function(x) nrow(x)) ), decreasing = T )




# Plot of time series of each ROI foreach patient -------------------------

library(manipulate)

manipulate(plot(asd_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "skyblue", lwd = 2, 
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(asd_data)), region = slider(1,116))


manipulate(plot(td_data[[patient]][[region]], type='l', main = paste0("Patient-", patient, '\nROI-', region), col = "darkred" , lwd = 2 ,  
                xlab = 'observations', ylab = 'value'),
           patient = slider(1,length(td_data)), region = slider(1,116))



length(asd_data)

# Take summaries of ROI's patients -----------------------------------------

# a <- function(list, fun){
#   len <- length(list)
#   tab <- matrix(NA, nrow = 0, ncol = 116)
#   
#   for(i in 1:len){
#     print((sapply(list[i], mean)))
#     tab <- rbind(tab, unname(sapply(list[i], fun)))
#   }
#   
#   # return(tab)
# }
# 
# a(td_data, mean)
# 
# o <- unname(sapply(td_scale$caltech_0051478, mean))
