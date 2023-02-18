load('hw3_data.RData')



# Put all data together.
super_data <- c()
for(el in 1:length(asd_data))  super_data <- c(super_data, asd_data[el])
for(el in 1:length(td_data))  super_data <- c(super_data, td_data[el])


# Take unique names of labs.
lab <- names(super_data)
init_lab <- c()
for(el in lab) init_lab <- c(init_lab, substr(el,1,2))
u_init_lab <- unique(init_lab)


# Scale every data frame by lab.
for(l in u_init_lab){
  pos <- which(grepl(l, init_lab))
  d <- matrix(NA, 0, 116)
  for(i in pos) d <- rbind(d, super_data[[i]])
  d <- as.matrix(d)
  m <- mean(d)
  s <- sd(d)
  for(i in pos) super_data[[i]] <- (super_data[[i]] - m) / s
}

# Get new scaled data frames.
asd_data_lab_scale <- super_data[1:85]
td_data_lab_scale <- super_data[86:178]




# Get summaries from all patient's ROIs.
list_summary <- function(list, fun){
  
  len <- length(list)
  tab <- matrix(NA, nrow = len, ncol = 116)
  
  for(i in 1:len){
    tab[i,] <- unname(sapply(list[[i]], fun, na.rm = T))  }
  
  return(tab)
}

td_ROI_mean <- list_summary(td_data_lab_scale, mean)


