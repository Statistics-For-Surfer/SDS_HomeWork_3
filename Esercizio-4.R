rm(list = ls())






load('hw3_data.RData')
list_summary <- function(list, fun){
  len <- length(list)
  tab <- matrix(NA, nrow = len, ncol = 116)
  
  for(i in 1:len){
    tab[i,] <- unname(sapply(list[[i]], fun, na.rm = T))  }
  
  return(tab)
}

td_ROI_mean <- list_summary(td_data, mean)
td_ROI_median <- list_summary(td_data, median)
td_ROI_sd <- list_summary(td_data, sd)
x_data <- c()

x_data <- as.data.frame(cbind(td_ROI_mean, td_ROI_median, td_ROI_sd))

x_data <- as.data.frame(cbind(x_data,rep(0, length(td_data))))




names(x_data)[ncol(x_data)] <- 'label'
#####
asd_ROI_mean <- list_summary(asd_data, mean)
asd_ROI_median <- list_summary(asd_data, median)
asd_ROI_sd <- list_summary(asd_data, sd)
z_data <- c()
z_data <- as.data.frame(cbind(asd_ROI_mean, asd_ROI_median, asd_ROI_sd))
z_data <- as.data.frame(cbind(z_data, rep(1, length(asd_data))))
names(z_data)[ncol(z_data)] <- 'label'


#####



u <- as.matrix(rbind(x_data,z_data))
k <- dim(x_data)[2]

mod_1 <- glmnet(
  x = u[,1:(k-1)],
  y = u[,k],
  alpha = 0.2,
  maxit = 10,
  family = "gaussian",
)

scores <- predict(mod_1,u[,1:(k-1)])

hist(scores[1:93] ,  xlim = c(0,1))
hist(scores[94:178],add = T)
ks.test(scores[1:93] , scores[94:178])$statistic

#####
x_data$label <- c(rep(0,50),rep(1,43))
x_data <- as.data.frame(x_data)
mod <- glm(label~., data = x_data)
scores <- predict(mod,x_data[,1:(k-1)])

scores_2 <- predict(mod,x_data[,1:(k-1)] , s = .005)
hist(scores_2[1:50],xlim = c(0,1))
hist(scores_2[51:93],xlim = c(0,1), add  = T)
