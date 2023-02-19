
rm(list = ls())
load('hw3_data.RData')

Friedman_procedure <- function(P,x_data , z_data , permut = FALSE){
  x_fri <- x_data       # Copy the data
  z_p <- z_data         # copy the z data
  kolm_t <- rep(NA, P)  # Pre-set the Kolmogorov-Smirnov statistic
  labels <- c(rep(0,n0),rep(1,n1))   # Set of the label
  for(i in 1:P){                     # Loop 
    if(permut == F){                 # Check goodnees of fit
      z_p <- Take_sample_normal(n1, mu, sigma, label =1)  # Sample under the null F
    } 
    if(permut == T){                 # Two sample test
      idx <- sample(x = 1:(n0+n1), n0+n1)    # Shuffle the label
      
      x_fri$label <- labels[idx[1:n0]]            # Permuted label
      z_p$label <- labels[idx[(n0+1):(n0+n1)]]
    }   # Permuted label
    u_p <- as.matrix(rbind(x_fri,z_p))      # Row bind the two data frame
    glm_f <-glmnet(
      x = u_p[,1:(k-1)],
      y = u_p[,k],
      family = "gaussian",
      alpha = 0.2,
      maxit = 100
    )
    
    # Train the model
    scores <- predict(glm_f ,u_p[,1:(k-1)], s = 0.05)         # Compute the scores
    u_p <- as.data.frame(u_p)
    kolm_t[i] <- ks.test(scores[u_p$label == 0] , scores[u_p$label == 1])$statistic  
  }
  return(kolm_t)
  
}



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





# Function to get summaries from all patient's ROIs.
list_summary <- function(list, fun){
  len <- length(list)
  tab <- matrix(NA, nrow = len, ncol = 116)
  
  for(i in 1:len){
    tab[i,] <- unname(sapply(list[[i]], fun, na.rm = T))  }
  
  return(tab)
}


asd_corr <- asd_data_lab_scale
for(p in 1:length(asd_corr)){
  patient <- asd_corr[[p]]
  for(i in 1:ncol(patient)){
    auto_corr <- acf(asd_corr[[p]][,i], plot=F, type = "correlation", lag = length(asd_corr[[p]][,i]))$acf[,,1]
    if(NaN %in% auto_corr) auto_corr <- rep(1, length(asd_corr[[p]][,i]))
    asd_corr[[p]][,i] <- auto_corr
  }
}



td_corr <- td_data_lab_scale
for(p in 1:length(td_corr)){
  patient <- td_corr[[p]]
  for(i in 1:ncol(patient)){
    auto_corr <- acf(td_corr[[p]][,i], plot=F, type = "correlation", lag = length(td_corr[[p]][,i]))$acf[,,1]
    if(NaN %in% auto_corr) auto_corr <- rep(1, length(td_corr[[p]][,i]))
    td_corr[[p]][,i] <- auto_corr
  }
}


td_ROI_mean <- list_summary(td_corr, mean)
td_ROI_median <- list_summary(td_corr, median)
td_ROI_sd <- list_summary(td_corr, sd)
x_data <- c()
x_data <- as.data.frame(cbind(td_ROI_mean, td_ROI_median, td_ROI_sd, rep(0, length(td_data))))
names(x_data)[ncol(x_data)] <- 'label'


asd_ROI_mean <- list_summary(asd_corr, mean)
asd_ROI_median <- list_summary(asd_corr, median)
asd_ROI_sd <- list_summary(asd_corr, sd)
z_data <- c()
z_data <- as.data.frame(cbind(asd_ROI_mean, asd_ROI_median, asd_ROI_sd, rep(1, length(asd_data))))
names(z_data)[ncol(z_data)] <- 'label'

#####
u <- as.matrix(rbind(x_data[1:50,],z_data))
k <- dim(u)[2]
n0 <- 50
n1 <- 85


my_model <- glmnet(
  x = u[,1:(k-1)],
  y = u[,k],
  family = "gaussian",
  alpha = 0,
  maxit = 10000
)

sc <- predict(my_model, u[,1:(k-1)] , s = 0.05)
hist(sc[1:50], xlim = c(0,1), freq = F)
hist(sc[51:135], add = T, col = "red", freq = F)

aa <- Friedman_procedure(100,x_data = x_data[1:50,], z_data = z_data, permut = T )
k_obs_1 <- ks.test(sc[1:50] , sc[51:135])$statistic

k_obs_1 < quantile(aa, 0.95)


##############
set.seed(261198)
dim(x_data)
u_2 <- as.matrix(x_data)
label <- c(rep(0, 50), rep(1, 43))

n0 <- 50
n1 <- 43

model_2 <- glmnet(
  x = u_2[, 1:(k - 1)],
  y = label,
  family = "gaussian",
  alpha = 0.1,
  maxit = 100
)
sc_2 <- predict(model_2, u_2[, 1:(k - 1)] , s = 0.7)
hist(sc_2[1:50], xlim = c(0, 1))
hist(sc_2[51:93], add = T, col = "red")
k2 <- ks.test(sc_2[1:50] , sc_2[51:93])$statistic

bb <-
  Friedman_procedure(
    P = 100,
    x_data = x_data[1:50, ],
    z_data = x_data[51:93, ],
    permut = T
  )
k2 < quantile(bb, 0.95)
save(aa, file = "data/test_4_different.Rdata")
save(k_obs_1, file = "data/k_obs_different.Rdata")
save(bb, file = "data/test_4_same.Rdata")
save(k2, file = "data/k_obs_same.Rdata")

save(sc , file = "data/score_different.RData")
save(sc_2 , file = "data/score_same.RData")

#######
