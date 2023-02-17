```{r power info , warning = FALSE}
power_info <- function(P, k, increment){
  
  prop_ks <- rep(NA, P)
  for(i in 1:P){
    x <- Take_sample_normal(n = n0, mu= mu, sigma = sigma, label = 0)
    z <- Take_sample_normal(n = n1, mu=  mu2, sigma = sigma2, label = 1)
    u <- rbind(x,z) 
    p_model <- glm(label ~ ., data = u)
    x_scores <- predict(p_model, x)
    z_scores <- predict(p_model, z) 
    true_kolm <- ks.test(x_scores, z_scores, alternative = "two.sided")$statistic
    kk <- Friedman_procedure(P,x_data = x, z_data = z)
    prop_ks[i] <- true_kolm > quantile(kk,1-alpha)
  }
  data <-  mean(prop_ks)
  
  
  return(data)
}

```



ll <- seq(from = 0 , to = max_in , by = 0.3 )
powers <- c()
distance <- c()
for (i in (ll)){
  mu2 <- mu + i
  print(mu2)
  sigma2 <- sigma
  distance <- c(distance,Wasserstein_distance(mu, mu2, sigma, sigma2))
  powers <- c(powers,power_info(100, k = length(mu),increment = i))
}
plot(distance,powers,col = 'darkblue' , type = "o" , lwd = 2)s