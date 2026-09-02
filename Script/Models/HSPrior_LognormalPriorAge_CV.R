HS_lognormalPriorAge <- nimbleCode({

  for(i in 1:n) {
    mu[i] <- beta0 + inprod(X[i, 1:nLoc], beta[1:nLoc])
    # Main likelihood
    y[i] ~ dnorm(mu[i], sd = sigma)
    
    # Transform in years, and use a lognorm prior if trans==TRUE 
    y_yr[i] <- llin3.inv(y[i], avMat, trans) 
    
    # Soft prior info 
    m_prior[i] ~ dlnorm(log(y_yr[i]), sd = sd_prior[i])
    
    # Inferred age for uncertain samples
    y_rep[i] ~ dnorm(mu[i], sd=sigma)
  }
  
  # predict values for test data
  for(t in 1:nTest) {
    mu_pred[t] <- beta0 + inprod(Xtest[t, 1:nLoc], beta[1:nLoc])
  }
  
  # Priors
  for(j in 1:nLoc) {
    beta[j] ~ dnorm(0, sd = tau*lambda[j])
    lambda[j] ~ T(dt(mu=0, sigma=1, df=1), 0, Inf)
  }
  
  beta0 ~ dnorm(0, sd = 10)
  sigma ~ dexp(1)
  tau ~ T(dt(mu=0, sigma=1, df=1), 0, Inf)
})
