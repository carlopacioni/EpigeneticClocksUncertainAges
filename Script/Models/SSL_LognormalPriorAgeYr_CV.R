SSL_LogNPA <- nimbleCode({
  
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
  
  # means predicted values for test data
  for(t in 1:nTest) {
    mu_pred[t] <- beta0 + inprod(Xtest[t, 1:nLoc], beta[1:nLoc])
  }
  
  # Priors
  for(j in 1:nLoc) {
    gamma[j] ~ dbern(theta) # Binary variable
    # Definition of lambda
    # lambda1 << lambda0
    lambda[j] <- gamma[j] * lambda1 + (1 - gamma[j]) * lambda0 # Mixing two lambdas
    tau_sq[j] ~ dexp(lambda[j]^2/2)
    beta[j] ~ dnorm(0, var = sigma^2 * tau_sq[j])
  }
  
  beta0 ~ dnorm(0, sd = 10)
  sigma ~ T(dnorm(0, sd=10),0, )   # half-normal prior on sigma
  theta ~ dbeta(a, b)
})