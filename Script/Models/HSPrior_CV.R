lmodHS <- nimbleCode({
  
  for (i in 1:n) {
    mu[i] <- inprod(X[i, 1:nLoc], beta[1:nLoc]) + beta0
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
  
  # mean predicted values for test data
  for(t in 1:nTest) {
    mu_pred[t] <- beta0 + inprod(Xtest[t, 1:nLoc], beta[1:nLoc])
  }
  
  for (j in 1:nLoc) {
    beta[j] ~ dnorm(0, sd = tau*lambda[j])
    lambda[j] ~ T(dt(mu=0, sigma=1, df=1), 0, Inf)
  }
  
  beta0 ~ dnorm(0, sd = 10)
  tau ~ T(dt(mu=0, sigma=1, df=1), 0, Inf)
  sigma   ~ dexp(1)
})
