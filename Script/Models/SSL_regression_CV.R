SSLglm <- nimbleCode({ 
  for (i in 1:n) {
    mu[i] <- beta0 + inprod(X[i, 1:nLoc], beta[1:nLoc])
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
  
  # mean predicted values for test data
  for(t in 1:nTest) {
    mu_pred[t] <- beta0 + inprod(Xtest[t, 1:nLoc], beta[1:nLoc])
  }
  
  # Priors for beta and reparameterization for Laplace prior
  for (j in 1:nLoc) {
    gamma[j] ~ dbern(theta) # Binary variable
    # Definition of lambda
    # lambda1 << lambda0
    lambda[j] <- gamma[j] * lambda1 + (1 - gamma[j]) * lambda0 # Mixing two lambdas
    tau_sq[j] ~ dexp(lambda[j]^2/2)
    beta[j] ~ dnorm(0, var = sigma^2 * tau_sq[j])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd = 10)
  theta ~ dbeta(a, b)
  sigma   ~ dexp(1)
  
})
