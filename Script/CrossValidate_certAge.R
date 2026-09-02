#' Fit a SSL and HS model to a train dataset and extract the mean predicted value 
#' 
#'   of a test dataset (usually for cross validation)
#' k index of the k fold 
#' kfold_sel selection of k fold (outout from rsample::vfold_cv). It is expecting 
#'    IDs in column 1, Age in column 2 and strata variable in Column 3
#' PCvar proportion of variance from the principal components retained 
#' lambda1 value of landba1 parameter for the SSL
#' #' SSL Logical: whether to fit the SSL model 
#' HS  Logical: whether to fit the HS model 
#' niter, nburnin, nchains vectors of length=2 for mcmc settings, where 
#'       the first value is used for SSL and the second for HS
#' pathModels path to the text file where the models are stored
#' pathConvergence where to write the coefficient values (and convergence diagnostics)

cvKfold_certAge <- function(k, kfold_sel, PCvar, lambda1,
                            SSL, HS, PCR, niter, nburnin, nchains, 
                              pathModels, pathConvergence) {
  require(data.table)
  require(nimble)
  require(ggplot2)
  require(MCMCvis)
  require(glmnet)
  require(data.table)
  require(tidyverse)
  require(rsample)
  
  #### Prep data ####
  train <- analysis(kfold_sel$splits[[k]])
  test <- assessment(kfold_sel$splits[[k]])
  
  y <- train[, AgeYrs]
  n <- length(y)
  
  ssize_train <- train[, table(AgeClass)]
  ssize_test <- test[, table(AgeClass)]
  
  write.csv(ssize_train, 
            file = file.path(pathConvergence, paste0("Cert_AgeClass_N_fold", k,"_train.csv"))) 
  
  write.csv(ssize_test, 
            file = file.path(pathConvergence, paste0("Cert_AgeClass_N_fold", k,"_test.csv"))) 
  
  if(PCR) {
    # Prep train data
    pca <- prcomp(train[, -(1:3)])
    prVar <- pca$sdev^2/sum(pca$sdev^2)
    cutoff <- which(cumsum(prVar)>PCvar)[1]
    pcax <- pca$x[, seq_len(cutoff)]
    
    # Prep test data
    testx <- predict(object=pca, test[, -(1:3)])
    testx <- testx[, seq_len(cutoff)]
    testy <- test[, AgeYrs]
  } else {
    # Prep train data
    # variable names are unfortunate, but I was too lazy to edit the code...
    pcax <- train[, -(1:3)]
    
    # Prep test data
    testx <- test[, -(1:3)]
    testy <- test[, AgeYrs]
  }
  
  
  #### PC SSL regression ####
  if(SSL) {
  source(file.path(pathModels, "SSL_regression_CV.R"))
  dataList <- list(X=pcax, y=y, Xtest=testx)
  constList <- list(n=n, nLoc=ncol(pcax), a=2, b=2, nTest=length(testy),
                    lambda0=lambda1 * 20, lambda1=lambda1)
  initList <- list(beta=rep(0, ncol(pcax)), sigma=1)
  
  SSLModelpca <- nimbleModel(code=SSLglm, constants=constList, 
                             data=dataList, inits=initList)
  
  # Compile model
  CSSLModel <- compileNimble(SSLModelpca)
  
  # Configure and compile MCMC
  mcmcConf <- configureMCMC(CSSLModel, enableWAIC=TRUE)
  mcmcConf$addMonitors(c('beta0', 'beta', 'sigma', 'theta', 'mu_pred'))
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project=CSSLModel)
  
  # Run MCMC
  SSL_time <- system.time(
    SSLpcaSamples <- runMCMC(Cmcmc, niter=niter[1], nburnin=nburnin[1], nchains=nchains[1])
  )
  
  res <- MCMCsummary(SSLpcaSamples)
  message(paste("Ran SSL in", round(SSL_time[3], 1), "secs"))
  write.csv(res[grep("^beta", rownames(res)),], 
            file=file.path(pathConvergence, paste0("SSL_PC_cert_betas_fold", k,".csv")))
  MCMCtrace(SSLpcaSamples, filename = file.path(pathConvergence, 
                                    paste0("SSL_PC_cert_betas_fold", k,".pdf")), 
            open_pdf = FALSE)
  
  # res
  # max(res[grep("^beta", rownames(res)), "Rhat"])
  # res[grep("^beta", rownames(res)),]
  # y_mod <- res[grep("^y\\[", rownames(res)), "50%"][s]
  # y_inf <- res[grep("^y_rep", rownames(res)), "50%"][s]
  mu_pred <- res[grep("^mu_pred", rownames(res)), "50%"]
  SSL_res <- test[, .(RNR, AgeYrs)]
  SSL_res[, ':='(AgePred=mu_pred, PC=TRUE, Model="SSL", Data="Test", UncertAge=FALSE)]
  } else {
    SSL_res <- test[, .(RNR, AgeYrs)]
    SSL_res[, ':='(AgePred=NA, PC=TRUE, Model="SSL", Data="Test", UncertAge=FALSE)]
  }
  
  
  #### PC regression Horseshoe prior ####
  if(HS) {
  source(file.path(pathModels, "HSprior_CV.R"))
  dataList <- list(X=pcax, y=as.vector(y), Xtest=testx)
  constList <- list(n=n, nLoc=ncol(pcax), nTest=length(testy))
  initList <- list(beta=rep(0, ncol(pcax)), sigma=1, lambda=rep(0.1, ncol(pcax)), 
                   tau=0.1)
  
  
  HSModelPca <- nimbleModel(code=lmodHS, constants=constList, 
                            data=dataList, inits=initList)
  
  # Compile model
  CHSModel <- compileNimble(HSModelPca)
  params <- c('beta0', 'beta', 'sigma', 'tau', #'lambda', 
              'mu_pred')
  
  # Configure and compile MCMC
  mcmcConf <- configureMCMC(CHSModel, enableWAIC=TRUE)
  mcmcConf$addMonitors(params)
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project=CHSModel)
  
  # Run MCMC
  system.time(
    HSpcaSamples <- runMCMC(Cmcmc, niter=niter[2], nburnin=nburnin[2], nchains=nchains[2])
  )
  
  res <- MCMCsummary(HSpcaSamples)
  # max(res$Rhat)
  # res[grep("^beta", rownames(res)),]
  # #MCMCtrace(HSpcaSamples)
  write.csv(res[grep("^beta", rownames(res)),], 
            file=file.path(pathConvergence, paste0("HS_PC_cert_betas_fold", k, ".csv")))
  MCMCtrace(HSpcaSamples, filename = file.path(pathConvergence, 
                                               paste0("HS_PC_cert_betas_fold", k,".pdf")), 
            open_pdf = FALSE)
  
  mu_pred <- res[grep("^mu_pred", rownames(res)), "50%"]
  HS_res <- test[, .(RNR, AgeYrs)]
  HS_res[, ':='(AgePred=mu_pred, PC=TRUE, Model="HS", Data="Test", UncertAge=FALSE)]
  } else {
    HS_res <- test[, .(RNR, AgeYrs)]
    HS_res[, ':='(AgePred=NA, PC=TRUE, Model="HS", Data="Test", UncertAge=FALSE)]
  }
  
  return(rbind(SSL_res, HS_res))
}