#' Fit a SSL and HS model to a train dataset and extract the mean predicted value 
#' 
#' The FUN will select p samples to select a min(y) <= y_min <= y_true which is 
#'   used to put a lognormal prior on age as explained below (see probq). 
#' 
#' k index of the k fold 
#' kfold_sel selection of k fold (outout from rsample::vfold_cv) It is expecting 
#'        IDs in column 1, Age in column 2 and strata variable in Column 3
#' p the proportion of samples with uncertain age
#' std standard deviation of the lognorm prior distribution on age. 
#' 
cvKfold_uncertAge <- function(k, kfold_sel, p, probq, std, 
                              y_min, qInterval,
                              trans, maturity, PCvar, lambda1,
                              SSL, HS, PCR,
                              niter, nburnin, nchains,
                              pathModels, pathConvergence) {
  
  require(data.table)
  require(nimble)
  require(ggplot2)
  require(MCMCvis)
  require(glmnet)
  require(data.table)
  require(tidyverse)
  require(rsample)
  require(MammalMethylClock)
  source("./Script/Nimble_FUNs/LogLineraAgeInvTransformation.R")
  
  #### Helper FUN lognormal params ####
  lognormal_mean <- function(qmin, q975, probq) {
    z1 <- qnorm(probq)
    z2 <- qnorm(0.975)
    
    sigma <- (log(q975) - log(qmin)) / (z2 - z1)
    mu <- log(q975) - z2 * sigma
    
    logmu <- exp(mu + sigma^2 / 2)
    return(logmu)
  }
  
  lognormal_sigma <- function(qmin, q975, probq) {
    
    z1 <- qnorm(probq)
    z2 <- qnorm(0.975)
    
    sigma <- (log(q975) - log(qmin)) / (z2 - z1)
    return(sigma)
  }
  
  #### Prep data ####
  train <- analysis(kfold_sel$splits[[k]])
  test <- assessment(kfold_sel$splits[[k]])
  
  y <- train[, AgeYrs]
  n <- length(y)
  
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
  
  # Select samples with uncertain y
  s <- rep(FALSE, length(y))
  if(trans) y_yr <- fun_llin3.inv(y, maturity = maturity) else y_yr <- y
  set.seed(22)
  # sample with age < 1 yr are known
  for(i in seq_along(y)) {
    s[i] <- ifelse(y_yr[i]<=1, s[i], as.logical(rbinom(1, size=1, p)))
  }
  
  # Get the sample sizes and the uncertain ages
  ssize_train <- train[, table(AgeClass)]
  tmp <- cbind(train[, .(AgeYrs, AgeClass)], s)
  ssize_train <- merge(ssize_train, tmp[, .(Uncert=sum(s)), by=AgeClass], all.x=TRUE)
  #ssize_test <- test[, table(AgeClass)]
  
  write.csv(ssize_train, 
            file = file.path(pathConvergence, 
                  paste0("Uncert_AgeClass_N_fold", k,"_train_", p,".csv")), 
            row.names=FALSE) 
  
  # write.csv(ssize_test, 
  #           file = file.path(pathConvergence, paste0("Cert_AgeClass_N_fold", k,"_test.csv"))) 
  
  # Prep to fit the model
  y_obs <- y
  y_min_Nim <- rep(-Inf, length=length(y))
  
  # find out which value is 1 yr if transformed
  if(trans) one_yr <- fun_llin3.trans(1, maturity) else one_yr <- 1
  
  # randomly take a value between 1 and the actual age
  for(i in seq_along(y[s])) {
    if(y_min=="informative") {
      y_min_Nim[s][i] <-  seq(one_yr, y[s][i], length=10)[9]
    } else {
      if(y_min=="random") {
          y_min_Nim[s][i] <- sample(seq(one_yr, y[s][i], length=10), size=1)
      } else {
        stop("The argument 'y_min' can only be either 'random' or 'informative'")
      }
    }
  }
  
  y_obs[s] <- NA
  y_init <- y
  y_init[!s] <- NA
  
  # Start with uninformative prior (use actual age yr for mean and large sd)
  m_prior <- if(trans) fun_llin3.inv(y, maturity) else y
  s_prior <- rep(20, n)
  
  # Now give informative prior for missing y.
  if(std == "wide") { # the q0.975 is the max(age) + 1
    if(trans) maxAge <- fun_llin3.inv(max(y), maturity) + 1 else
      maxAge <- max(y) + 1
    maxAge <- rep(maxAge, sum(s))
  } else { # # the q0.975 is the min between qInterval and (y_min - max(age) + 1)
    if(trans) maxAge <- fun_llin3.inv(y_min_Nim[s], maturity) + 
                            min(qInterval, fun_llin3.inv(max(y), maturity) + 1 -
                                  fun_llin3.inv(y_min_Nim[s], maturity)) else
      maxAge <- y_min_Nim[s] + min(qInterval, max(y) + 1 - y_min_Nim[s])
  }
  
  # Find the mean and sigma of a lognorm distribution using where y_min is the 
  #       probq (0.025) and maxAge is q0.975 
  # z <- qnorm(probq)
  if(trans) {
    m_prior[s] <- mapply(lognormal_mean, qmin=fun_llin3.inv(y_min_Nim[s], maturity), 
                         q975=maxAge, MoreArgs = list(probq=probq))
    s_prior[s] <- mapply(lognormal_sigma, qmin=fun_llin3.inv(y_min_Nim[s], maturity), 
                         q975=maxAge, MoreArgs = list(probq=probq))
  } else {
    m_prior[s] <- mapply(lognormal_mean, qmin=y_min_Nim[s], q975=maxAge, 
                         MoreArgs = list(probq=probq))
    s_prior[s] <- mapply(lognormal_sigma, y_min_Nim[s], q975=maxAge, 
                         MoreArgs = list(probq=probq))
  }
  
  #### SSL LogNormPA ####
  if(SSL) {
  source(file.path(pathModels, "SSL_LognormalPriorAgeYr_CV.R")) # (SSL_LogNPA)
  
  dataList <- list(X=pcax, y=y_obs, m_prior=m_prior, Xtest=testx)
  constList <- list(n=n, nLoc=ncol(pcax), nTest=length(testy), a=2, b=2, 
                    sd_prior=s_prior, trans=trans, avMat=maturity,
                    lambda0=lambda1 * 20, lambda1=lambda1)
  initList <- list(beta=rep(0, ncol(pcax)), sigma=1, y= y_init)
  
  
  SSL_LNPA_Modelpca <- nimbleModel(code=SSL_LogNPA, constants=constList, 
                                   data=dataList, inits=initList)
  
  # Compile model
  CSSL_LNPA_Modelpca <- compileNimble(SSL_LNPA_Modelpca)
  
  # Configure and compile MCMC
  mcmcConf <- configureMCMC(CSSL_LNPA_Modelpca, enableWAIC=TRUE)
  mcmcConf$addMonitors(c('beta0', 'beta', 'sigma', 'theta', 'y_rep', 'mu_pred'))
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project=CSSL_LNPA_Modelpca)
  
  # Run MCMC
  SSL_LNPA_time <- system.time(
    CSSL_LNPA_pcaSamples <- runMCMC(Cmcmc, niter=niter[1], nburnin=nburnin[1], nchains=nchains[1])
  )
  message(paste("Ran SSL_LNPA in", round(SSL_LNPA_time[3], 1), "secs"))
  res <- MCMCsummary(CSSL_LNPA_pcaSamples)
  # res
  # max(res[grep("^beta", rownames(res)), "Rhat"])
  # res[grep("^beta", rownames(res)),]
  #y_mod <- res[grep("^y\\[", rownames(res)), "50%"][s]
  
  write.csv(res[grep("^beta", rownames(res)),], 
            file=file.path(pathConvergence, 
                      paste0("SSL_LNPA_PC_uncert_betas_fold", k, "_p", p,".csv")))
  MCMCtrace(CSSL_LNPA_pcaSamples, filename = file.path(pathConvergence, 
                      paste0("SSL_LNPA_PC_uncert_betas_fold", k,"_p", p,".pdf")), 
            open_pdf = FALSE)
  
  y_inf <- res[grep("^y_rep", rownames(res)), "50%"][s]
  mu_pred <- res[grep("^mu_pred", rownames(res)), "50%"]
  tmp <- test[, .(RNR, AgeYrs)]
  tmp[, ':='(AgePred=mu_pred, PC=TRUE, Model="SSL", Data="Test", 
             UncertAge=TRUE, PropUncertain=p)]
  SSL_LNPA_res <- train[s, .(RNR, AgeYrs)]
  SSL_LNPA_res[, ':='(AgePred=y_inf, PC=TRUE, Model="SSL", Data="UncertainAge", 
                      UncertAge=TRUE, PropUncertain=p)]
  SSL_LNPA_res <- rbind(SSL_LNPA_res, tmp)
  } else {
    tmp <- test[, .(RNR, AgeYrs)]
    tmp[, ':='(AgePred=NA, PC=TRUE, Model="SSL", Data="Test", 
               UncertAge=TRUE, PropUncertain=p)]
    SSL_LNPA_res <- train[s, .(RNR, AgeYrs)]
    SSL_LNPA_res[, ':='(AgePred=NA, PC=TRUE, Model="SSL", Data="UncertainAge", 
                        UncertAge=TRUE, PropUncertain=p)]
    SSL_LNPA_res <- rbind(SSL_LNPA_res, tmp)
  }
  
  #### HS LogNormPA ####
  if(HS) {
  source(file.path(pathModels, "HSPrior_LognormalPriorAge_CV.R")) # (HS_lognormalPriorAge)
  
  dataList <- list(X=pcax, y=y_obs, m_prior=m_prior, Xtest=testx)
  constList <- list(n=n, nLoc=ncol(pcax), y_min=y_min_Nim, sd_prior=s_prior, 
                    trans=trans, avMat=maturity,
                    nTest=length(testy))
  initList <- list(beta=rep(0, ncol(pcax)), sigma=1, lambda=rep(0.1, ncol(pcax)), 
                   tau=0.1, y=y_init)
  
  HS_LNPA_Modelpca <- nimbleModel(code=HS_lognormalPriorAge, constants=constList, 
                                  data=dataList, inits=initList)
  
  # Compile model
  CHS_LNPA_Modelpca <- compileNimble(HS_LNPA_Modelpca)
  
  # Configure and compile MCMC
  mcmcConf <- configureMCMC(CHS_LNPA_Modelpca, enableWAIC=TRUE)
  mcmcConf$addMonitors(c('beta0', 'beta', 'sigma', 'tau', #'lambda', 
                         'y_rep', 'mu_pred'))
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project=CHS_LNPA_Modelpca)
  
  # Run MCMC
  HS_LNPA_time <- system.time(
    CHS_LNPA_pcaSamples <- runMCMC(Cmcmc, niter=niter[2], nburnin=nburnin[2], nchains=nchains[2])
  )
  message(paste("Ran HS_LNPA in", round(SSL_LNPA_time[3], 1), "secs"))
  res <- MCMCsummary(CHS_LNPA_pcaSamples)
  
  # max(res[grep("^beta", rownames(res)), "Rhat"])
  # res[grep("^beta", rownames(res)),]
  # res[grep("^y", rownames(res)),]
  # res[grep("^y", rownames(res)), "Rhat"][s]
  
  write.csv(res[grep("^beta", rownames(res)),], 
            file=file.path(pathConvergence, 
                           paste0("HS_LNPA_PC_uncert_betas_fold", k, "_p", p, ".csv")))
  MCMCtrace(CHS_LNPA_pcaSamples, filename = file.path(pathConvergence, 
                            paste0("HS_LNPA_PC_uncert_betas_fold", k, "_p", p,".pdf")), 
            open_pdf = FALSE)
  
  y_inf <- res[grep("^y_rep", rownames(res)), "50%"][s]
  mu_pred <- res[grep("^mu_pred", rownames(res)), "50%"]
  tmp <- test[, .(RNR, AgeYrs)]
  tmp[, ':='(AgePred=mu_pred, PC=TRUE, Model="HS", Data="Test", 
             UncertAge=TRUE, PropUncertain=p)]
  HS_LNPA_res <- train[s, .(RNR, AgeYrs)]
  HS_LNPA_res[, ':='(AgePred=y_inf, PC=TRUE, Model="HS", Data="UncertainAge", 
                     UncertAge=TRUE, PropUncertain=p)]
  HS_LNPA_res <- rbind(HS_LNPA_res, tmp)
  } else {
    tmp <- test[, .(RNR, AgeYrs)]
    tmp[, ':='(AgePred=NA, PC=TRUE, Model="HS", Data="Test", 
               UncertAge=TRUE, PropUncertain=p)]
    HS_LNPA_res <- train[s, .(RNR, AgeYrs)]
    HS_LNPA_res[, ':='(AgePred=NA, PC=TRUE, Model="HS", Data="UncertainAge", 
                       UncertAge=TRUE, PropUncertain=p)]
    HS_LNPA_res <- rbind(HS_LNPA_res, tmp)
  }
  
  return(rbind(SSL_LNPA_res, HS_LNPA_res))
  }