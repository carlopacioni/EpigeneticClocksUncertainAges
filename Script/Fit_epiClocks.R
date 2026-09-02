#' Fit a SSL and HS model to a train dataset and extract the mean predicted value
#' of a test dataset using k folds. 
#' 
#' The FUN will 
#' 1.	Fit the data with known ages using glmnet 
#' 2.	Fit the data with known ages using glmnet using Principal Component (PC) 
#'    as predictors using alpha=0.5 and lambda.min
#' 3.	Fit the data with known ages using a SSL and/or HS (optionally using PC as predictors)
#' 4.	Fit the data with simulated uncertain ages (for each propUncert value) 
#'    using SSL and HS (optionally using PC as predictors). In this analysis, 
#'    the model-inferred age of the individuals with uncertain age is also returned.
#'   
#' For point 4, the FUN selects propUncert samples to select y_min so that 
#'   1 <= y_min <= y_true where y is the age. y_min is 
#'   used to put a lognormal prior on age as explained below (see probq). 
#'   
#' y_min  How the value y_min should be selected. Options include 
#' - "random" A 10 equal interval vector between 1 year old and the actual age 
#'     of the individuals sampled is generated. One value is randomly selected 
#'     and used as y_min to construct a lognormal prior (see below)
#' - "informative" as above but the before-last
#'     of the 10 equal interval values is always selected (i.e. y_min is always 
#'     close to the true value)
#'     
#'   
#' Convergence diagnostics: 
#'  1. tablessize_p: Sample size for each p value, where 
#'                   .id is the fold 
#'                   AgeClass is the age class
#'                   Freq is the total sample size in the training dataset
#'                   Uncert is the number of samples with uncertain age
#'  2. Rhat:  Beta coefficients with Rhat and n.eff. Lat two cols report the fold
#'            and p value (NA when not applicable)
#'  3. Rhat_short: sum of number of instances where Rhat is >1.05 (models are SSL or HS.
#'                 LNPA stands for LogNormal Prior on Age. That is the fit with 
#'                 uncertain ages)
#'   4. Cert_AgeClass_N_fold_k_t*: sample size for the train and test dataset for 
#'                                 fold k
#'   5. *(un)cert_betas_fold_k.csv: Beta coefficient estimates and diagnostics with 
#'                             (un)certain age for fold k. * either SSL or HS model
#'   6. *(un)cert_betas_fold_k.pdf: as above but trace and density plots
#'   
#' nFolds Number of folds 
#' kfold_sel selection of k fold (outout from rsample::vfold_cv) It is expecting 
#'        IDs in column 1, Age in column 2 and strata variable in Column 3
#' ageTrans Whether age transformation should be applied. Currently only 
#' log-linear v3  transformation is applied (via MammalMethylClock::fun_llin3.trans 
#' or its inverse in Nimble).
#' SSL Logical: whether to fit the SSL model 
#' HS  Logical: whether to fit the HS model 
#' PCR whether a Principal Component Regression should be used
#' lambda1 value of lambda1 parameter for the SSL. lambda0 is 20 x lambda1
#' PCvar The cumulative proportion of variance explained by the PC that needs to 
#'   be retained
#' maturity Age (in years) of when sexual maturity is reached
#' lnstd Whether the std of the lognorm prior distribution should be 'wide' (a 
#'       value such as the interval between y_min and the max(Age) + 1 yrs 
#'       corresponds to the 
#'       0.025 and 0.975 quantiles. 'narrow' selects a std in such a way that 
#'       the 0.975 quantiles is qInterval from y_min.
#' propUncert the proportion of samples with uncertain age
#' probq the probability of the quantile of the left hand side of the lognorm distribution.
#'   The y_min value selected will be considered the quantile with probq. 
#' y_min  How the minimum age should be selected for uncertain age. Options 
#'     include "random" or "informative".  
#' qInterval The interval (in years) between Q(probq) and Q(0.975). This is only 
#'          considered when lnstd="narrow"
#' glmnet_alpha, glmnet_lambda The alpha fit for the glmnet elastic net fit when 
#'    PCR=FALSE. The fit wtih PC is alsways glmnet_alpha=0.5, glmnet_lambda="min"
#' cert Logical. Whether the analysis with certain ages should be conducted. When
#'       multiple analysis are condicted to experiment with different setting of
#'       the uncertain age analysis, this can be set to FALSE to save computation time.
#' niter_x, nburnin_x, nchains vectors of length=2 for mcmc settings, where 
#'       the first value is used for SSL and the second for HS (respectively for 
#'       the data with certain and uncretain ages)
#' pathModels path to the text file where the models are stored
#' pathConvergence where to write the coefficient values (and convergence diagnostics)
#' pathModels path to the text file where the models are stored


fit_clocks <- function(nFolds, kfolds, ageTrans=FALSE, SSL=TRUE, HS=TRUE, PCR=TRUE,
                       lambda1=0.05, PCvar=0.999, maturity=1, lnstd=c("wide", "narrow"),
                       propUncert=c(0.3, 0.4, 0.5), probq=0.025, 
                       y_min=c("random", "informative"), qInterval=3,
                       glmnet_alpha=0.5, glmnet_lambda="1se", cert=TRUE,
                       niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                       niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                       nchains=c(2, 2),
                       pathResults, pathConvergence, pathModels="./Script/Models") {
  
  require(data.table)
  require(nimble)
  require(ggplot2)
  require(MCMCvis)
  require(glmnet)
  require(tidyverse)
  require(rsample)
  require(parallel)
  require(MammalMethylClock)
  source("./Script/cvkfold_glmnet.R")
  source("./Script/CrossValidate_certAge.R")
  source("./Script/CrossValidate_uncertAge.R")
  source("./Script/Nimble_FUNs/LogLineraAgeInvTransformation.R")
  
  dir.create(pathResults, showWarnings=FALSE)
  dir.create(pathConvergence, showWarnings=FALSE)
  
  #### Helper FUN ####
  collectRhat <- function(lf, p=NA) {
    ltables <- lapply(lf, fread)
    nms <- sapply(lapply(basename(lf), str_split, pattern="_", simplify=TRUE), 
                  function(x) paste(x[,1:2], collapse="_"))
    if(is.na(p)) {
      folds <- sapply(lapply(basename(lf), 
                             function(x) str_split(x, pattern="_", simplify=TRUE)), 
                      function(x) sub(".cvs$", "", x[,ncol(x)]))
    } else {
      folds <- sapply(lapply(basename(lf), 
                             function(x) str_split(x, pattern="_", simplify=TRUE)), 
                      function(x) x[,ncol(x) - 1])
    }
    
    ltables <- lapply(seq_along(ltables), function(i, tab, fold) tab[[i]][, Fold:=fold[i]], 
                      tab=ltables, fold=folds)
    names(ltables) <- nms
    table <- rbindlist(ltables, idcol=TRUE)
    table[, p:=p]
    subtable <- table[Rhat>1.05,]
    return(subtable)
  }
  
#### Fit glmnet ####
glmnet_res <- lapply(seq_len(nFolds), cvKfold_glmnet, kfold_sel=kfolds, PCvar=PCvar,
                     glmnet_alpha=glmnet_alpha, glmnet_lambda=glmnet_lambda)
glmnet_res <- do.call(rbind, glmnet_res)
if(ageTrans) {
  glmnet_res[, AgeYrs:=fun_llin3.inv(AgeYrs, maturity=maturity)]
  glmnet_res[, AgePred:=fun_llin3.inv(AgePred, maturity=maturity)]
}

dt <- glmnet_res[, .(r=round(cor(AgeYrs, AgePred), 2),
                 MAE=round(median(abs(AgeYrs - AgePred)), 2)), 
                 by=PC]

ggplot(glmnet_res, aes(AgeYrs, AgePred)) + geom_point() + geom_abline(slope=1) + 
  geom_smooth(method="lm") + facet_grid(PC~.) +
  geom_text(data=dt, col="red", aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                         label=paste0("r = ", r, "\nMAE = ", MAE)))

ggsave(file.path(pathResults, "glmnet_plot.png"), width=15, height=10, units='cm')

#### Fit SSL and HS when ages are certain ####
# debug
# debug(cvKfold_certAge)
# cv4fold_certAge <- cvKfold_certAge(k=1, kfold_sel=kfolds)

cl <- makeCluster(nFolds)
on.exit(stopCluster(cl))
if(cert) {
lcv4fold_certAge <- parLapply(cl=cl, seq_len(nFolds), fun=cvKfold_certAge, lambda1=lambda1,
                              kfold_sel=kfolds, SSL=SSL, HS=HS, PCR=PCR, PCvar=PCvar,
                              niter=niter_cert, nburnin=nburnin_cert, nchains=nchains,
                              pathModels=pathModels, pathConvergence=pathConvergence)

cv4fold_certAge <- do.call(rbind, lcv4fold_certAge)
if(ageTrans) {
  cv4fold_certAge[, AgeYrs:=fun_llin3.inv(AgeYrs, maturity=maturity)]
  cv4fold_certAge[, AgePred:=fun_llin3.inv(AgePred, maturity=maturity)]
}

dt_cert <- cv4fold_certAge[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                               MAE=round(median(abs(AgeYrs - AgePred)), 2)), 
                           by=Model]

ggplot(cv4fold_certAge[!is.na(AgePred),], aes(AgeYrs, AgePred)) + geom_point() + geom_abline(slope=1) + 
  geom_smooth(method="lm") + facet_grid(Model~.) +
  geom_text(data=dt_cert, col="red", aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                             label=paste0("r = ", r, "\nMAE = ", MAE))) 

ggsave(file.path(pathResults, "Cert_plot.png"), width=15, height=10, units='cm')
} else {
  cv4fold_certAge <- NULL
}
#### Fit SSL and HS when ages are uncertain ####
# debug
# debug(cvKfold_uncertAge)
# cv4fold_uncertAge <- cvKfold_uncertAge(k=8, kfold_sel=kfolds)

cv4fold_uncertAge_res <- vector("list", length = length(propUncert))
lsubtable <- vector("list", length = length(propUncert) + 1)

for(p in propUncert) {
  lcv4fold_uncertAge <- parLapply(cl=cl, seq_len(nFolds), fun=cvKfold_uncertAge, 
                                  kfold_sel=kfolds, p=p, probq=probq, std=lnstd,
                                  y_min=y_min, qInterval=qInterval,
                                  trans=ageTrans, maturity=maturity, PCvar=PCvar, 
                                  lambda1=lambda1, SSL=SSL, HS=HS, PCR=PCR,
                                  niter=niter_uncert, nburnin=nburnin_uncert, 
                                  nchains=nchains, pathModels=pathModels,
                                  pathConvergence=pathConvergence)
  
  cv4fold_uncertAge <- do.call(rbind, lcv4fold_uncertAge)
  if(ageTrans) {
    cv4fold_uncertAge[, AgeYrs:=fun_llin3.inv(AgeYrs, maturity=maturity)]
    cv4fold_uncertAge[, AgePred:=fun_llin3.inv(AgePred, maturity=maturity)]
  }
  
  dt_uncert <- cv4fold_uncertAge[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                                     MAE=round(median(abs(AgeYrs - AgePred)), 2)),
                                 by=.(Model, Data)]

  ggplot(cv4fold_uncertAge[!is.na(AgePred), ], aes(AgeYrs, AgePred)) + geom_point() + 
    geom_abline(slope=1) + 
    geom_smooth(method="lm") + facet_grid(Model~Data) +
    geom_text(data=dt_uncert, col="red", aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                                  label=paste0("r = ", r, "\nMAE = ", MAE)))
  
  ggsave(file.path(pathResults, paste0("uncert_p_", p, "_plot.png")), 
         width=15, height=15, units='cm')
  
  # Diagnostics 
  lf <- list.files(pathConvergence, paste0("fold[0-9]+","_p", p,".csv$"), full.names=TRUE)
  lsubtable[[which(propUncert == p)]] <- collectRhat(lf, p)
  
  lss <- list.files(pathConvergence, paste0("train_", p, ".csv"), full.names=TRUE)
  ltablessize <- lapply(lss, fread)
  tablessize <- rbindlist(ltablessize, idcol=TRUE)
  write.csv(tablessize, file.path(pathConvergence, paste0("tablessize_", p, ".csv")), 
            row.names=FALSE)
  
  cv4fold_uncertAge_res[[which(propUncert == p)]] <- cv4fold_uncertAge
}

names(cv4fold_uncertAge_res) <- propUncert
if(cert){
lf <- list.files(pathConvergence, paste0("fold[0-9]+.csv$"), full.names=TRUE)
lsubtable[[length(lsubtable)]] <- collectRhat(lf, p=NA)
}
subtable <- rbindlist(lsubtable)
write.csv(subtable, file.path(pathConvergence, "Rhat.csv"), 
          row.names=FALSE)

write.csv(subtable[, .N, by=.(.id, Fold, p)], file.path(pathConvergence, 
                                                     "Rhat_short.csv"), 
          row.names=FALSE)

return(list(glmnet_res, cv4fold_certAge, cv4fold_uncertAge_res))
}