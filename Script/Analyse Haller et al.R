library(data.table)
library(nimble)
library(ggplot2)
library(MCMCvis)
library(glmnet)
library(tidyverse)
library(rsample)
library(parallel)
source("./Script/cvkfold_glmnet.R")
source("./Script/CrossValidate_certAge.R")
source("./Script/CrossValidate_uncertAge.R")
source("./Script/Fit_epiClocks.R")

#### Prelim data prep ####
meta <- fread(file="../Haller_et_al_2025_supp_mat_doi_10_5061_dryad_wstqjq2wq__v20250207/Extern_Samples_AgeingClock.txt")
meth <- fread(file="../Haller_et_al_2025_supp_mat_doi_10_5061_dryad_wstqjq2wq__v20250207/MethylationData_AgeingClock.txt")
finx <- as.matrix(t(meth[, -1]))
colnames(finx) <- meth[, CpG]
keep <- apply(finx, 2, function(x) mean(x)>=0.05 & mean(x)<=0.95 & sd(x)>=0.05)
sum(keep)
finx <- finx[, keep]
hist(finx)
meta[, AgeYrs := Age_in_Days/365]
meta <- meta[RNR %in% rownames(finx),]
hist(meta[, AgeYrs])

pathResults="./HallerResultsRandom"
dir.create(pathResults, showWarnings = FALSE)
ggplot(meta, aes(AgeYrs)) + 
  geom_histogram(breaks=seq(floor(meta[, min(AgeYrs)]), ceiling(meta[, max(AgeYrs)]), by=1))
ggsave(file.path(pathResults, "AgeDistributionHaller.png"), 
       width = 12, height = 7, units = "cm")

meta[, AgeClass := cut(AgeYrs, breaks=0:ceiling(max(AgeYrs)))]
meta[, nAgeClass := .N, by=AgeClass]
meta[, table(AgeClass)]

nFolds <- 8
set.seed(13)
kfolds <- rsample::vfold_cv(data=cbind(meta[, .(RNR, AgeYrs, AgeClass)], finx), 
                            v=nFolds, strata="AgeClass")

ls_fitrandom <- fit_clocks(nFolds, kfolds, ageTrans=FALSE, maturity=1, 
                     glmnet_alpha=0.2, glmnet_lambda="min", lnstd="wide",
                     #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                     #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                     niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                     niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                     propUncert=c(0.3, 0.4, 0.5), y_min="random", 
                     pathResults=pathResults, 
                     pathConvergence="./HallerConvergenceRandom")
save(list = c("kfolds", "ls_fitrandom"), 
     file = file.path(pathResults, "list_random.rda"))

ls_fit_informative <- fit_clocks(nFolds, kfolds, ageTrans=FALSE, maturity=1, 
                     glmnet_alpha=0.2, glmnet_lambda="min", lnstd="wide",
                     #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                     #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                     niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                     niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                     propUncert=c(0.3, 0.4, 0.5), y_min="informative",
                     pathResults="./HallerResultsInformative", 
                     pathConvergence="./HallerConvergenceInformative")
save(list = c("kfolds", "ls_fit_informative"), 
     file = file.path("./HallerResultsInformative", "list_informative.rda"))

ls_fitrandomNarrow <- fit_clocks(nFolds, kfolds, ageTrans=FALSE, maturity=1, 
                           glmnet_alpha=0.2, glmnet_lambda="min", lnstd="narrow",
                           qInterval = 1.5, cert=FALSE,
                           #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                           #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                           niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                           niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                           propUncert=c(0.3, 0.4, 0.5), y_min="random", 
                           pathResults="./HallerResultsRandomNarrow", 
                           pathConvergence="./HallerConvergenceRandomNarrow")
save(list = c("kfolds", "ls_fitrandomNarrow"), 
     file = file.path("./HallerResultsRandomNarrow", "list_randomNarrow.rda"))

ls_fit_informativeNarrow <- fit_clocks(nFolds, kfolds, ageTrans=FALSE, maturity=1, 
                                 glmnet_alpha=0.2, glmnet_lambda="min", lnstd="narrow",
                                 qInterval = 1.5, cert=FALSE,
                                 #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                                 #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                                 niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                                 niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                                 propUncert=c(0.3, 0.4, 0.5), y_min="informative",
                                 pathResults="./HallerResultsInformativeNarrow", 
                                 pathConvergence="./HallerConvergenceInformativeNarrow")
save(list = c("kfolds", "ls_fit_informativeNarrow"), 
     file = file.path("./HallerResultsInformativeNarrow", "list_informativeNarrow.rda"))

#### Raw predictors ####
ls_fitrandomRaw <- fit_clocks(nFolds, kfolds, ageTrans=FALSE, maturity=1, PCR=FALSE,
                           glmnet_alpha=0.2, glmnet_lambda="min", lnstd="wide",
                           #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                           #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                           niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                           niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                           propUncert=c(0.3), y_min="random", HS=FALSE,
                           pathResults="./HallerResultsRandomRaw", 
                           pathConvergence="./HallerConvergenceRandomRaw")
save(list = c("kfolds", "ls_fitrandomRaw"), 
     file = file.path(pathResults, "list_randomRaw.rda"))