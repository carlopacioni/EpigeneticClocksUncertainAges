library(data.table)
library(ggplot2)
library(rsample)
library(tidyverse)
library(MammalMethylClock)
source("./Script/Fit_epiClocks.R")

#### Prelim data prep ####
# using the R package subset of the data. >300 samples are available from Lu et al paper (data in 'mydata_GitHub.Rds')
# This dataset is probably a more realistic rapresentation of what available for wildlife
load("../mammalian-methyl-clocks-v1.1.0/jazoller96-mammalian-methyl-clocks-a9425df/TutorialData/methdatConsortium_subset.rda")
datAllSamp <- methdatConsortium_subset
load("../mammalian-methyl-clocks-v1.1.0/jazoller96-mammalian-methyl-clocks-a9425df/TutorialData/infoConsortium.rda")
infoAllSamp <- infoConsortium
anAge <- getAnAgeTable() %>%
  dplyr::filter(profiled == T) %>%
  dplyr::mutate(gestationYears = Gestation.Incubation..days. / 365) %>%
  dplyr::select(SpeciesLatinName, averagedMaturity.yrs, maxAgeCaesar, gestationYears)
infoAllSamp <- base::merge(infoAllSamp, anAge, by = "SpeciesLatinName", all.x = T, sort = F)

yxs.list <- alignDatToInfo(infoAllSamp, datAllSamp, "SID", "SID")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
ys$Tissue <- factor(ys$Tissue)
ys$SpeciesLatinName <- factor(ys$SpeciesLatinName)
hist(ys$Age)
ys$AgeLL <- fun_llin3.trans(ys$Age, ys$averagedMaturity.yrs)
hist(ys$AgeLL)

pathResults="./BaboonResultsRandom"
dir.create(pathResults, showWarnings = FALSE)
ggplot(ys, aes(Age)) + geom_histogram(breaks=seq(floor(min(ys$Age)), ceiling(max(ys$Age)), by=1))
ggsave(file.path(pathResults, "AgeDistributionBaboons.png"), 
       width = 12, height = 7, units = "cm")

# Note keeping the label 'AgeYrs' to get the functions to work, but this is LL transformed
meta <- data.table(RNR=seq_along(ys$AgeLL), AgeYrs=ys$AgeLL, 
                   AgeClass=cut(ys$Age, breaks=0:ceiling(max(ys$Age)))) 

nFolds <- 8
set.seed(13)
kfolds <- rsample::vfold_cv(data=cbind(meta, xs), v=nFolds, strata="AgeClass")

# debug(fit_clocks)
ls_fitRandom <- fit_clocks(nFolds, kfolds, ageTrans=TRUE, maturity=ys$averagedMaturity.yrs[1], 
                     lnstd="wide", glmnet_alpha=0.5, glmnet_lambda="1se",
                     #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                     #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                     niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                     niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                     propUncert=c(0.3, 0.4, 0.5), y_min="random", 
                     pathResults=pathResults, 
                     pathConvergence="./BaboonConvergenceRandom")
save(list=c("kfolds", "ls_fitRandom"), file = file.path(pathResults, "list_fitRandom.rda"))

ls_fit_informative <- fit_clocks(nFolds, kfolds, ageTrans=TRUE, maturity=ys$averagedMaturity.yrs[1], 
                     lnstd="wide", glmnet_alpha=0.5, glmnet_lambda="1se", cert=FALSE,
                     #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                     #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                     niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                     niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                     propUncert=c(0.3, 0.4, 0.5), y_min="informative",
                     pathResults="./BaboonResultsInformative", 
                     pathConvergence="./BaboonConvergenceInformative")
save(list=c("kfolds", "ls_fit_informative"), 
     file = file.path("./BaboonResultsInformative", "list_informative.rda"))

ls_fitRandomNarrow <- fit_clocks(nFolds, kfolds, ageTrans=TRUE, maturity=ys$averagedMaturity.yrs[1], 
                           lnstd="narrow", glmnet_alpha=0.5, glmnet_lambda="1se",
                           qInterval = 6, cert=FALSE,
                           #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                           #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                           niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                           niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                           propUncert=c(0.3, 0.4, 0.5), y_min="random", 
                           pathResults="./BaboonResultsRandomNarrow", 
                           pathConvergence="./BaboonConvergenceRandomNarrow")
save(list=c("kfolds", "ls_fitRandomNarrow"), 
     file = file.path("./BaboonResultsRandomNarrow", "list_fitRandomNarrow.rda"))

ls_fit_informativeNarrow <- fit_clocks(nFolds, kfolds, ageTrans=TRUE, maturity=ys$averagedMaturity.yrs[1], 
                               lnstd="narrow", glmnet_alpha=0.5, glmnet_lambda="1se",
                               qInterval = 4, cert=FALSE,
                               #niter_cert=c(250, 500), nburnin_cert=c(50, 200), # testing
                               #niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                               niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                               niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                               propUncert=c(0.3, 0.4, 0.5), y_min="informative",
                               pathResults="./BaboonResultsInformativeNarrow", 
                               pathConvergence="./BaboonConvergenceInformativeNarrow")
save(list=c("kfolds", "ls_fit_informativeNarrow"), 
     file = file.path("./BaboonResultsInformativeNarrow", "list_informativeNarrow.rda"))

#### Raw predictors ####
ls_fitRandomRaw <- fit_clocks(nFolds, kfolds, ageTrans=TRUE, maturity=ys$averagedMaturity.yrs[1], 
                           lnstd="wide", glmnet_alpha=0.5, glmnet_lambda="1se", PCR = FALSE,
                           niter_cert=c(250, 50), nburnin_cert=c(50, 20), # testing
                           niter_uncert=c(100, 100), nburnin_uncert=c(50, 50), 
                           #niter_cert=c(25000, 50000), nburnin_cert=c(5000, 20000),
                           #niter_uncert=c(100000, 100000), nburnin_uncert=c(50000, 50000), 
                           propUncert=0.3, y_min="random", HS=FALSE,
                           pathResults="./BaboonResultsRandomRaw", 
                           pathConvergence="./BaboonConvergenceRandomRaw")
save(list=c("kfolds", "ls_fitRandomRaw"), 
     file = file.path("./BaboonResultsRandomRaw", "list_fitRandomRaw.rda"))
