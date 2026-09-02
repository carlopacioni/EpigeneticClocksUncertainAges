install.packages(c('glmnet', 'WGCNA', 'tidyr', 'snakecase'))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("preprocessCore")

install.packages("tidyverse")
install.packages("rsample")

install.packages("WGCNA")
install.packages("./mammalian-methyl-clocks-v1.1.0/jazoller96-mammalian-methyl-clocks-a9425df/MammalMethylClock-Package/MammalMethylClock_1.0.0.tar.gz", 
                 repos = NULL, type = "source")
