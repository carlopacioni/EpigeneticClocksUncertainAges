library(data.table)
library(ggplot2)
library(cowplot)

path <- "C:/Users/cp0d/OneDrive - Department of Energy, Environment and Climate Action/Waterfowl/Duck genetics 2024/Epigenetics/Writing"

#### Histograms ####
bab_plot <- readRDS("./BaboonResultsRandom/AgeDistributionBaboons.rds")
hal_plot <- readRDS("./HallerResultsRandom/AgeDistributionHaller.rds")
histoGrop <- plot_grid(bab_plot, hal_plot, ncol=2, labels="AUTO")
save_plot(file.path(path, "histoGrop.png"), histoGrop)

#### Scatter plots ####

# Haller
load(file=file.path("./HallerResultsInformativeNarrow", "list_informativeNarrow.rda"))
haller_inf_narrow_0.5 <- ls_fit_informativeNarrow[[3]][["0.5"]]

dt_uncert <- haller_inf_narrow_0.5[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                                                  MAE=round(mean(abs(AgeYrs - AgePred)), 2)),
                               by=.(Model, Data)]

hal_info_narr_plot <- ggplot(haller_inf_narrow_0.5[!is.na(AgePred) & Model == "SSL", ], 
                             aes(AgeYrs, AgePred)) + geom_point() + 
  geom_abline(slope=1) + 
  geom_smooth(method="lm") + 
  facet_grid(Model~Data, labeller=as_labeller(c("Test"="Test data", "SSL"="SSL",
                                            "UncertainAge"= "Uncertain ages"))) +
  labs(x="True age (years)", y="Predicted age") +
  geom_text(data=dt_uncert[Model == "SSL",], col="red", aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                                           label=paste0("r=", r, "\nMAE=", MAE)))

ggsave("./HallerResultsInformativeNarrow/uncert_p_05_plot_final.png", 
       width=15, height=8, units='cm')


load(file=file.path("./HallerResultsRandomNarrow", "list_randomNarrow.rda"))
haller_rand_narrow_0.5 <- ls_fitrandomNarrow[[3]][["0.5"]]

dt_uncert <- haller_rand_narrow_0.5[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                                                      MAE=round(mean(abs(AgeYrs - AgePred)), 2)),
                                   by=.(Model, Data)]

hal_rand_narr_plot <- ggplot(haller_rand_narrow_0.5[!is.na(AgePred) & Model == "SSL", ], 
                             aes(AgeYrs, AgePred)) + geom_point() + 
  geom_abline(slope=1) + 
  geom_smooth(method="lm") + 
  facet_grid(Model~Data, labeller=as_labeller(c("Test"="Test data", "SSL"="SSL",
                                                  "UncertainAge"= "Uncertain ages"))) +
  labs(x="True age (years)", y="Predicted age") +
  geom_text(data=dt_uncert[Model == "SSL",], col="red", 
            aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                label=paste0("r=", r, "\nMAE=", MAE)))

ggsave("./HallerResultsRandomNarrow/uncert_p_05_RandNar_plot_final.png", 
       width=15, height=8, units='cm')

hallGrop <- plot_grid(hal_info_narr_plot, hal_rand_narr_plot, ncol=1, labels="AUTO")
save_plot(file.path(path, "hallGrop.png"), hallGrop, #base_asp = 1
            base_width=5.9, base_height=5.9)

# Baboon
rm(list=ls())

load(file=file.path("./BaboonResultsInformativeNarrow", "list_informativeNarrow.rda"))

Baboon_inf_narrow_0.5 <- ls_fit_informativeNarrow[[3]][["0.5"]]

dt_uncert <- Baboon_inf_narrow_0.5[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                                                       MAE=round(mean(abs(AgeYrs - AgePred)), 2)),
                                    by=.(Model, Data)]

bab_info_narr_plot <- ggplot(Baboon_inf_narrow_0.5[!is.na(AgePred) & Model == "SSL", ], aes(AgeYrs, AgePred)) + geom_point() + 
  geom_abline(slope=1) + 
  geom_smooth(method="lm") + 
  facet_grid(Model~Data, labeller=as_labeller(c("Test"="Test data", "SSL"="SSL",
                                                  "UncertainAge"= "Uncertain ages"))) +
  labs(x="True age (years)", y="Predicted age") +
  geom_text(data=dt_uncert[Model == "SSL",], col="red", 
            aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
                label=paste0("r=", r, "\nMAE=", MAE)))

ggsave("./BaboonResultsInformativeNarrow/uncert_p_05_infNar_plot_final.png", 
       width=15, height=8, units='cm')

load( file=file.path("./BaboonResultsRandomNarrow", "list_fitRandomNarrow.rda"))

Baboon_rand_narrow_0.5 <- ls_fitRandomNarrow[[3]][["0.5"]]

dt_uncert <- Baboon_rand_narrow_0.5[!is.na(AgePred), .(r=round(cor(AgeYrs, AgePred), 2),
                                                      MAE=round(mean(abs(AgeYrs - AgePred)), 2)),
                                   by=.(Model, Data)]

bab_rand_narr_plot <- ggplot(Baboon_rand_narrow_0.5[!is.na(AgePred) & Model == "SSL", ], aes(AgeYrs, AgePred)) + geom_point() + 
  geom_abline(slope=1) + 
  geom_smooth(method="lm") + 
  facet_grid(Model~Data, labeller=as_labeller(c("Test"="Test data", "SSL"="SSL",
                                                  "UncertainAge"= "Uncertain ages"))) +
  labs(x="True age (years)", y="Predicted age") +
  geom_text(data=dt_uncert[Model == "SSL",], col="red", 
            aes(x=Inf, y=-Inf, vjust=-1, hjust=1.1, 
            label=paste0("r=", r, "\nMAE=", MAE)))

ggsave("./BaboonResultsRandomNarrow/uncert_p_05_randNar_plot_final.png", 
       width=15, height=8, units='cm')

babGrop <- plot_grid(bab_info_narr_plot, bab_rand_narr_plot, ncol=1, labels="AUTO")
save_plot(file.path(path, "babGrop.png"), hallGrop, #base_asp = 1
          base_width=5.9, base_height=5.9)