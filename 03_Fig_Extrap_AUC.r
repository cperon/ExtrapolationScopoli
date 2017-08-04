##--------------------------------------------------------------------------------------------------------
## SCRIPT : Extrapolation vs AUC plots
##
## Author : Matthieu Authier & Clara P?ron
## Last update : 2017-04-30
## R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
##--------------------------------------------------------------------------------------------------------

lapply(c("ggplot2", "reshape"), library, character.only = TRUE)

theme_set(theme_bw(base_size = 12))

rm(list = ls())

WorkDir <- "C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/Tracking/results/GAM_models/2011/RawLocSim/_v2/extrapolation_analysis/v_2/"
OutDir <- paste(WorkDir, "output", sep = "/")
DataDir <- paste(WorkDir, "data", sep = "/")

extra <- read.csv(file = paste(DataDir, "AUC_extra.csv", sep = "/"), header = TRUE,
                  sep = ";", dec= ",")
extra$calib <- as.character(extra$site1) ; extra$calib <- ifelse(extra$calib == "Riou", "Marseille", extra$calib)
extra$pred <- as.character(extra$site2) ; extra$pred <- ifelse(extra$pred == "Riou", "Marseille", extra$pred)
extra$calib <- factor(extra$calib, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
extra$pred <- factor(extra$pred, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
extra$nom <- substr(as.character(extra$pred), 1, 1)

#plot with ggplot
theme_set(theme_bw(base_size = 12))

g <- ggplot(data = extra,
            aes(x = Extra, y = AUC, label = nom) #, color = pred
            ) +
  geom_point(size = 2, color = grey(0.6)) +
  geom_text(nudge_y = 0.025) +
  facet_wrap( ~ calib, nrow = 2) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  xlab(label = "% extrapolation")+ ylab(label = "AUC") +
  scale_y_sqrt(breaks = seq(0.4, 1.0, 0.1), limits = c(0.4, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 100.0, 20), limits = c(0, 100)) +
  # scale_color_manual(name = "", values = c(")) +
  theme(legend.position = "top", 
        legend.key.width = unit(2, "cm"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(lineheight = 0.8, face = "bold"), 
        axis.text = element_text(size = 10)
        )
g
ggsave(plot = g, 
       filename = paste(OutDir, "extra.png", sep = "/"), dpi = 600, units = "cm", 
       height = 10, width = 12
       )
