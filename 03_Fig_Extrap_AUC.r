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

WorkDir <- "C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/Tracking/results/GAM_models/2011/RawLocSim/_v2/extrapolation_analysis/v_2/ExtrapolationScopoli"
OutDir <- paste(WorkDir, "output", sep = "/")
DataDir <- paste(WorkDir, "data", sep = "/")

# Load extrapolation file
  extra <- read.csv(paste(OutDir, "Extrap_sites.csv", sep="/"), sep=';', header=T)

  
# Load predictive ability table : TAB_val
  load('C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/Tracking/results/GAM_models/2011/RawLocSim/_v2/validation/TAB_val_inter_sites_2011.Rdata')

  extra$Calib_site <- as.character(extra$Calib_site) ; extra$Calib_site <- ifelse(extra$Calib_site == "Riou", "Marseille", extra$Calib_site)
  extra$Test_site <- as.character(extra$Test_site) ; extra$Test_site <- ifelse(extra$Test_site == "Riou", "Marseille", extra$Test_site)

  
  TAB_val$Site1 <- as.character(TAB_val$Site1) ; TAB_val$Site1 <- ifelse(TAB_val$Site1 == "Riou", "Marseille", TAB_val$Site1)
  TAB_val$Site2 <- as.character(TAB_val$Site2) ; TAB_val$Site2 <- ifelse(TAB_val$Site2 == "Riou", "Marseille", TAB_val$Site2)
  names(TAB_val) <- c('Calib_site', 'Test_site', 'AUC')
  
  TAB_val$Calib_site <- factor(TAB_val$Calib_site, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
  TAB_val$Test_site <- factor(TAB_val$Test_site, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
  TAB_val$nom <- substr(as.character(TAB_val$Test_site), 1, 1)

  extra$Calib_site <- factor(extra$Calib_site, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
  extra$Test_site <- factor(extra$Test_site, levels = c("Marseille", "Porquerolles", "Giraglia", "Lavezzi"))
  extra$nom <- substr(as.character(extra$Test_site), 1, 1)
  
  TAB_val$comb <- paste(TAB_val$Calib_site, TAB_val$Test_site, sep='_')
  extra$comb <- paste(extra$Calib_site, extra$Test_site, sep='_')
  
  TAB <- merge(TAB_val, extra, by='comb')
  
#plot with ggplot
theme_set(theme_bw(base_size = 12))

colo <- c('Riou', 'Giraglia', 'Porquerolles', 'Lavezzi')
colsite <- c('#91bfdb', '#fc8d59', '#1b7837', '#d73027')

g <- ggplot(data = TAB,
            aes(x = Extrapolation, y = AUC, label = Test_site.x) #, color = pred
            ) +
  geom_point(size = 2, color = 'black') +
  geom_text(nudge_y = 0.03, size=3.4) +
  facet_wrap( ~ Calib_site.x, nrow = 2) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  xlab(label = "% extrapolation")+ ylab(label = "AUC") +
  scale_y_continuous(breaks = seq(0.3, 1.0, 0.1), limits = c(0.3, 1.0)) +
  scale_x_continuous(breaks = seq(0.0, 100.0, 20), limits = c(-6, 100)) +
  # scale_color_manual(name = "", values = c(")) +
  theme(legend.position = "top", 
        legend.key.width = unit(2, "cm"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(lineheight = 0.8, face = "bold"), 
        axis.text = element_text(size = 10),
        strip.text = element_text(size=12)
        )
g
ggsave(plot = g, 
       filename = paste(OutDir, "extra.png", sep = "/"), dpi = 600, units = "cm", 
       height = 17, width = 22
       )
