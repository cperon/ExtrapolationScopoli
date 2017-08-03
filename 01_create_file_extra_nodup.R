rm(list=ls())

# Load libraries
libr <- c("raster", "mgcv", "gdata", "CLmapping", "ROCR", "plyr", "ggplot2", "dplyr", "gridExtra", "grid")
lapply(libr, require, character.only=T)  

source('results/GAM_models/2011/RawLocSim/_v2/imputation_fun.R')

# Load data
  # Training dataset 2011
  load('results/GAM_models/2011/RawLocSim/_v2/data/track_oceano_2011-cleaned.Rdata')
  ALL2011 <- track_oceano
  #ALL2011 <- ALL
  ALL2011$Site <- as.character(ALL2011$Site)    
  ALL2011$Site <- ifelse(ALL2011$Site=='Frioul', 'Riou', ALL2011$Site)      
  ALL2011$Site <- as.factor(ALL2011$Site)  
  table(ALL2011$Site)

# Loop on site

shrw_2011 <- NULL
for(site in levels(ALL2011$Site)){
  print(site)
  
  X <- ALL2011[ALL2011$Site %in% site,]
  X1 <- X
  
  # plotCorse()
  # with(X1[X1$OccFo==0,], points(Longitude, Latitude, cex=0.3, pch=16))   
  # with(X1[X1$OccFo==1,], points(Longitude, Latitude, cex=0.3, pch=16, col='red'))   
  
  X1 <- drop.levels(X1)
  
  # Remove NAs by data imputation
  X1 <- data_imputation(X1)
  
  # Format the training dataset
  if(site=='Riou')
  {
    X1$Bathy <- ifelse(X1$Bathy>2000, 2000, X1$Bathy)
    X1$logCHLA7 <- ifelse(X1$logCHLA7<0.2, 0.2, X1$logCHLA7)
    X1$logCHLA1 <- ifelse(X1$logCHLA1<0.2, 0.2, X1$logCHLA1)
    X1$logCHLA28 <- ifelse(X1$logCHLA28<0.2, 0.2, X1$logCHLA28)
    
  }
  
  if(site=='Porquerolles')
  {
    X1$Bathy <- ifelse(X1$Bathy>2000, 2000, X1$Bathy)
    X1$logCHLA7 <- ifelse(X1$logCHLA7<0.2, 0.2, X1$logCHLA7)
    X1$logCHLA1 <- ifelse(X1$logCHLA1<0.2, 0.2, X1$logCHLA1)
    X1$logCHLA28 <- ifelse(X1$logCHLA28<0.2, 0.2, X1$logCHLA28)
  }
  
  if(site=='Giraglia')
  {
    X1$Bathy <- ifelse(X1$Bathy>1000, 1000, X1$Bathy)
    X1$SLA1 <-  ifelse(X1$SLA1>0.15, 0.15, X1$SLA1)
  }
  
  if(site=='Lavezzi')
  {
    X1$Bathy <- ifelse(X1$Bathy>1000, 1000, X1$Bathy)
    X1$VEL28 <- ifelse(X1$VEL28>0.3, 0.3, X1$VEL28)
    X1$SLA7 <-  ifelse(X1$SLA7>0.15, 0.15, X1$SLA7)
    #X1 <- X1[X1$logCHLA28<0.07,]
    #X1 <- X1[X1$ID != 'ID24L',]
  }
  
  varkeep <- c('Longitude', 'Latitude', 'year', 'Site', 'OccFo', 
               'Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')
  X2 <- X1[, varkeep]
  X2$rownames <- 1:nrow(X2)
  X2$Dup <- duplicated(X2[,c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')])
  X3 <- X2[X2$Dup==FALSE,]
  shrw_2011 <- rbind(shrw_2011, X3)
}

shrw_2011$Dup <- NULL
shrw_2011$rownames <- NULL

save(shrw_2011, file='results/GAM_models/2011/RawLocSim/_v2/extrapolation_analysis/v_2/data/shrw_2011_nodup.Rdata')
