##--------------------------------------------------------------------------------------------------------
## SCRIPT : Analysis de l'extrapolation/interpolation entre covariables
##
## Authors : Matthieu Authier & Clara Peron
## Last update : 2017-09-04
##
## R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
## Copyright (C) 2016 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

lapply(c("WhatIf", "lpSolve"), library, character.only = TRUE)

rm(list = ls())

WorkDir <- getwd()
DataDir <- paste(WorkDir, "data", sep = "/")
OutDir <- paste(WorkDir, "output", sep = "/")

# function to prepare data
make_cfact_2 <- function(calibration_data, 
                         test_data, 
                         var_name = NULL,
                         howmany = 1,  
                         eps = 6
                         ) {
  
  ### custom code
  if(is.null(var_name)) { var_name = names(calibration_data) }
  
  ## standardize new data to predict from
  ### useful functions
  rescale <- function(x) { return((x - mean(x)) / sd(x)) }
  rescale2 <- function(ynew, y) { return((ynew - mean(y, na.rm = TRUE)) /(sd(y, na.rm = TRUE))) }
  # this simplifies computation A LOT!
  make_X <- function(calibration_data, test_data, var_name){
    X <- sapply(var_name,
                function(k) { rescale2(ynew = test_data[, k],
                                       y = calibration_data[, k]
                )}
    )
    X <- as.data.frame(X)
    names(X) <- var_name
    return(X)
  }
  ### standardize
  Xcal = make_X(calibration_data = calibration_data, test_data = calibration_data, var_name)
  Xtest = make_X(calibration_data = calibration_data, test_data = test_data, var_name)
  
  # Round the standardized oceano values
  Xcal <- round(Xcal, eps) ; Xtest <- round(Xtest, eps)
  
  # Remove duplicates 
  dup <- duplicated(Xcal[, var_name])
  Xcal <- Xcal[dup == FALSE, ] ; rm(dup)
  
  # rename rows
  row.names(Xcal) <- 1:nrow(Xcal)
  row.names(Xtest) <- 1:nrow(Xtest)
  
  # compute counterfactuals
  cf <- whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "hull", mc.cores = 1)
  return(cf)
}

# Load data
load(paste(DataDir, "shrw_2011_nodup.Rdata", sep = "/"))
load(paste(DataDir, "shrw_2012_nodup.Rdata", sep = "/"))

# Check correlations
sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id){nrow(subset(shrw_2011, Site == id))})
round(cor(subset(shrw_2011, Site == "Riou")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
round(cor(subset(shrw_2011, Site == "Porquerolles")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
round(cor(subset(shrw_2011, Site == "Giraglia")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
round(cor(subset(shrw_2011, Site == "Lavezzi")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)

### Run cfact for Riou 
riou <- make_cfact_2(calibration_data = subset(shrw_2011, Site == "Riou"),
                     test_data = subset(shrw_2012, Site == "Riou"),
                     var_name = c("Bathy", "Var_SST28", "logCHLA7", "SST28")
                     )
save.image(file = paste(OutDir, "ExtraTemporal_shrw2.RData.", sep = "/"), safe = "TRUE")

### Run cfact for Porquerolles
porque <- make_cfact_2(calibration_data = subset(shrw_2011, Site == "Porquerolles"),
                       test_data = subset(shrw_2012, Site == "Porquerolles"),
                       var_name = c("Bathy", "GBathy", "Var_SST28", "logCHLA28")
                       )
save.image(file = paste(OutDir, "ExtraTemporal_shrw2.RData.", sep = "/"), safe = "TRUE")

### Run cfact for Giraglia
giraglia <- make_cfact_2(calibration_data = subset(shrw_2011, Site == "Giraglia"),
                         test_data = subset(shrw_2012, Site == "Giraglia"),
                         var_name = c("Bathy", "GBathy", "logCHLA28", "SLA1")
                         )
save.image(file = paste(OutDir, "ExtraTemporal_shrw2.RData.", sep = "/"), safe = "TRUE")

### Run cfact for Lavezzi
lavezzi <- make_cfact_2(calibration_data = subset(shrw_2011, Site == "Lavezzi"),
                        test_data = subset(shrw_2012, Site == "Lavezzi"),
                        var_name = c("Bathy", "VEL28", "SLA7", "logCHLA28")
                        )
save.image(file = paste(OutDir, "ExtraTemporal_shrw2.RData.", sep = "/"), safe = "TRUE")

### concatenate results if a df
extra <- data.frame(Calibration = c("Riou 2011", "Porquerolles 2011", "Giraglia 2011", "Lavezzi 2011"),
                    Test = c("Riou 2012", "Porquerolles 2012", "Giraglia 2012", "Lavezzi 2012"),
                    Extrapolation = unlist(lapply(list(riou, porque, giraglia, lavezzi), function(x) { 100*round(mean(ifelse(x$in.hull, 0, 1)), 3) }))
                    )
extra
write.table(extra, file = paste(OutDir, "Scopoli_Temporalextrapolation.txt", sep = ";"), quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE
            )

save.image(file = paste(OutDir, "ExtraTemporal_shrw2.RData.", sep = "/"), safe = "TRUE")

