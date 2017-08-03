##--------------------------------------------------------------------------------------------------------
## SCRIPT : Analysis de l'extrapolation/interpolation entre covariables
##
## Authors : Matthieu Authier & Clara Peron
## Last update : 2017-08-03
##
##R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
## Copyright (C) 2016 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

lapply(c("WhatIf", "lpSolve", "CLmapping"), 
       library, character.only = TRUE)

rm(list = ls())

setwd(WorkDir <- "C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/Tracking/results/GAM_models/2011/RawLocSim/_v2/extrapolation_analysis/v_2/")
DataDir <- paste(WorkDir, "shearw/data", sep = "/")
OutDir <- paste(WorkDir, "shearw/output", sep = "/")

# Source Functions
## Gower's distance
# function to automate somewhat computations over 
# all possible combinations of variables
  source(paste(WorkDir, 'shearw/make_cfact_fun.R', sep=''))
  source(paste(WorkDir, 'shearw/make_cfact_2_fun.R', sep=''))
  source(paste(WorkDir, 'shearw/kappa_std_fun.R', sep=''))

# function to prepare data
  scopoli_data <- function(df, var_name, eps = 3) {
    X <- df[, var_name]
    # Round the standardized oceano values
    X <- round(X, eps)
    # Remove duplicates 
    X <- X[duplicated(X[, var_name]) == FALSE, ]
    row.names(X) <- 1:nrow(X)
    return(X)
  }

# Load data
  load(paste(DataDir, "shrw_2011_nodup.Rdata", sep = "/"))

# Check correlations
  sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id){nrow(subset(shrw_2011, Site == id))})
  round(cor(subset(shrw_2011, Site == "Riou")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
  round(cor(subset(shrw_2011, Site == "Porquerolles")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
  round(cor(subset(shrw_2011, Site == "Giraglia")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)
  round(cor(subset(shrw_2011, Site == "Lavezzi")[, c('Bathy', 'GBathy', 'logCHLA7', 'SST28', 'Var_SST28', 'logCHLA28', 'VEL28', 'SLA7', 'SLA1')]), 2)

# RUN COUNTERFACTUALS  
  
  ###########################
  # Marseilles Is.          #
  ###########################

  # Check correlations
    round(cor(scopoli_data(df = subset(shrw_2011, Site == "Riou"), var_name = c("Bathy", "Var_SST28", "logCHLA7", "SST28"))), 2)
    #           Bathy Var_SST28 logCHLA7  VEL7
    # Bathy      1.00     -0.32    -0.62  0.43
    # Var_SST28 -0.32      1.00     0.37  0.04
    # logCHLA7  -0.62      0.37     1.00 -0.24
    # VEL7       0.43      0.04    -0.24  1.00
    sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id) {
      kappa_std(calibration_data = subset(shrw_2011, Site == "Riou"), 
                test_data = subset(shrw_2011, Site == id), 
                var_name = c("Bathy", "Var_SST28", "logCHLA7", "SST28")
                )
      }
      )
    # Riou Porquerolles     Giraglia      Lavezzi
    #  2.5          4.9          9.3          8.0
  
  # # Run cfact : OLD version 
  #   riou <- function() {
  #     sapply(c("Porquerolles", "Giraglia", "Lavezzi"), function(site) {
  #     make_cfact(calibration_data = subset(shrw_2011, Site == "Riou"),
  #                test_data = subset(shrw_2011, Site == site),
  #                var_name = c("Bathy", "Var_SST28", "logCHLA7", "SST28")
  #                )
  #   })
  #   }
  # 
  # # Results
  #   riou()
  #   # Porquerolles     Giraglia      Lavezzi 
  #   #    0.1868958    0.7724466    0.8188816 

    # NEW version
    riou <- function() {
      sapply(c("Porquerolles", "Giraglia", "Lavezzi"), function(site) {
        testdata <- subset(shrw_2011, Site == site)
        testdata$extra <- make_cfact_2(calibration_data = subset(shrw_2011, Site == "Riou"),
                                       test_data =testdata,
                                       var_name = c("Bathy", "Var_SST28", "logCHLA7", "SST28")
                                      )
        return(testdata)
      })
    }
    extra <- riou()
    
    # % extrapolation
    summary(extra[[1]]) ; hist(extra[[1]])
    
    
  #############################
  # Porquerolles Is.          #
  #############################

    # Check correlations
      round(cor(scopoli_data(df = subset(shrw_2011, Site == "Porquerolles"), var_name = c("Bathy", "Var_SST28", "logCHLA7", "logGBathy"))), 2)
      #           Bathy Var_SST28 logCHLA7 logGBathy
      # Bathy      1.00      0.30    -0.35     -0.50
      # Var_SST28  0.30      1.00     0.17     -0.30
      # logCHLA7  -0.35      0.17     1.00     -0.11
      # logGBathy -0.50     -0.30    -0.11      1.00
      sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id) {
        kappa_std(calibration_data = subset(shrw_2011, Site == "Porquerolles"), 
                  test_data = subset(shrw_2011, Site == id), 
                  var_name = c("Bathy", "Var_SST28", "logCHLA7", "logGBathy")
                  )
        }
        )
      # Riou Porquerolles     Giraglia      Lavezzi 
      #  6.5          3.2          9.2         20.9
      
    # Run cfact  
      porquerolles <- function() {
        sapply(c("Riou", "Giraglia", "Lavezzi"), function(site) {
        make_cfact(calibration_data = subset(shrw_2011, Site == "Porquerolles"),
                   test_data = subset(shrw_2011, Site == site),
                   var_name = c("Bathy", "Var_SST28", "logCHLA7", "logGBathy")
                   )
          })
      }

    porquerolles()
    #       Riou   Giraglia    Lavezzi 
    # 0.07057949 0.47669553 0.63837638 

  #############################
  # Giraglia Is.              #
  #############################
    
    # CHeck corraltions
      round(cor(scopoli_data(df = subset(shrw_2011, Site == "Giraglia"), var_name = c("Bathy", "logGBathy", "SLA7", "GSST7"))), 2)
      #           Bathy logGBathy  SLA7 GSST7
      # Bathy      1.00      0.15 -0.25 -0.23
      # logGBathy  0.15      1.00  0.05 -0.02
      # SLA7      -0.25      0.05  1.00  0.14
      # GSST7     -0.23     -0.02  0.14  1.00
      sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id) {
        kappa_std(calibration_data = subset(shrw_2011, Site == "Giraglia"), 
                  test_data = subset(shrw_2011, Site == id), 
                  var_name = c("Bathy", "logGBathy", "SLA7", "GSST7")
                  )
        }
        )
      # Riou Porquerolles     Giraglia      Lavezzi 
      #  2.4          2.9          1.7          4.1
    
    # Run cfact  
      giraglia <- function() {
        sapply(c("Riou", "Porquerolles", "Lavezzi"), function(site) {
        make_cfact(calibration_data = subset(shrw_2011, Site == "Giraglia"),
                   test_data = subset(shrw_2011, Site == site),
                   var_name = c("Bathy", "logGBathy", "SLA7", "GSST7")
                   )
          })
      }
      
      giraglia()
      #      Riou Porquerolles      Lavezzi 
      # 0.5247207    0.5906122    0.2433088

  #############################
  # Lavezzi Is.          #
  #############################

      # Check correlation
        round(cor(scopoli_data(df = subset(shrw_2011, Site == "Lavezzi"), var_name = c("Bathy", "logGBathy", "SLA7", "logCHLA7"))), 2)
        #           Bathy logGBathy SLA7 logCHLA7
        # Bathy      1.00      0.40 0.06    -0.24
        # logGBathy  0.40      1.00 0.05    -0.24
        # SLA7       0.06      0.05 1.00     0.12
        # logCHLA7  -0.24     -0.24 0.14     1.00
        
        sapply(c("Riou", "Porquerolles", "Giraglia", "Lavezzi"), function(id) {
          kappa_std(calibration_data = subset(shrw_2011, Site == "Giraglia"), 
                    test_data = subset(shrw_2011, Site == id), 
                    var_name = c("Bathy", "logGBathy", "SLA7", "logCHLA7")
          )
        }
        )
        #   Riou Porquerolles     Giraglia      Lavezzi 
        #  304.1         90.8          2.5          6.0
        
        # Run cfact
          lavezzi <- function() {
            sapply(c("Riou", "Porquerolles", "Giraglia"), function(site) {
            make_cfact(calibration_data = subset(shrw_2011, Site == "Lavezzi"),
                       test_data = subset(shrw_2011, Site == site),
                       var_name = c("Bathy", "logGBathy", "SLA7", "logCHLA7")
                       )
          })
          }
          
          lavezzi()
          #      Riou Porquerolles     Giraglia 
          # 0.8019495    0.4631418    0.1193783 

          
          