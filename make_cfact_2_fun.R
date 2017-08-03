make_cfact_2 <- function(calibration_data, test_data, var_name = NULL, howmany = 1, batch_size = 10000, eps = 6, brute_force = FALSE) {
  
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
  ### break the pb is small pieces
  wIf_bypart <- function(test, calib) {
    ### index of rows for the testing dataset
    outside <- matrix(NA, nrow = nrow(test), ncol = 5)
    ### loop to take several samples to approximate the convex hull
    for(k in 1:5){
      writeLines(paste("\t approximate convex hull:", k, " pass", sep = " "))
      ### take a subsample of the calibration dataset
      sub <- sample(1:nrow(calib), size = batch_size, replace = FALSE)
      ### this loop is needed to exhaust the rows of Xtest
      cf <- whatif(formula = NULL, 
                   data = calib[sub, ], 
                   cfact = test, 
                   choice = "distance", 
                   nearby = howmany
      )
      outside[, k] <- cf$sum.stat
    }
    return(apply(outside, 1, mean))
  }
  ### loop over the test data row
  if(!brute_force) {
    zz <- seq(0, nrow(Xtest), by = 5000)
    if(zz[length(zz)] < nrow(Xtest)) { zz <- c(zz, nrow(Xtest)) }
    return(
      do.call(what = "c",
              args = lapply(lapply(2:length(zz), function(j){ (zz[j-1] + 1):zz[j] }), 
                            function(x) { 
                              writeLines(paste("Processing rows", min(x), ":", max(x), sep = " "))
                              return(wIf_bypart(test = Xtest[x, ], calib = Xcal))
                            }
              )
      )
    )
  }
  ### brute force
  else {
    return(whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "distance", nearby = howmany))
  }
}
