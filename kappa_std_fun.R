kappa_std <- function(calibration_data, test_data, var_name, eps = 3) {
  rescale2 <- function(ynew, y) { return((ynew - mean(y, na.rm = TRUE)) /(sd(y, na.rm = TRUE))) }
  X <- sapply(var_name,
              function(k){ rescale2(ynew = test_data[, k], y = calibration_data[, k]) }
  )
  X <- as.data.frame(X)
  names(X) <- var_name
  # Round the standardized oceano values
  X <- round(X, eps)
  
  # Remove duplicates 
  X <- X[duplicated(X[, var_name]) == FALSE, ]
  row.names(X) <- 1:nrow(X)
  
  # condition number
  round(kappa(X), 1)
}