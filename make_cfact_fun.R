

### Gower's distance
# function to automate somewhat computations over 
# all possible combinations of variables
make_cfact <- function(calibration_data, test_data, var_name = NULL, batch_size = 5000, brute_force = FALSE, howmany = 3, index = FALSE, eps = 3) {
  whatif <- function (formula = NULL, data, cfact, range = NULL, freq = NULL, 
                      nearby = 1, distance = "gower", miss = "list", choice = "hull", 
                      return.inputs = FALSE, return.distance = FALSE, ...) 
  {
    #print("Preprocessing data ...")
    if (!((is.character(cfact) && is.vector(cfact) && length(cfact) == 
           1) || is.data.frame(cfact) || (is.matrix(cfact) && !is.character(cfact)))) {
      stop("'cfact' must be either a string, a R data frame, or a R non-character matrix")
    }
    if (is.character(cfact)) {
      cfact <- read.table(cfact)
    }
    if (dim(cfact)[1] == 0) {
      stop("no counterfactuals supplied: 'cfact' contains zero rows")
    }
    if (!any(complete.cases(cfact))) {
      stop("there are no cases in 'cfact' without missing values")
    }
    if ("(Intercept)" %in% dimnames(cfact)[[2]]) {
      cfact <- cfact[, -(which(dimnames(cfact)[[2]] == "(Intercept)"))]
    }
    if (is.list(data) && !(is.data.frame(data))) {
      if (!((("formula" %in% names(data)) || ("terms" %in% 
                                              names(data))) && (("data" %in% names(data)) || ("model" %in% 
                                                                                              names(data))))) {
        stop("the list supplied to 'data' is not a valid output object")
      }
      tt <- terms(data)
      attr(tt, "intercept") <- rep(0, length(attr(tt, "intercept")))
      if ("data" %in% names(data)) {
        if (is.data.frame(data$data)) {
          data <- model.matrix(tt, model.frame(tt, data = data$data, 
                                               na.action = NULL))
        }
        else {
          data <- model.matrix(tt, model.frame(tt, data = eval(data$data, 
                                                               envir = .GlobalEnv), na.action = NULL))
        }
      }
      else {
        data <- model.matrix(tt, data = data$model)
      }
      if (!(is.matrix(data))) {
        stop("observed covariate data could not be extracted from output object")
      }
      rm(tt)
    }
    else {
      if (!((is.character(data) && is.vector(data) && length(data) == 
             1) || is.data.frame(data) || (is.matrix(data) && 
                                           !is.character(data)))) {
        stop("'data' must be either a string, a R data frame, a R non-character matrix, or an output object")
      }
      if (is.character(data)) {
        data <- read.table(data)
      }
    }
    if (dim(data)[1] == 0) {
      stop("no observed covariate data supplied: 'data' contains zero rows")
    }
    if (!any(complete.cases(data))) {
      stop("there are no cases in 'data' without missing values")
    }
    if (!(is.null(formula))) {
      if (identical(class(formula), "formula")) {
        if (!(is.data.frame(as.data.frame(data)))) {
          stop("'data' must be coercable to a data frame in order to use 'formula'")
        }
        if (!(is.data.frame(as.data.frame(cfact)))) {
          stop("'cfact' must be coercable to a data frame in order to use 'formula'")
        }
        formula <- update.formula(formula, ~. - 1)
        ttvar <- all.vars(formula)
        for (i in 1:length(ttvar)) {
          if (!(ttvar[i] %in% dimnames(data)[[2]])) {
            stop("variable(s) in 'formula' either unlabeled or not present in 'data'")
          }
          if (!(ttvar[i] %in% dimnames(cfact)[[2]])) {
            stop("variable(s) in 'formula' either unlabeled or not present in 'cfact'")
          }
        }
        rm(ttvar)
        data <- model.matrix(formula, data = model.frame(formula, 
                                                         as.data.frame(data), na.action = NULL))
        cfact <- model.matrix(formula, data = model.frame(formula, 
                                                          as.data.frame(cfact), na.action = NULL))
      }
      else {
        stop("'formula' must be of class 'formula'")
      }
    }
    if (!(identical(complete.cases(cfact), rep(TRUE, dim(cfact)[1])))) {
      cfact <- na.omit(cfact)
      print("Note:  counterfactuals with missing values eliminated from cfact")
    }
    if (is.data.frame(data)) {
      if (is.character(as.matrix(data))) {
        stop("observed covariate data not coercable to numeric matrix due to character column(s)")
      }
      data <- suppressWarnings(data.matrix(data))
    }
    else {
      data <- data.matrix(as.data.frame(data))
    }
    if (is.data.frame(cfact)) {
      if (is.character(as.matrix(cfact))) {
        stop("counterfactual data not coercable to numeric matrix due to character column(s)")
      }
      cfact <- suppressWarnings(data.matrix(cfact))
    }
    else {
      cfact <- data.matrix(as.data.frame(cfact))
    }
    if (!(is.matrix(data) && is.numeric(data))) {
      stop("observed covariate data not coercable to numeric matrix")
    }
    if (!(is.matrix(cfact) && is.numeric(cfact))) {
      stop("counterfactual data not coercable to numeric matrix")
    }
    na.fail(cfact)
    if (!identical(ncol(cfact), ncol(data))) {
      stop("number of columns of 'cfact' and 'data' are not equal")
    }
    if (!(is.null(range))) {
      if (!(is.vector(range) && is.numeric(range))) {
        stop("'range' must be a numeric vector")
      }
      if (!identical(length(range), ncol(data))) {
        stop("length of 'range' does not equal number of columns of 'data'")
      }
    }
    if (!(is.null(freq))) {
      if (!(is.vector(freq) && is.numeric(freq))) {
        stop("'freq' must be a numeric vector")
      }
      na.fail(freq)
    }
    if (!(is.null(nearby))) {
      if (!(is.numeric(nearby) && is.vector(nearby) && length(nearby) == 
            1 && nearby >= 0)) {
        stop("'nearby' must be numeric, greater than or equal to 0, and a scalar")
      }
    }
    if (!(identical(miss, "list") || identical(miss, "case"))) {
      stop("'miss' must be either ''case'' or ''list''")
    }
    if (!(identical(distance, "gower") || identical(distance, 
                                                    "euclidian"))) {
      stop("'distance' must be either ''gower'' or ''euclidian''")
    }
    if (!(identical(choice, "both") || identical(choice, "hull") || 
          identical(choice, "distance"))) {
      stop("'choice' must be either ''both'', ''hull'', or ''distance''")
    }
    if (!(is.logical(return.inputs))) {
      stop("'return.inputs' must be logical, i.e. either TRUE or FALSE")
    }
    if (!(is.logical(return.distance))) {
      stop("'return.distance' must be logical, i.e. either TRUE or FALSE")
    }
    n = nrow(data)
    convex.hull.test <- function(x, z) {
      n <- nrow(x)
      k <- ncol(x)
      m <- nrow(z)
      A <- rbind(t(x), rep(1, n))
      C <- c(rep(0, n))
      D <- c(rep("=", k + 1))
      hull = rep(0, m)
      for (i in 1:m) {
        B <- c(z[i, ], 1)
        lp.result <- lp(objective.in = C, const.mat = A, 
                        const.dir = D, const.rhs = B)
        if (lp.result$status == 0) 
          hull[i] <- 1
      }
      hull <- as.logical(hull)
      return(hull)
    }
    calc.gd <- function(dat, cf, range) {
      n <- nrow(dat)
      m <- nrow(cf)
      dat = t(dat)
      dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
      for (i in 1:m) {
        temp <- abs(dat - cf[i, ])/range
        if (any(range == 0)) {
          temp[is.nan(temp)] <- 0
          temp[temp == Inf] <- NA
        }
        dist[i, ] <- colMeans(temp, na.rm = T)
      }
      return(t(dist))
    }
    calc.ed <- function(dat, cf) {
      n <- nrow(dat)
      m <- nrow(cf)
      dat <- t(dat)
      dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
      for (i in 1:m) {
        temp <- (dat - cf[i, ])^2
        dist[i, ] <- (colSums(temp))
      }
      return(t(dist))
    }
    geom.var <- function(dat, rang) {
      n <- nrow(dat)
      dat <- t(dat)
      ff <- function(x) {
        temp <- abs(dat - x)/rang
        if (any(rang == 0)) {
          temp[is.nan(temp)] <- 0
          temp[temp == Inf] <- NA
        }
        tmp <- sum(colMeans(temp, na.rm = TRUE))
        return(tmp)
      }
      sum.gd.x <- sum(apply(dat, 2, ff), na.rm = TRUE)
      gv.x <- (0.5 * sum.gd.x)/(n^2)
      return(gv.x)
    }
    calc.cumfreq <- function(freq, dist) {
      m <- length(freq)
      n <- ncol(dist)
      res <- matrix(0, n, m)
      for (i in 1:m) res[, i] <- (colSums(dist <= freq[i]))/nrow(dist)
      return(res)
    }
    if (identical(miss, "list")) {
      data <- na.omit(data)
      n <- nrow(data)
    }
    if ((choice == "both") | (choice == "hull")) {
      #print("Performing convex hull test ...")
      test.result <- convex.hull.test(x = na.omit(data), z = cfact)
    }
    if ((choice == "both") | (choice == "distance")) {
      #print("Calculating distances ....")
      if (identical(distance, "gower")) {
        samp.range <- apply(data, 2, max, na.rm = TRUE) - 
          apply(data, 2, min, na.rm = TRUE)
        if (!is.null(range)) {
          w <- which(!is.na(range))
          samp.range[w] <- range[w]
        }
        if (identical(TRUE, any(samp.range == 0))) {
          print("Note:  range of at least one variable equals zero")
        }
        dist <- calc.gd(dat = data, cf = cfact, range = samp.range)
      }
      else {
        dist <- calc.ed(dat = na.omit(data), cf = cfact)
      }
      #print("Calculating the geometric variance...")
      if (identical(distance, "gower")) {
        gv.x <- geom.var(dat = data, rang = samp.range)
      }
      else {
        gv.x <- 0.5 * mean(calc.ed(dat = na.omit(data), cf = na.omit(data)))
      }
      if (identical(miss, "case") && identical(distance, "euclidian")) {
        summary <- colSums(dist <= nearby * gv.x) * (1/nrow(na.omit(data)))
      }
      else {
        summary <- colSums(dist <= nearby * gv.x) * (1/n)
      }
      #print("Calculating cumulative frequencies ...")
      if (is.null(freq)) {
        if (identical(distance, "gower")) {
          freqdist <- seq(0, 1, by = 0.05)
        }
        else {
          min.ed <- min(dist)
          max.ed <- max(dist)
          freqdist <- round(seq(min.ed, max.ed, by = (max.ed - 
                                                        min.ed)/20), 2)
        }
      }
      else {
        freqdist <- freq
      }
      cumfreq <- calc.cumfreq(freq = freqdist, dist = dist)
      dimnames(cumfreq) <- list(seq(1, nrow(cfact), by = 1), 
                                freqdist)
    }
    #print("Finishing up ...")
    if (return.inputs) {
      if (choice == "both") {
        if (return.distance) {
          out <- list(call = match.call(), inputs = list(data = data, 
                                                         cfact = cfact), in.hull = test.result, dist = t(dist), 
                      geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
        }
        else {
          out <- list(call = match.call(), inputs = list(data = data, 
                                                         cfact = cfact), in.hull = test.result, geom.var = gv.x, 
                      sum.stat = summary, cum.freq = cumfreq)
        }
      }
      if (choice == "distance") {
        if (return.distance) {
          out <- list(call = match.call(), inputs = list(data = data, 
                                                         cfact = cfact), dist = t(dist), geom.var = gv.x, 
                      sum.stat = summary, cum.freq = cumfreq)
        }
        else {
          out <- list(call = match.call(), inputs = list(data = data, 
                                                         cfact = cfact), geom.var = gv.x, sum.stat = summary, 
                      cum.freq = cumfreq)
        }
      }
      if (choice == "hull") {
        out <- list(call = match.call(), inputs = list(data = data, 
                                                       cfact = cfact), in.hull = test.result)
      }
    }
    else {
      if (choice == "both") {
        if (return.distance) {
          out <- list(call = match.call(), in.hull = test.result, 
                      dist = t(dist), geom.var = gv.x, sum.stat = summary, 
                      cum.freq = cumfreq)
        }
        else {
          out <- list(call = match.call(), in.hull = test.result, 
                      geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
        }
      }
      if (choice == "distance") {
        if (return.distance) {
          out <- list(call = match.call(), dist = t(dist), 
                      geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
        }
        else {
          out <- list(call = match.call(), geom.var = gv.x, 
                      sum.stat = summary, cum.freq = cumfreq)
        }
      }
      if (choice == "hull") {
        out <- list(call = match.call(), in.hull = test.result)
      }
    }
    class(out) <- "whatif"
    return(invisible(out))
  }
  
  
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
  Xcal <- round(Xcal, eps)
  Xtest <- round(Xtest, eps)
  
  # Remove duplicates 
  dup <- duplicated(Xcal[, var_name])
  Xcal <- Xcal[dup == F, ] ; rm(dup)
  row.names(Xcal) <- 1:nrow(Xcal)
  dup <- duplicated(Xtest[, var_name])
  Xtest <- Xtest[dup == F, ] ; rm(dup)
  row.names(Xtest) <- 1:nrow(Xtest)
  
  # compute counterfactuals
  ### break the pb is small pieces
  wIf_bypart <- function(Xtest, Xcal, batch_size) {
    ### index of rows for the testing dataset
    outside <- 1:nrow(Xtest)
    
    ### loop to take several samples to approximate the convex hull
    for(k in 1:howmany){
      ### take a subsample of the calibration dataset
      sub <- sample(1:nrow(Xcal), size = batch_size, replace = FALSE)
      ### this loop is needed to exhaust the rows of Xtest
      if(nrow(Xtest[outside, ]) != 0) {
        cf <- whatif(formula = NULL, 
                     data = Xcal[sub, ], 
                     cfact = Xtest[outside, ]
        )
        outside <- outside[which(cf$in.hull == FALSE)]
      }
    }
    return(outside)
  }
  
  ### loop over the test data row
  if(!brute_force) {
    # extra <- NULL
    # zz <- seq(0, nrow(test_data), by = 1000)
    # zz <- lapply(2:length(zz), function(j){ c((zz[j-1]+1):zz[j])} )
    # for(j in 1:length(zz)) {
    #   extra <- c(extra, zz[[j]][wIf_bypart(Xtest = Xtest[zz[[j]], ], Xcal = Xcal)])
    # }
    exfoo <- function() {
      zz <- seq(0, nrow(Xtest), by = 1000)
      return(
        do.call(what = "c",
                args = lapply(lapply(2:length(zz), function(j){ c((zz[j-1]+1):zz[j])} ), 
                              function(x) { x[wIf_bypart(Xtest = Xtest[x, ], Xcal = Xcal)] }
                )
        )
      )
    }
    extra <- exfoo()
  }
  ### brute force
  else {
    extra <- whatif(formula = NULL, data = Xcal, cfact = Xtest)
    extra <- which(extra$in.hull == FALSE)
  }
  if(!index) { extra <- length(extra)/nrow(Xtest) }
  return(extra)
}