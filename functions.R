# functions before running Shiny App ----

install_load <- function (package1, ...)  {   
  # function to check whether packages are installed and load them
  # packages that are not installed will be installed and loaded
  # packages that are already installed will be loaded
  
    # args
    #  vector of packages names

  # convert arguments to vector
  packages <- c(package1, ...)

  # start loop to determine if each package is installed
  for(package in packages){

    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))

    # if package is not installed locally, download, then load
      else {
        install.packages(package)
        do.call("library", list(package))
      }
   }
}



# functions for the Shiny App ----

my_mselect <- function (object, fctList = NULL, nested = FALSE,
                        sorted = c("IC", "Res var", "Lack of fit", "no"), linreg = FALSE, icfct = AIC) {
  # this is an adapted version of the function drc::mselect
  # the difference is to add the source of the data in order to make the whole selection work
  sorted <- match.arg(sorted)
  if (!is.logical(nested)) {
    stop("'nested' argument takes only the values: FALSE, TRUE")
  }
  contData <- identical(object$type, "continuous")
  nestedInd <- 3 + contData + nested
  mc <- match.call()
  lenFL <- length(fctList)
  retMat <- matrix(0, lenFL + 1, 3 + contData + nested)
  retMat[1, 1] <- logLik(object)
  retMat[1, 2] <- icfct(object)
  retMat[1, 3] <- modelFit(object)[2, 5]
  if (contData) {
    tryRV <- try(summary(object)$resVar, silent = TRUE)
    if (!inherits("tryRV", "try-error")) {
      retMat[1, 4] <- tryRV
    }
    else {
      retMat[1, 4] <- NA
    }
  }
  if (nested) {
    retMat[1, nestedInd] <- NA
  }
  fctList2 <- rep("", lenFL + 1)
  fctList2[1] <- object$fct$name
  if (!is.null(fctList)) {
    prevObj <- object
    for (i in 1:lenFL) {
      tempObj <- try(update(object, fct = fctList[[i]],
                            data = object$origData), # <--- line added
                     silent = TRUE)
      fctList2[i + 1] <- fctList[[i]]$name
      if (!inherits(tempObj, "try-error")) {
        retMat[i + 1, 1] <- logLik(tempObj)
        retMat[i + 1, 2] <- icfct(tempObj)
        retMat[i + 1, 3] <- modelFit(tempObj)[2, 5]
        if (contData) {
          tryRV2 <- try(summary(tempObj)$resVar, silent = TRUE)
          if (!inherits("tryRV2", "try-error")) {
            retMat[i + 1, 4] <- tryRV2
          }
          else {
            retMat[i + 1, 4] <- NA
          }
        }
        if (nested) {
          retMat[i + 1, nestedInd] <- anova(prevObj, 
                                            tempObj, details = FALSE)[2, 5]
        }
      }
      else {
        retMat[i + 1, ] <- NA
      }
      prevObj <- tempObj
    }
  }
  rownames(retMat) <- as.vector(unlist(fctList2))
  cnames <- c("logLik", "IC", "Lack of fit")
  if (contData) {
    cnames <- c(cnames, "Res var")
  }
  if (nested) {
    cnames <- c(cnames, "Nested F test")
  }
  colnames(retMat) <- cnames
  if (linreg) {
    drcData <- as.data.frame(object$data[, c(2, 1)])
    names(drcData) <- c("yVec", "xVec")
    linFitList <- list(lm(yVec ~ xVec, data = drcData), lm(yVec ~ 
                                                             xVec + I(xVec * xVec), data = drcData), lm(yVec ~ 
                                                                                                          xVec + I(xVec * xVec) + I(xVec * xVec * xVec), data = drcData))
    linModMat <- matrix(unlist(lapply(linFitList, function(listObj) {
      c(logLik(listObj), icfct(listObj), NA, (summary(listObj)$sigma)^2)
    })), 3, 4, byrow = TRUE)
    rownames(linModMat) <- c("Lin", "Quad", "Cubic")
    colnames(linModMat) <- cnames[1:4]
    if (nested) {
      retMat <- retMat[, 1:4]
    }
    retMat <- rbind(retMat, linModMat)
  }
  if (sorted != "no") {
    return(retMat[order(retMat[, sorted]), ])
  }
  else {
    return(retMat)
  }
}


first_function <- function(df){
  # This is the fist drm function to be calculated.
  # Tests whether a model can be fit
  # if yes, it returns a first model which can then be compared
  # to similar models but using other function
  #
  # Args: 
  #  dataframe
  #
  # Output:
  #  a model of class 'drm'
  #
  
  # Check if the function can be fitted
  fcheck <- try(drm(value ~ conz,
                    data = df,
                    fct = LL.4()),
                silent = TRUE)
  if (inherits(fcheck, 'try-error')) {
    # abort if no function can be fitted
    stop("No function could be fitted please check the dataset")
  }

    
  # build the model
  model1 <- drm(value ~ conz,
                data = df,
                fct = LL.4())
  
  return(model1)
}

best_model <- function(model_ranks, df, rank = 1) {
  # fits the "best" model according to the choice of the user
  # The "best" is chosen by the function mselect() below in the server() part of the App
  # There might be clearer way to run it but this is the way I could make it work
  #
  # Args:
  #  model_ranks: an object of class 'matrix' as output from the mselect() function
  #  df: the input dataframe
  #  rank (default 1): the rank of the model (according to AIC) chosen by the user.
  #
  # Output:
  #  a model of class 'drm'
  
  # define the parameter names
  arg_3 <- "(names = c('slope', 'high', 'IC50'))"
  arg_4 <- "(names = c('slope', 'low', 'high', 'IC50'))"
  arg_5 <- "(names = c('slope', 'low', 'high', 'IC50', 'f'))"
  
  # return model name according to the choosen rank
  best <- row.names(model_ranks)[rank]
  
  # select the model function and parameter list based on the selected one
  if (best %in% c("LL.3", "W1.3", "W2.3")){
    funct <- paste0(best, arg_3)
  } else if (best %in% c("LL.4", "W1.4", "W2.4")){
    funct <- paste0(best, arg_4)
  } else if (best %in% c("LL.5")){
    funct <- paste0(best, arg_5)
  } else {
    stop("invalid model function selected")
  }
  
  # fit the model with corresponding model function and parameters
  fit <- drm(value ~ conz,
             data = df,
             fct = eval(parse(text = funct)),
             type = "continuous")
  return(fit)
}


rankTable <- function(mat){
  # function to reorder the output matrix of mselect and print a table of models
  # from best to worst
  # Only if a model could be fitted or not
  #
  # Args:
  #   mat: an object of class 'matrix' to reorder the model names
  #
  # Output
  #  a text table of model ranks
  

  # matrix of names to be finally sorted
  mod_struct <- data.frame(
    cbind(c("LL.3", "LL.4", "LL.5", "W1.3", "W1.4", "W2.3", "W2.4"),
          c("3 param. Logistic", "4 param. Logistic", "5 param. Logistic",
            "3 param. Weibull typ 1", "4 Param. Weibull typ 1",
            "3 param. Weibull typ 2", "4 Param. Weibull typ 2")
          )
    )
  abbrsort <- as.vector(rownames(mat))
  mod_struct$X1 <- factor(mod_struct$X1,
                          levels = abbrsort)
  tmp <- mod_struct[order(mod_struct$X1), ]
  #modelSorted <- data.frame("model_list" = tmp[,2],
  #                          row.names = 1:7)
  modelSorted <- matrix(tmp[,2],
                        nrow = 7, ncol = 1,
                        dimnames = list(c(1:7),
                                        c("ranked model list")))
  print(modelSorted)
  #print(as.matrix(modelSorted))
  
}


best_model_CI <- function(model, ci = 0.95){
  # function to extract IC50 values from the 'best' model
  # and calculate confidence interval
  #
  # Args:
  #  model: a model of class 'drm'
  #  ci: level of alpha error chosen by the user
  
  # calculate the CI of ic50
  fit.confint <- ED(model, 50,
                    interval = "delta",
                    level = ci,
                    display = FALSE)
  return(fit.confint)
}

# function: plot data ----
plot.ic50 <- function(fit, CImat, nam="ic50") {
  # plot ic50 fit
  #
  # Args:
  #  fit: list form ic50.calc()
  #  raw.d: data.frame from import.raw()
  #  nam: string (name for pdf file)
  
  #  Returns:
  #    None
  
  # get min, max for plot margins
  min.v <- min(fit$origData[, 1], CImat[1, 3], na.rm = TRUE)
  max.v <- max(fit$origData[, 1], CImat[1, 4], na.rm = TRUE)
  min.h <- min(fit$origData[, 2], na.rm = TRUE)
  max.h <- max(fit$origData[, 2], na.rm = TRUE)
  
  # set figure margins
  par(mar = c(5, 5, 3, 8) + 0.1) # set figure margins
  
  # create the plot
  plot(fit, main = nam,
       xlim = c(min.v - 1, max.v + 1), 
       ylim = c(min.h - 1, max.h + 1),
       xlab = "concentration",
       ylab = "inhibition")
  points(fit$origData[, 1], fit$origData[, 2], col="blue", pch = "x")
  abline(v = CImat[1, 1], col = "red")
  abline(v = CImat[1, 3])
  abline(v = CImat[1, 4])
  
  # add figure legend topright, outside plot region
  legend("topleft",
         inset = c(1, 0),
         xpd = TRUE,
         legend = c("raw data","fitted values", "fitted curve", "IC50", "IC50 CI"),
         col = c("blue", "black", "black", "red", "black"),
         lty = c(NA, NA, 1 , NA, NA),
         pch = c("x", "o", NA, "|" , "|"),
         bty = "n")
}


# export results
print.res <- function(fit.d, CImat, raw.d, CI=0.95, res.nam="result") {
  # function to build the model summary and the CI
  
  # Args:
  #  fit: list from ic50.calc()
  #  raw.d: data.frame from import.raw()
  #  CI: Num (value btwen 0 and 1) def= 0.95
  #  res.nam: string (name for pdf file)
  
  #  Returns:
  #    a (long) character string
  
  #  TODO: 
  
  #sink(paste(res.nam, ".txt"))
  
  cat(paste("IC50 report: ", Sys.time(), " for: ", res.nam))
  cat("\n")
  cat(paste("R version: ", getRversion()), " and packages: ")
  print(names(sessionInfo()$otherPkgs))

  print(summary(fit.d))
  
  cat("")
  cat("###################################################", "\n")
  cat(paste("IC50:     ", round(CImat[1, 1], 4)), "\n")
  cat(paste("IC50 LCL: ", round(CImat[1, 3], 4), " for CI =", CI), "\n")
  cat(paste("IC50 UCL: ", round(CImat[1, 4], 4), " for CI =", CI), "\n")
  cat("")
  
  if (CImat[1, 1] < min(raw.d[, 1])) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 smaller than lowest concentration")
    cat("###################################################")
  } 
  
  if (CImat[1, 1] > max(raw.d[, 1])) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 bigger than highest concentration")
    cat("###################################################")
  } 
  #sink()
}
