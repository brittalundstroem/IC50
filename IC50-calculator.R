################################################################################
#
# Calculate IC50 from imported data
# Author: Dominic Ritler 
# Version: 0.0.1
#
################################################################################

library(drc)
#library(ggplot2)
#library(ic50)

# for test
#library(PharmacoGx) 
#library(grid)
#library(gridExtra)

# Import data
 import.raw <- function(raw.path) {
   # Import raw data in condensed format
   #
   # Args:
   #  raw.path: The path to the file containing the data in the format
   #            col1: concentrations, col2: value
   #  
   #  Returns:
   #    data.frame of the raw data
   #
   #  TODO: include import tests and handling of different data formats 
   #
   
  # as alternative, one could use read.csv2(), I generally had good rexperience with the latter when using csv exported from MS-Excel,
  # I believe the decimal point is another issur, and that depends on the system language, especially GER/FR vs. ENG.
  # We might test, at a later stage i this is a source of bug.
   read.csv(raw.path, sep=";")
 }

# evaluate the best model function
  best_fit <- function(df, rank = 1){
   # function to get the best fitting model from a list
   # a first model is  fit using drm()
   # then the funciton mselect() is used to select the best model
   # by default the best fit model is return (rank = 1), we can select a second best using another rank
   #
   # the first part of the function is complicated and I don't know why it has to be like this
   # otherwise the mselect() function doesn't work within a function
   #
   # Args:
   #   df: dataframe with 2 col an n rows, col1: concentrations, col2: value
   #   rank: integer to select wich model to select. default = best model (rank = 1)
   #
   # Returns:
   #  Vector of parameters ordered by the best to the worst fit
   #
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("df"), names(mf), 0L)
   data.name <- as.character(mf[m])
   model1 <- eval(parse(text = paste0('drm(value ~ conz, data = ', data.name, ', fct = LL.4())')
   )
   )
   
   # compare different models
   M1 <- drc::mselect(model1, list(LL.3(), LL.5(), W1.3(), W1.4(), W2.3(), W2.4()
   )
   )
   # return model names according to their performance
   best <- row.names(M1)[rank]
   return(best)
 }
 
 
# run analysis 
ic50.calc <- function(df, bestmodel, ci = 0.95) {
  # calculate IC50 using drm
  # evaluate the best model returned with the best_fit() function
  # and builds the final formula for the best model
  #
  # Args:
  #  df: dataframe with 2 col an n rows, col1: concentrations, col2: value
  #  bestmodel: a character string corresponding to the best fitting function
  #  
  #  Returns:
  #    list with drm object, and a matrix of IC50 with CI
  #
  #  TODO: 
  #
  
  # define the parameter names
  arg_3 <- "(names = c('slope', 'high', 'IC50'))"
  arg_4 <- "(names = c('slope', 'low', 'high', 'IC50'))"
  arg_5 <- "(names = c('slope', 'low', 'high', 'IC50', 'f'))"
  
  # select the model function and parameter list based on the selected one
  if (bestmodel %in% c("LL.3", "W1.3", "W2.3")){
    funct <- paste0(bestmodel, arg_3)
  } else if (bestmodel %in% c("LL.4", "W1.4", "W2.4")){
    funct <- paste0(bestmodel, arg_4)
  } else if (bestmodel %in% c("LL.5")){
    funct <- paste0(bestmodel, arg_5)
  }
  else {
    stop("invalid model function selected")
  }
  
  # fit the model with corresponding model function and parameters
  fit <- drm(value ~ conz, data = df, 
             fct = eval(parse(text = funct)),
             type = "continuous")
  
  # Extract IC50 and calculate confidence interval
  fit.confint <- ED(fit, 50,
                    interval = "delta",
                    level = ci,
                    display = FALSE)
  
  # put together the model and IC in a list
  ret.val <- list(fit, fit.confint)

  return(ret.val)
  
}

# plot data
plot.ic50 <- function(fit, raw.d, nam="ic50") {
  # plot ic50 fit to pdf file
  #
  # Args:
  #  fit: list form ic50.calc()
  #  raw.d: data.frame from import.raw()
  #  nam: string (name for pdf file)
  #  
  #  Returns:
  #    None
  #
  #  TODO: adapt axes limit depending on CI min and max values
  #
  
  # get min, max for plot margins
  min.v <- min(raw.d$conz, fit[[2]][1, 3])
  max.v <- max(raw.d$conz, fit[[2]][1, 4])
  
  # plot
  pdf(paste(nam, ".pdf"))
  plot(fit[[1]], main = nam,
       xlim = c(min.v - 1, max.v + 1), 
       ylim = c(0, 100),
       xlab = "concentration")
  points(raw.d$conz, raw.d$value, col="blue", pch = "x")
  abline(v = fit[[2]][1, 1], col = "red")
  abline(v = fit[[2]][1, 3])
  abline(v = fit[[2]][1, 4])
  dev.off()
  
}

# export results
print.res <- function(fit.d, raw.d, CI=0.95, res.nam="result") {
  # export ic50 fit to results
  #
  # Args:
  #  fit: list from ic50.calc()
  #  raw.d: data.frame from import.raw()
  #  CI: Num (value btwen 0 and 1) def= 0.95
  #  res.nam: string (name for pdf file)
  #  
  #  Returns:
  #    None
  #
  #  TODO: 
  #
  
  sink(paste(res.nam, ".txt"))
  
  cat(paste("IC50 report: ", Sys.time(), " for: ", res.nam))
  cat("\n")
  cat(paste("R version: ", getRversion()), " and packages: ")
  print(names(sessionInfo()$otherPkgs))
  #a <- Sys.time()
  
  print(summary(fit.d[[1]]))
  
  cat("")
  cat("###################################################", "\n")
  cat(paste("IC50: ", round(fit.d[[2]][1, 1], 6)), "\n")
  cat(paste("IC50 CI l: ", round(fit.d[[2]][1, 3], 6), " for CI = ", CI), "\n")
  cat(paste("IC50 CI h: ", round(fit.d[[2]][1, 4], 6), " for CI = ", CI), "\n")
  cat("")
   
  if (fit.d[[2]][1, 1] < min(raw.d$conz)) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 smaller than lowest concentration")
    cat("###################################################")
  } 
  
  if (fit.d[[2]][1, 1] > max(raw.d$conz)) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 bigger than highest concentration")
    cat("###################################################")
  } 
  sink()
}


ic50calc <- function(f.path, nam="ic50", CI=0.95, rank = 1) {
  # main function
  #
  # Args:
  #  f.path: String (path to rawdata file)
  #  nam: string (name for pdf and result file)
  #  CI: Num (value btwen 0 and 1) def= 0.95  
  #
  #  Returns:
  #    None
  #
  #  TODO: 
  #
  
  # import data
  iraw.d <<- import.raw(f.path)
  
  # find the best fitting formula
  best.fit <- best_fit(iraw.d, rank)
  
  # make fit and ci
  fit.d <- ic50.calc(iraw.d, best.fit, CI)
  
  # plot fit
  plot.ic50(fit.d, iraw.d, nam)
  
  #export data
  print.res(fit.d, iraw.d, CI, nam)
}

################################################################################
# Debug and test 

# run script
f.path <- "Brittas_Samples/a1.csv"
nam <- "ic50"
CI <- 0.95

ic50calc(f.path, nam, CI, 1)
