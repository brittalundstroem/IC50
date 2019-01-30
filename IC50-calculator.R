################################################################################
#
# Calculate IC50 from imported data
# Author: Dominic Ritler 
# Version: 0.0.1
#
################################################################################

library(drc)
library(ggplot2)
library(ic50)

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

# run analysis 
ic50.calc <- function(df, ci = 0.95) {
  # calculate IC50 using drm
  #
  # Args:
  #  df: dataframe with 2 col an n rows, col1: concentrations, col2: value
  #  
  #  Returns:
  #    list with drm object, confint, low and high
  #
  #  TODO: 
  #
  
  fit <- drm(value~conz, data=df, 
               # currently, the change from LL.4 to LL.3 is my major change, which affects the rest of the code
             # since a concentration of 0% of any compound is expected to have 0% inhibition,
             # I believe a 3 parameter model (forcing lower assymptote to 0) is adequate.
             # Indeed, using LL.4, the test dataset gave a large CI for lower asymptote which included absurd negative values.
             # One might even argue, wheter a 2 parameter model would be adequate, fixing maximum value at 100,
             # given that a large concentration of any compound would lead to 100% inhibition.
             # The latter is rather theoretical (toxicology dogma: "Any compound is toxic, it's the dosage that's important").
             # However I'm not sure its of practical use, therefore I thought the 3 parameter model would better suit the needs.
             fct=LL.3(names=c("slope","high","IC50")),
             type="continuous")
  
  fit.confint <- confint(fit, level = ci)
  fit.conf.low <- fit.confint[,1] # conf.low
  fit.conf.hig <- fit.confint[,2] # conf.high
  
  ret.val <- list(fit, fit.confint, fit.conf.low, fit.conf.hig)

  # test and debug
  #summary(fit.d)
  
  
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
  #  TODO: 
  #
  
  # get min, max for plot margins
  min.v <- min(raw.d$conz, fit[[3]][3]) # we might change raw.d$conz to raw.d[, 1] since it depends on the table, if someone changes the column name say to Konz, the formula won't work
  max.v <- max(raw.d$conz, fit[[4]][3]) # same as line 93
  
  # plot
  pdf(paste(nam, ".pdf"))
  plot(fit[[1]], main=nam, xlim = c(min.v - 1, max.v + 1), 
       ylim = c(0, 100),  # not sure if you like it, I like to have an overview of the curve. With the test dataset, we see for instance, that more lower concentrations would have helped in better defining the curve.
       xlab = "concentration")
  #plot(fit[[1]], type = "confidence", confidence.level = 0.95, add = TRUE)
  points(raw.d$conz, raw.d$value, col="blue", pch = "x")
  abline(v=fit[[1]]$parmMat[3], col = "red")
  abline(v=fit[[3]][3])
  abline(v=fit[[4]][3])
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
  cat("###################################################")
  cat(paste("IC50: ", round(fit.d[[1]]$parmMat[3], 6)))
  cat(paste("IC50 CI l: ", round(fit.d[[3]][3], 6), " for CI = ", CI))
  cat(paste("IC50 CI h: ", round(fit.d[[4]][3], 6), " for CI = ", CI))
  cat("")
   
  if (fit.d[[1]]$parmMat[3] < min(raw.d$conz)) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 smaller than lowest concentration")
    cat("###################################################")
  } 
  
  if (fit.d[[1]]$parmMat[3] > max(raw.d$conz)) {
    cat("###################################################")
    cat("Extrapolation: ")
    cat("IC50 bigger than highest concentration")
    cat("###################################################")
  } 
  sink()
}


ic50calc <- function(f.path, nam="ic50", CI=0.95) {
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
  iraw.d <- import.raw(f.path)
  
  # make fit and ci
  fit.d <- ic50.calc(iraw.d, CI)
  
  # plot fit
  plot.ic50(fit.d, iraw.d, nam)
  
  #export data
  print.res(fit.d, iraw.d, CI, nam)
}

################################################################################
# Debug and test 

# run script
f.path <- "Brittas_Samples/a1.csv"
nam="ic50"
CI <- 0.95

ic50calc(f.path, nam, CI)
