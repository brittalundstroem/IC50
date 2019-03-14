





install.r.packages <- function(r.install.names) {
  # install r packages
  #
  # Args:
  #   install.names: (string) names of the R packages to be installed from CRAN
  #
  # Return:
  #   none (install packages if not installed)
  
  install.packages(r.install.names)
  
}


test.installed.r.packages <- function(r.package.names) {
  # test if r package are installed and if not install it.
  #
  # Args:
  #   package.names: (vector of string) packages to be needed for the analysis
  #
  # Return:
  #   none (install packages if not installed)
  
  # find not installed packages
  installed <- r.package.names %in% rownames(installed.packages())
  
  # install package if not installed but first test if list is empty
  if (length(r.package.names[!installed]) != 0) {
    install.r.packages(r.package.names[!installed])
  } else {
    #print("all R packages installed")
  }
}


# packages needed
pack.ned <- c("drc", "MASS", "markdown")
test.installed.r.packages("drc")



#


