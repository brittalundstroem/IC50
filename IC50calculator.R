# load necessary functions
source("functions.R")

# Check if necessary packages are installed and load them
pack_ned <- c("shiny", "drc", "MASS", "markdown", "zip")
install_load(pack_ned)

# run the App
runApp(launch.browser = TRUE)