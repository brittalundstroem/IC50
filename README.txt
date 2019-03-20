IC50 calculation tool

In order to work locally, you need at least base R installation and an internet browser.
For the application R also needs packages "shiny", "drc", and "zip".
However, the script automatically checks whether these packages are installed.
If not, it installs the missing packages and loads them automatically.
You might be asked to select a CRAN mirror to download the packages binaries.


In order to use the tool, you need to run following commands in your console:

setwd("[path to your folder]")
source("IC50calculator.R")


If you are using R studio, open the script and click on 'Source'.