################################################################
# The commands parts
# In this part, the different functions necessary for the IC 50 calculator
# are defined.

# functions to use before launching the Shiny App ----
install_load <- function(pkg){
  # function to check whether packages are installed and load them
  # packages that are not installed will be installed and loaded
  # packages that are already installed will be loaded

  # args
    #  vector of packages names
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg)
  sapply(pkg, require, character.only = TRUE)
}

# Check if necessary packages are installed and load them
install_load(c("shiny", "drc", "zip"))


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

plot.ic50 <- function(fit, CImat, nam="ic50") {
  # function to draw the model function, the IC50 and confidence interval
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



################################################################
# The actual Shiny App part
# App is embedded within the runApp command
# this allows keeping a single file for the IC50 calculator
#
# as shown in the help document, the "UI" and "server" arguments
# are entered in a list.


runApp(list(
  # User interface part----
  ui = fluidPage(
    navbarPage(
      # App title ----
      title = 'IC 50 Calculation',
      
      # Tab 1 Upload Dataset ----
      tabPanel('Upload File',
               sidebarLayout(
                 
                 # Sidebar panel: Load inputs ----
                 sidebarPanel(
                   
                   # heading ----
                   h4("Upload your dataset"),
                   p("Use the 'Browse' function below to load your dataset in .csv format:"),
                   tags$ul(
                     tags$li("The first column must contain the concentration data,"),
                     tags$li("The second column must contain the inhibition values.")
                     ),

                   # Input: Select a file ----
                   fileInput("file1", NULL,
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                   
                   # Horizontal line ----
                   tags$hr(),

                   # heading ----
                   h4("Table Header"),
                   tags$p("If your table has headers, leave the checkmark below.
                         If your data table has no headers, remove the checkmark"),
                   
                   # Input: Checkbox if file has header ----
                   checkboxInput("header", "Header", TRUE),
                   
                   # Horizontal line ----
                   tags$hr(),
                   
                   # heading ----
                   h4("Table separator"),
                   tags$p("csv data exported from MS-Excel usually have semicolons (;) as separators.
                          If your data uses another separator, select the correct one below."),
                   
                   # Input: Select separator ----
                   radioButtons("sep", NULL,
                                choices = c(Comma = ",",
                                            Semicolon = ";",
                                            Tab = "\t"),
                                selected = ";")
                 ),
                 
                 # Main Panel: Display inputs ----
                 mainPanel(
                   # Output: Data file ----
                   tableOutput("data")
                 )
               )
      ),
      
      # Tab 2 Generate IC 50 ----
      tabPanel('Generate IC50',
               sidebarLayout(
                 
                 # Sidebar panel: Select options ----
                 sidebarPanel(
                   
                   # Action: Check the dataset ---
                   actionButton("Suitability", label = "Check the dataset"),
                   
                   # line break ----
                   tags$p(),
                   
                   # Output: Model Ranks ---
                   verbatimTextOutput("RanksTable"),
                   
                   # line break ----
                   tags$hr(),
                   
                   # Input: Model Rank selection ---
                   selectInput("RankSelect", "Select model rank (1 should be best)", 
                               choices = list("1" = 1, "2" = 2, "3" = 3,
                                              "4" = 4, "5" = 5, "6" = 6, "7" = 7), 
                               selected = 1),
                   
                   # Input: Confidence interval ---
                   numericInput("CI", "Confidence level (%)", value = 95,
                                min = 50, max = 99, step = 1),
                   
                   # Horizontal line ----
                   tags$hr(),
                   
                   # Action: Run the analysis ---
                   actionButton("RunIC50", label = "(Re-)Caclulate IC50"),
                   
                   # Horizontal line ----
                   tags$hr(),
                   
                   uiOutput("buttons")
                   
                   
                 ),
                 
                 # Main panel: display outputs ----
                 mainPanel(
                   
                   # Output: Plot of the IC50 ---
                   plotOutput("IC50plot"),
                   
                   # Horizontal line ----
                   tags$hr(),
                   
                   # Output: Summary Output IC50 calculation ---
                   verbatimTextOutput("IC50summary")
                 )
               )
      ),
      
      # Tab 3 About ----
      tabPanel('About',
               tags$h1("An IC 50 Calculator"),
               tags$h4("Calculate IC50 and confidence interval."),
               tags$p("The data need to be imported as a .csv file.
                      The first column should contain concentration values
                      and the second column should contain the inhibition values.
                      Make sure your table follows this structure before you upload the file."),
               tags$ol(
                 tags$li("Go to the",
                         tags$code("Upload File"),
                         "tab and browse your files to upload your dataset."),
                 tags$li("Check the printed table to ensure
                         your dataset has the right structure.
                         It is possible to use different column names, but the",
                         tags$strong("order of the columns"),
                         "should remain the same.",
                         tags$ol(
                           tags$li("If your dataset has no header,
                                   remove the checkmark from the",
                                   tags$code("Heather"),
                                   "checkbox."),
                           tags$li("If the cell separator is different than a semicolon (;),
                                   choose the corresponding character.")
                           )
                         ),
                 tags$li("Go to the",
                         tags$code("Upload File"),
                         "tab",
                         tags$ol(
                           tags$li("Click on",
                                   tags$code("Check the dataset"),
                                   ".",
                                   tags$ol(
                                     tags$li("If a model can be fitted a list of models is printed
                                             below the button. First model on the list
                                             is the best according to the Akaike information criterion.
                                             The second model on the list is the second best and so on."),
                                     tags$li("If no model can be fitted,
                                             a message asking to check the dat aset is printed.")
                                     )
                                   ),
                           tags$li("Select which model you want to fit from the selection box."),
                           tags$li("Choose an",
                                   HTML("&alpha;"),
                                   "error level for generating the confidence interval, default is 95%.")
                           )
                         ),
                 tags$li("Click on",
                         tags$code("Re-)Calculate IC50"),
                         ".")
               ),
               tags$p("The model will be plotted with the raw data. A summary of the model will be printed below the plot."),
               tags$p("I strongly recommended to compare several models,
                      since the AIC is no more than an automated selection process.
                      After changes are made, the plot and ouptut can be updated
                      by clicking again on",
                      tags$code("(re-)calculate IC50"),
                      "."),
               tags$p("The visible plot and model output can be downloaded in a .zip archive if desired."),
               tags$br(),
               tags$p("2019, Dominic Ritler, Nelson Marreros, Britta Lundström-Stadelmann, Institute of Parasitology, Vetsuisse Faculty, University of Bern, Switzerland")
               )
      )
  ),
  ## Server part # ----------------
  server = function(input, output) {
    
    # Reactive function to load the dataset -----
    mydat <- reactive({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      req(input$file1)
      
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df <- read.csv(input$file1$datapath,
                         header = input$header,
                         sep = input$sep,
                         quote = input$quote)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      return(df)
    })
    
    # functiont to show the dataset -----
    output$data <- renderTable({
      mydat()
    })
    
    # functiont to get a cached dataset with desirable column names -----
    mydata <- reactive({
      tmp <- mydat()
      names(tmp) <- c("conz", "value")
      tmp
    })
    
    
    # functions to check feasability of model -----
    firstmodel <- eventReactive(input$Suitability, {
      first_function(mydata())
    })
    
    # functions to build rank list of models -----
    ModelRanks <- reactive({
      my_mselect(firstmodel(), list(LL.3(), LL.5(), W1.3(), W1.4(), W2.3(), W2.4()))
    })
    
    # functiont to generate the model ranks table and output-----
    myRanksTable <- reactive({
      rankTable(ModelRanks())
    })
    output$RanksTable <- renderPrint(myRanksTable())
    
    # functions to calculate the desired model and CI from dataset -----
    bestmodel <- eventReactive(input$RunIC50, {
      myrank <- as.numeric(input$RankSelect)
      best_model(ModelRanks(), mydata(), myrank)
    })
    
    bestmodelCI <- eventReactive(input$RunIC50, {
      best_model_CI(bestmodel(), ci = input$CI/100)
    })
    
    
    # functiont to generate the plot and output-----
    IC50plot_t <- reactive({
      plot.ic50(fit = bestmodel(), CImat = bestmodelCI() )
    })
    
    output$IC50plot <- renderPlot({
      IC50plot_t()
    })
    
    # functions to generate the textoutput -----
    mysummary <- eventReactive(input$RunIC50, {
      print.res(bestmodel(), bestmodelCI(), mydat(), CI = input$CI/100, res.nam = input$file1$name)
    })
    output$IC50summary <- renderPrint(mysummary() )
    
    
    # function to show downloadplots only if a model was selected
    output$buttons <- renderUI({
      if(is.null(bestmodel())){return()}
      list(
        # report download button
        downloadButton(outputId = "down", label = "Download the zip output")
      )
    })
    
    # Downloadable output in zip file ----
    output$down <- downloadHandler(
      filename = "IC50_export.zip",
      content = function(fname) {
        fs <- c("IC50.txt", "IC50.pdf")
        tmpdir <- tempdir()
        setwd(tempdir())
        sink("IC50.txt")
        print.res(bestmodel(), bestmodelCI(), mydat(), CI = input$CI/100, res.nam = input$file1$name)
        sink()
        pdf("IC50.pdf")
        plot.ic50(fit = bestmodel(), CImat = bestmodelCI() )
        dev.off()
        zip::zip(zipfile = fname, files = fs)
      },
      contentType = "application/zip"
    )
    
  }
),
launch.browser = TRUE)

