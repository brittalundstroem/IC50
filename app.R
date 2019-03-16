# User interface ----

ui <- fluidPage(
  navbarPage(
    # App title ----
    title = 'IC 50 Calculation',
    
    # Tab 1 Upload Dataset ----
    tabPanel('Upload File',
             sidebarLayout(
               
               # Sidebar panel: Load inputs ----
               sidebarPanel(
                 
                 # Input: Select a file ----
                 fileInput("file1", "Choose CSV File",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 # Input: Select separator ----
                 radioButtons("sep", "Separator",
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
             includeMarkdown("about.md")
             )
    )
  )



# Server part # ----------------

# Define server logic to read selected file ----
server <- function(input, output) {
  
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


# The Shiny app function # ----------------

shinyApp(ui, server)