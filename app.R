# UI part
########


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
                       selected = ";"),
          # Input: Select quotes ----
          radioButtons("quote", "Quote",
                       choices = c(None = "",
                                   "Double Quote" = '"',
                                   "Single Quote" = "'"),
                       selected = '"'),
          
          # Horizontal line ----
          tags$hr(),
          
          # Input: Select number of rows to display ----
          radioButtons("disp", "Display",
                       choices = c(Head = "head",
                                   All = "all"),
                       selected = "head")
          ),
        
        # Main Panel: Display inputs ----
        mainPanel(
          # Output: Data file ----
          tableOutput("dataset")
          )
        )
      ),
    
    # Tab 2 Generate IC 50 ----
    tabPanel('Generate IC50',
      sidebarLayout(
 
        # Sidebar panel: Select options ----
        sidebarPanel(),

        # Main panel: display outputs ----
        mainPanel()
        )
      )
    )
  )



########
# Server part
########


# Define server logic to read selected file ----
server <- function(input, output) {

    output$dataset <- renderTable({
    
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
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
}



########
# The Shiny app function
########

shinyApp(ui, server)