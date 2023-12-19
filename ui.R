# TANDANG, Bernard Jezua R.
# 2021-09992
# CMSC 150 - B1L

# This is the UI code for the CMSC 150 final project. 
# You can run the application by clicking the 'Run App' button.

library(readr)
library(shiny)
library(shinyjs)
library(shinyMatrix)
library(shinythemes)
source("QuadraticSpline.R")
source("PolynomialReg.R")
# source("SimplexMethod.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Google+Sans:wght@400;700&display=swap"
    ),
    tags$style(
      HTML(
        "
      body { font-family: 'Google Sans', sans-serif; }
      .navbar { background-color: #3498db; color: white; }
      .navbar-default .navbar-nav > .active > a, 
      .navbar-default .navbar-nav > .active > a:hover, 
      .navbar-default .navbar-nav > li > a:hover, 
      .navbar-default .navbar-nav > li > a:focus { background-color: #2980b9; color: white; }
      .navbar-default .navbar-brand { color: white; }
      .navbar-default .navbar-brand:hover { color: white; }
      .tab-content { background-color: #ecf0f1; }
      h1 { font-size: 36px; font-weight: bold; color: #3498db; margin-bottom: 20px; }
      h2 { font-size: 28px; font-weight: bold; color: #3498db; margin-bottom: 15px; }
      h4 { font-size: 20px; font-weight: bold; color: #3498db; margin-bottom: 10px; }
      p { font-size: 16px; color: #555; }
      .btn-primary { background-color: #3498db; color: white; border: none; }
      .btn-primary:hover { background-color: #2980b9; color: white; }
      "
      )
    ) # CSS Styles for the UI
  ),
  shinyjs::useShinyjs(), # Used for some styles
  h1("CMSC 150 Project", align = "center"),
  navbarPage(
    "",
    tabPanel(
      "QSI Calculator",
      wellPanel(
        h2("Quadratic Spline Interpolation Calculator", align = "center"),
      ),
      sidebarPanel(
        h2("Parameters"),
        numericInput("Value", "Enter Value to be Evaluated", 0),
        fileInput("file", "Upload CSV File"),
        p("Note: CSV file must have two columns for x and y values. Otherwise, it will print as nothing or NULL values."),
        actionButton("Generate", label = "Generate", class = "btn-primary")
      ),  # sidebar panel
      mainPanel(
        h4("Function Per Interval"),
        p("Generated Functions Per Interval"),
        verbatimTextOutput("qsi.fxns"),
        h4("Estimated Value"),
        p("Estimated Value of Target"),
        verbatimTextOutput("est"),
      )  # main panel
    ), # tab panel qsi
    tabPanel(
      "PR Calculator",
      wellPanel(
        h2("Polynomial Regression Calculator", align = "center"),
      ),
      sidebarPanel(
        h2("Parameters"),
        numericInput("polyDegree", "Polynomial Degree", 0),
        numericInput("polyXVal", "Enter X Value", 0),
        fileInput("polyFile", "Upload CSV File"),
        p("Note: CSV file must have two columns for x and y values. Otherwise, it will print as nothing or NULL values."),
        actionButton("RunPolyReg", label = "Generate", class = "btn-primary")
      ),
      mainPanel(
        h4("Polynomial Function"),
        p("Generated Polynomial Function"),
        verbatimTextOutput("polyFunction"),
        
        h4("Estimated Value"),
        p("Estimated Value of Target"),
        verbatimTextOutput("polyEstimate"),
      )
    ), #tab panel pr
    tabPanel(
      "Simplex Method",
      wellPanel(
        h2("Simplex For Diet Problem", align = "center"),
      ),
      sidebarPanel(
        h2("Parameters"),
        # Add UI elements for Simplex Method parameters
        selectInput("foodItems", "Select Food Items", choices = NULL, multiple = TRUE),
        actionButton("checkAll", "Check All"),
        actionButton("resetAll", "Reset All"),
        actionButton("solveSimplex", label = "Solve Simplex", class = "btn-primary")
      ),
      mainPanel(
        # Add outputs for the Simplex Method
        h4("Optimal Diet Plan"),
        verbatimTextOutput("simplexIterations"),
        h4("Simplex Iterations"),
        tableOutput("optimalDiet")
      )
    )
  )
)

# Server function for all actions
server <- function(input, output, session) {
  # Create reactive values to store matrix data, set as NULL to avoid bugs
  matrixDataQSI <- reactiveVal(NULL)
  matrixDataPoly <- reactiveVal(NULL)
  
  # ================== QUADRATIC SPLINE ==================
  # QSI Reactive expression when "Generate" button is clicked
  quadraticSInterpolation <- eventReactive(input$Generate, {
    req(matrixDataQSI()) # checks if matrixDataQSI has a value
    x <- matrixDataQSI()[, 1] # first column
    y <- matrixDataQSI()[, 2] # second column
    
    if (length(x) != length(y) || any(is.na(x)) || any(is.na(y))) {
      stop("Data for QSI must have exactly two columns with matching non-missing rows.")
    }
    
    poly.qsi(list(x, y), input$Value) # calls out poly.qsi() function in QuadraticSpline.R
  })
  
  # Updates uploaded CSV data for Quadratic Spline Interpolation
  observeEvent(input$file, {
    req(input$file) # checks if CSV file has been uploaded
    # Try-Catch method to handle error and avoid app crash
    tryCatch({
      data <- read.csv(input$file$datapath, header = FALSE, stringsAsFactors = FALSE)
      # Checks for number of cols and data points if there are any missing values
      if (ncol(data) == 2 && !any(is.na(data[[1]])) && !any(is.na(data[[2]]))) { 
        matrixDataQSI(data)
      } else {
        matrixDataQSI(NULL)  # Set matrix data to NULL to avoid using incomplete data
        stop("CSV file must have exactly two columns and no missing values.")
      }
    }, error = function(e) {
      warning(paste("Error in reading CSV file:", e))
    })
  })
  
  # ================== POLYNOMIAL REGRESSION ==================
  # PR Reactive expression when "Generate" button is clicked
  polynomialRegression <- eventReactive(input$RunPolyReg, {
    req(matrixDataPoly(), input$polyFile)
    data <- read.csv(input$polyFile$datapath, header = FALSE, stringsAsFactors = FALSE)
    
    if (ncol(data) != 2) {
      return(list(polyFunction = NA, polyEstimate = NA))
    } else if (any(is.na(data))) {
      return(list(polyFunction = NA, polyEstimate = NA))
    }
    
    if (input$polyDegree < 1) {
      return(list(polyFunction = NA, polyEstimate = NA))
    }
    
    PolynomialRegression(input$polyDegree, list(data[[1]], data[[2]]), input$polyXVal)
  })
  
  # Updates uploaded CSV data for Polynomial Regression
  observeEvent(input$polyFile, {
    req(input$polyFile)
    # Try-Catch method to handle error and avoid app crash
    tryCatch({
      data <- read.csv(input$polyFile$datapath, header = FALSE, stringsAsFactors = FALSE)
      # Checks for number of cols and data points if there are any missing values
      if (ncol(data) == 2 && !any(is.na(data[[1]])) && !any(is.na(data[[2]]))) { 
        matrixDataPoly(data)
      } else {
        matrixDataPoly(NULL)  # Set matrix data to NULL to avoid using incomplete data
        stop("CSV file must have exactly two columns and no missing values.")
      }
    }, error = function(e) {
      warning(paste("Error in reading CSV file:", e))
    })
  })
  
  # ================== SIMPLEX METHOD ==================
  foods_data <- read.csv("FoodItem.csv") # Read the foods.csv file
  food_items <- foods_data$Foods # Extract food items from column 1 (excluding the header)
  food_items <- food_items[-c((length(food_items) - 1):length(food_items))]
  
  # Create reactiveValues to store original and selected food items
  foodItemsData <- reactiveValues(
    original = foods_data,  # Store the original data
    selected = NULL
  )
  
  # Update choices in selectInput
  updateSelectInput(session, "foodItems", choices = food_items)
  
  # Check all button
  observeEvent(input$checkAll, {
    updateSelectInput(session, "foodItems", selected = food_items)
  })
  
  # Reset all button
  observeEvent(input$resetAll, {
    updateSelectInput(session, "foodItems", choices = character(0), selected = NULL) # deletes all selected
    updateSelectInput(session, "foodItems", choices = food_items) # updates food items from csv
  })
  
  # Create a reactiveValues to store the Solve Simplex button state
  solveSimplexClicked <- reactiveValues(clicked = FALSE)
  
  # Reactive expression to determine whether to update optimalDiet
  updateOptimalDiet <- reactive({
    input$foodItems  # Include input$foodItems as a dependency
    
    if (solveSimplexClicked$clicked) {
      # Check if any food items are selected
      if (is.null(input$foodItems) || length(input$foodItems) == 0) {
        return(data.frame(Message = "No food items selected. Please check at least one item."))
      }
      
      # Print the selected values
      selected_values <- numeric()
      
      # Extract numeric values for selected food items
      for (food in input$foodItems) {
        selected_values <- c(selected_values, as.numeric(foods_data[foods_data$Foods == food, -1]))
      }
      
      # Changes it back to FALSE in order to be pressed again
      solveSimplexClicked <- reactiveValues(clicked = FALSE)
      
      # Display the selected values in the table
      return(data.frame(Selected_Values = selected_values))
    } else {
      # Solve Simplex button not clicked, return NULL
      return(NULL)
    }
  })
  
  # Output for optimalDiet
  output$optimalDiet <- renderTable({
    updateOptimalDiet()
  })
  
  # Observer for Solve Simplex button
  observeEvent(input$solveSimplex, {
    solveSimplexClicked$clicked <- TRUE
  })
  
  # Outputs for QSI
  output$qsi.fxns = renderPrint({ quadraticSInterpolation()$qsi.fxns })
  output$est = renderText({ quadraticSInterpolation()$y })
  
  # Outputs for polynomial regression
  output$polyFunction <- renderPrint({ polynomialRegression()$polynomial_string })
  output$polyEstimate <- renderText({ polynomialRegression()$estimate })
}

# Run the application 
shinyApp(ui = ui, server = server)
