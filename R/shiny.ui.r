

#  shinyUI(bootstrapPage())

shiny.ui <- fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      
      ##
      ## Data input
      ##
      fileInput('file1', 'Data (choose CSV file)',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
      ),
      
      ##
      ## Main controls
      ##
      hr(),
      selectInput(inputId = "method",
                  label = "Method of Isotonization:",
                  choices=c("QPE","PAVA") 
      ),
      selectInput(inputId = "tsfunc",
                  label = "Test Statistic:",
                  choices=c("LRT", "Williams") 
      ),
      checkboxGroupInput(inputId = "order",
                         label = "Order (select at least one):",
                         choices=c("Simple","Umbrella","Tree") 
      ),
      numericInput(inputId = "p1",
                   label = "Number of Groups (p1):",
                   min=2 , max=50 , value=2
      ),
      numericInput(inputId = "p2",
                   label = "Number of Covariates (p2):",
                   min=0 , max=500 , value=0
      ),
      numericInput(inputId = "q",
                   label = "Number of Random Effects (Q):",
                   min=0 , max=500 , value=0
      ),
      checkboxInput(inputId = "decreasing",
                    label = "Decreasing order (inverted):",
                    value=FALSE
      ),
      numericInput(inputId = "mboot",
                   label = "Number of Bootstraps:",
                   min=0 , max=50000 , value=1000
      ),
      
      ##
      ## Action buttons
      ##
      hr(),
      helpText("Click 'compute' to run model with selected inputs."),
      actionButton(inputId = "compute1",
                   label = "Run model"
      ),
      actionButton(inputId = "compute2",
                   label = "Update Plot"
      ),
      actionButton(inputId = "compute3",
                   label = "Update Summary"
      ),
      
      
      ##
      ## Output controls
      ##
      hr(),
      checkboxInput(inputId = "outcheck",
                    label = "Format Output:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.outcheck",
                       checkboxInput(inputId = "plotci",
                                     label = "CI on Plot:",
                                     value=FALSE
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "alpha",
                                   label = "Alpha level:",
                                   min=0 , max=0.15 , value=0.05 , step=0.01
                       )
      ),
      br(),
      conditionalPanel(condition = "input.outcheck",
                       helpText("Number of decimal places for:")
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "decp",
                                   label = "p-values:",
                                   min=0 , max=8 , value=4 , step=1
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "dect",
                                   label = "Theta coefficients:",
                                   min=0 , max=8 , value=2 , step=1
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "decw",
                                   label = "Test statistics:",
                                   min=0 , max=8 , value=3 , step=1
                       )
      ),
      conditionalPanel(condition = "input.outcheck",
                       sliderInput(inputId = "decv",
                                   label = "Variance estimates:",
                                   min=0 , max=8 , value=4 , step=1
                       )
      ),
      
      ##
      ## Extra parameters
      ##
      hr(),
      helpText("Other Parameters"),      
      checkboxInput(inputId = "varSSQ",
                    label = "Heterogeneity for Residuals:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.varSSQ",
                       textInput(inputId = "Nks",
                                 label = "Group N's (format: n1,n2,n3,...):" )
      ),
      checkboxInput(inputId = "varTSQ",
                    label = "Heterogeneity for Random Effects:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.varTSQ",
                       textInput(inputId = "Qs",
                                 label = "Group Q's (format: q1,q2,q3,...):" )
      ),
      
      
      
      ##
      ## Technical controls
      ##
      checkboxInput(inputId = "technical",
                    label = "Specify Technical Parameters:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "em.iter",
                                    label = "Max EM Iterations:",
                                    min=10 , max=50000 , value=500
                       )
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mq.iter",
                                    label = "Max MINQUE Iterations:",
                                    min=10 , max=50000 , value=500)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "em.eps",
                                    label = "EM Convergence Criteria:",
                                    min=0 , max=50000 , value=0.00001)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mq.eps",
                                    label = "MINQUE Convergence Criteria:",
                                    min=0 , max=50000 , value=0.00001)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "ranseed",
                                    label = "Set RNG Seed:",
                                    min=0 , max=Inf , value=42)
      )
    ),
    mainPanel(
      plotOutput(outputId = "fig1", height = "400px"),
      verbatimTextOutput(outputId = "summary" )
      
    )
  ))












