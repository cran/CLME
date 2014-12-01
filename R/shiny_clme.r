


shiny_clme <- function(){
  library("CLME")
  runApp(
    list(
      ui     = shinyUI_clme,
      server = shinyServer_clme
    )
  )
}




##############################################################################
##
## The user interface for the shiny app
##
##############################################################################
#  shinyUI(bootstrapPage())

shinyUI_clme <- fluidPage(
  
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
      selectInput(inputId = "solver",
                  label = "Type of Solver / fitting norm:",
                  choices=c("Least Squares (LS)", "Least Absolute Value (L1)",
                            "General LS (GLS)", "Asymmetrix LS" , 
                            "L1 approx", "Huber", "SILF", "Chebyshev",
                            "L-p Power Norm", "Quantile", "Poisson" ) 
      ),      
      selectInput(inputId = "tsfunc",
                  label = "Test Statistic:",
                  choices=c("LRT", "Williams") 
      ),
      checkboxGroupInput(inputId = "order",
                         label = "Order (select at least one):",
                         choices=c("Simple","Umbrella","Tree")
      ),
      helpText("Identify columns of data"),
      helpText("Use column letters or numbers"),
      helpText("e.g., 1-3 or A-C or a list: 1,4,6"),
      textInput(inputId = "yy",
                   label = "Column of response:"),
      textInput(inputId = "p1",
                   label = "Column of constrained effect:"),
      textInput(inputId = "p2",
                   label = "Column(s) of Covariates:", value="None"),
      textInput(inputId = "q",
                   label = "Column(s) of random effects:", value="None"),
      checkboxInput(inputId = "decreasing",
                    label = "Decreasing order (inverted):",
                    value=FALSE
      ),
      numericInput(inputId = "nsim",
                   label = "Number of Bootstraps:",
                   min=0 , max=50000 , value=1000
      ),
      
      ##
      ## Action buttons
      ##
      hr(),
      helpText("Click to run model or update output."),
      actionButton(inputId = "compute1",
                   label = "Run model"
      ),
      actionButton(inputId = "compute2",
                   label = "Update plot"
      ),
      actionButton(inputId = "compute3",
                   label = "Update summary"
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
      #br(),
      #conditionalPanel(condition = "input.outcheck",
      #                 helpText("Number of decimal places for:")
      #),
      #conditionalPanel(condition = "input.outcheck",
      #                 sliderInput(inputId = "digits",
      #                             label = "p-values:",
      #                             min=0 , max=8 , value=4 , step=1
      #                 )
      #),

      ##
      ## Extra parameters
      ##
      hr(),
      checkboxInput(inputId = "varssq",
                    label = "Heteroscedasticity:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.varssq",
                       textInput(inputId = "gfix",
                                    label = "Column of variance groups:")
      ),
      checkboxInput(inputId = "xlevel1",
                    label = "Define order of constrained groups:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.xlevel1",
                       textInput(inputId = "xlevels",
                                 label = "Column of ordered group levels:")
      ),
      ##
      ## Technical controls
      ##
      checkboxInput(inputId = "technical",
                    label = "Select Control Parameters:",
                    value=FALSE
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "emiter",
                                    label = "Max EM Iterations:",
                                    min=10 , max=50000 , value=500
                       )
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mqiter",
                                    label = "Max MINQUE Iterations:",
                                    min=10 , max=50000 , value=500)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "emeps",
                                    label = "EM Convergence Criteria:",
                                    min=0 , max=50000 , value=0.0001)
      ),
      conditionalPanel(condition = "input.technical",
                       numericInput(inputId = "mqeps",
                                    label = "MINQUE Convergence Criteria:",
                                    min=0 , max=50000 , value=0.0001)
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





##############################################################################
##
## The server for the shiny app
##
##############################################################################


shinyServer_clme <- function(input, output) {
  
  
  #library("CLME")
  
  
  clme.out <- reactive({
    
    compute1 <- input$compute1
    
    ## Put all the code to run CLME inside this
    if( compute1 > 0 ){
      isolate({
        
        file1   <- input$file1[4]
        solver  <- input$solver
        data1   <- as.matrix( read.csv( file=paste(file1) ) )
        
        yy      <- input$yy
        p1      <- input$p1
        p2      <- input$p2
        q       <- input$q
        nsim    <- input$nsim
        
        if( p2 == '' | p2==0 ){ p2 <- "None" }
        if(  q == '' |  q==0 ){  q <- "None" }
        
        ##
        mySolver <- switch( EXPR=solver ,
                            "Least Squares (LS)"        = "LS",
                            "Least Absolute Value (L1)" = "L1",
                            "General LS (GLS)"          = "GLS",
                            "Asymmetrix LS"             = "asyLS",
                            "L1 approx"                 = "L1eps",
                            "Huber"                     = "huber",
                            "SILF"                      = "SILF",
                            "Chebyshev"                 = "chebyshev",
                            "L-p Power Norm"            = "Lp",
                            "Quantile"                  = "quantile",
                            "Poisson"                   = "poisson")       
        
                
        ## Get the indexes for X2 and U
        parse_idx <- function( arg1 ){
          arg1 <- toupper(arg1)
          if( arg1=="NONE" ){
            p4 <- 0
          } else{
            p2s <- strsplit(arg1, ",")[[1]]
            if( length(p2s) > 1 ){
              p3 <- vector()
              for( ii in 1:length(p2s) ){
                next1 <- sapply( p2s[ii],  FUN=function(x){ strsplit(x, "-") }  )[[1]]
                
                if( length(next1)==1 ){
                  if(  next1 %in% LETTERS ){
                    p3 <- append( p3 , which(next1 == LETTERS) )
                  } else{
                    p3 <- append( p3 , as.numeric(next1) )
                    
                  }
                } else{
                  next1a <- next1[1]
                  next1b <- next1[2]
                  if( next1a%in%LETTERS & next1b%in%LETTERS ){
                    n1a <- which(next1a == LETTERS)
                    n1b <- which(next1b == LETTERS)
                    p3 <- append( p3 , n1a:n1b )
                  } else{
                    p3 <- append( p3 , as.numeric(next1a):as.numeric(next1b) )
                  }
                }
              }
            } else{
              next1 <- sapply( p2s,  FUN=function(x){ strsplit(x, "-") }  )[[1]]
              
              if( length(next1)==1 ){
                if(  next1 %in% LETTERS ){
                  p3 <-  which(next1 == LETTERS)
                } else{
                  p3 <-as.numeric(next1)
                }
              } else{
                next1a <- next1[1]
                next1b <- next1[2]
                if( next1a%in%LETTERS & next1b%in%LETTERS ){
                  n1a <- which(next1a == LETTERS)
                  n1b <- which(next1b == LETTERS)
                  p3  <- n1a:n1b
                } else{
                  p3 <- as.numeric(next1a):as.numeric(next1b)
                }
              }
            }
              p4 <- as.numeric(p3)
          }
          p4
        }
        
        idx_yy <- parse_idx(yy)
        yn <- colnames(data1)[idx_yy]
        
        idx_x1 <- parse_idx(p1)
        x1n <- colnames(data1)[idx_x1]
        
        cov <- ran <- FALSE
        if( p2 != "None" ){
          idx_x2 <- parse_idx(p2)
          x2n <- colnames(data1)[idx_x2]
          cov <- TRUE
        }
        if( q != "None" ){
          idx_u  <- parse_idx(q)
          uun <- colnames(data1)[idx_u]
          ran <- TRUE
        }
        
        ncon   <- length(idx_x1)
                
        ## Create the formula
        
        ## Reorder the levels of X1 if need be
        if( ncon==1 & input$xlevel1 ){
          xlev    <- parse_idx(input$xlevels)
          nlev    <- length(levels(as.factor( data1[,xlev] )))
          xlevels <- data1[1:nlev,xlev]          
          data1[,idx_x1] <- factor( data1[,idx_x1], levels=xlevels )
        }
        
        ## Build the formula
        x1f <- paste( x1n , collapse=" + ")
        
        if( !cov & !ran ){
          ## No extra terms
          formula <- formula( paste( yn, "~", x1f ) )
        } else if( cov & !ran ){
          ## Yes covariates, no random effects
          x2f     <- paste( x2n , collapse=" + ")
          formula <- formula( paste( yn, "~", x1f, "+", x2f ) )
        } else if( !cov & ran ){
          ## No covariates, yes random effects
          uf      <-  paste( "(1|", uun , ")", collapse=" + " )
          formula <- formula( paste( yn, "~", x1f, "+", uf ) )
        } else if( cov & ran ){
          ## Covariates and random effects
          x2f     <- paste( x2n , collapse=" + ")
          uf      <-  paste( "(1|", uun , ")", collapse=" + " )
          formula <- formula( paste( yn, "~", x1f, "+", x2f, "+", uf ) )
        }
              
        ## Input control arguments
        tsf <- lrt.stat
        if( input$tsfunc=="Williams" ){
          tsf <- w.stat
        }
        
        constraints <- list()
        constraints$order      <- tolower(input$order)
        constraints$decreasing <- input$decreasing
        
        rep.idx <- grep( "tree" , constraints$order )
        constraints$order[rep.idx] <- "simple.tree"
        
        gfix <- NULL
        if( input$varssq ){
          idx_grp <- parse_idx(input$gfix)
          gfix <- data[,idx_grp]
        }
        
        emiter <- 500
        mqiter <- 500
        emeps  <- 0.00001
        mqeps  <- 0.00001
        seedvl <- NULL
                
        if( input$technical==TRUE ){
          emiter <- input$emiter
          mqiter <- input$mqiter
          emeps  <- input$emeps
          mqeps  <- input$mqeps
          seedvl <- input$ranseed
        }
              
        ## Run the model
        data2 <- as.data.frame( data1 )

        clme.results <- clme(
          formula=formula, data=data2, gfix=gfix, constraints=constraints, nsim=nsim, 
          verbose = c(FALSE, FALSE, FALSE), tsf=tsf, tsf.ind = w.stat.ind,
          mySolver=mySolver, seed=seedvl, ncon=ncon,
          em.eps=emeps, em.iter=emiter, mq.eps=mqeps, mq.iter=mqiter
        )
        
        clme.results
        
      })
      
      
    }   
  })
  
  
  output$fig1 <- renderPlot({
    
    compute2 <- input$compute2
    
    if( compute2 > 0 ){
      isolate({
        
        ci    <- FALSE
        alpha <- 0.05
        
        if( input$outcheck==TRUE ){
          ciwd  <- input$plotci
          alpha <- input$alpha
        }
        
        plot( clme.out() , ci=ci , alpha=alpha)

        
      })
    }
    
  })
  
  output$summary <- renderPrint({
    
    compute3 <- input$compute3
    
    if( compute3 > 0 ){
      isolate({
        
        alpha <- 0.05        
        
        if( input$outcheck==TRUE ){
          alpha <- input$alpha
        }
        
        summary( clme.out() , alpha=alpha )
      })
    }
    
  })
  
  
  ## Put all the code to run CLME inside this
  # output$data <- renderTable({ clme.out() })
  
}








