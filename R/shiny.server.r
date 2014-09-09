



shiny.server <- function(input, output) {
  
  
  library("CLME")

  
  clme.out <- reactive({
    
    compute1 <- input$compute1
    
    ## Put all the code to run CLME inside this
    if( compute1 > 0 ){
      isolate({
        
        
        file1  <- input$file1[4]
        method <- input$method
        data1  <- as.matrix( read.csv( file=paste(file1) ) )
        
        p1 <- input$p1
        p2 <- input$p2
        q  <- input$q
        Mb <- input$mboot
        
        nobs <- nrow(data1)
        
        
        ## Create the data matrices
        Y <- as.matrix( data1[,1] )
        
        x1s <- 2
        x1t <- x1s + p1 - 1
        X1  <- as.matrix( data1[,x1s:x1t] )
        colnames(X1) <- colnames(data1[x1s:x1t])
        
        if( p2 > 0 ){
          x2s <- x1t + 1
          x2t <- x2s + p2 - 1
          X2  <- as.matrix( data1[,x2s:x2t] )
          colnames(X2) <- colnames(data1[x2s:x2t])
        } else{
          X2  <- NULL
          x2t <- x1t
        }
        
        if( q > 0 ){
          us <- x2t + 1
          ut <- us  + q - 1
          U  <- as.matrix( data1[,us:ut] )
        } else{
          U  <- NULL
          Qs <- NULL
        }
        
        
        ## Input control arguments
        if( input$tsfunc=="Williams" ){
          tsf <- w.stat
        } else{
          tsf <- lrt.stat
        }
        
        constraints <- list()
        constraints$order      <- tolower(input$order)
        constraints$decreasing <- input$decreasing
        
        rep.idx <- grep( "tree" , constraints$order )
        constraints$order[rep.idx] <- "simple.tree"
        
        
        if( input$varSSQ==TRUE ){
          Nks <- eval(parse(text= paste("c(",input$Nks,")" , sep="")) )
        } else{
          Nks <- nrow(X1)
        }
        if( q>0 & input$varTSQ==TRUE ){
          Qs <- eval(parse(text= paste("c(",input$Qs,")" , sep="")) )
        } else{
          Qs <- ncol(U)
        }
        
        if( input$technical==TRUE ){
          em.iter <- input$em.iter
          mq.iter <- input$mq.iter
          em.eps  <- input$em.eps
          mq.eps  <- input$mq.eps
          seedvl  <- input$ranseed
        } else{
          em.iter <- 500
          mq.iter <- 500
          em.eps  <- 0.00001
          mq.eps  <- 0.00001
          seedvl  <- NULL
        }
        
        
        ## Run the model
        clme.results <- constrained.lme(
          method=method, Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs,
          constraints=constraints, nsim=Mb, em.eps=em.eps, em.iter=em.iter,
          mq.eps=mq.eps, mq.iter=mq.iter,  verbose = c(FALSE, FALSE, FALSE), 
          tsf=tsf, tsf.ind = w.stat.ind , seed=seedvl ) 
        
        
        clme.results
        
      })
      
      
    }   
  })
  
  
  output$fig1 <- renderPlot({ 
    
    compute2 <- input$compute2
    
    if( compute2 > 0 ){
      isolate({
        
        ciwd  <- 0
        alpha <- 0.05
        
        if( input$outcheck==TRUE ){
          if( input$plotci==TRUE ){
            ciwd  <- 0.05  
          }
          alpha <- input$alpha
        }
        
        plot( clme.out() , ci.wd=ciwd , alpha=alpha)
      })
    }
    
  })
  
  output$summary <- renderPrint({ 
    
    compute3 <- input$compute3
    
    if( compute3 > 0 ){
      isolate({
        
        alpha <- 0.05
        decp <- 4
        dect <- 2
        decw <- 3
        decv <- 4
        
        
        if( input$outcheck==TRUE ){
          alpha <- input$alpha
          decp  <- input$decp
          dect  <- input$dect
          decw  <- input$decw
          decv  <- input$decv
        }     
        
        summary( clme.out() , dec.theta=dect , dec.w=decw , dec.var=decv , dec.p=decp , alpha=alpha )
      })
    }
    
  })
    
  
  ## Put all the code to run CLME inside this
  # output$data <- renderTable({ clme.out() })
  
}

