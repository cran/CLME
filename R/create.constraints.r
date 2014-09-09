create.constraints <- function( X1, X2 = NULL, constraints ){
  
  P1 <- dim(X1)[2]
  X  <- cbind(X1, X2)
  
  order      <- tolower(constraints$order)
  node       <- constraints$node
  decreasing <- constraints$decreasing
    
  if( sum(order==c("simple", "simple.tree", "umbrella")) == 0 ){
    stop("'order' must be one or more of: simple, simple.tree, umbrella")
  }
  
  if( decreasing==FALSE ){ mult <-  1 }
  if( decreasing==TRUE  ){ mult <- -1 }
  
  A <- matrix( 0, nrow=(P1-1), ncol=dim(X)[2] )
  
  ## Simple order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_K
  if( order=="simple" ){
    for( a1 in 1:(P1-1) ){
      A[a1,a1]   <- -1*mult
      A[a1,a1+1] <-  1*mult    
    }
    B <- matrix( 0, nrow=1 , ncol=dim(X)[2] )
    B[1,1]  <- -1*mult
    B[1,P1] <-  1*mult
    node    <-  1
  }
  
  ## Simple tree order
  ## e.g. mu_1 <= mu_i ; i=2,...,P1
  if( order=="simple.tree" ){
    anti.nodes <- 1:P1
    anti.nodes <- anti.nodes[-node]
    for( a1 in 1:length(anti.nodes) ){
      A[a1,node]           <- -1*mult
      A[a1,anti.nodes[a1]] <-  1*mult    
    }
    B <- A    
  }
  
  ## Umbrella order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_b >= mu_{b+1} >= ... >= mu_K
  if( order=="umbrella" ){    
    for( a1 in 1:(P1-1) ){
      if( a1 < node ){
        A[a1,a1]     <- -1*mult
        A[a1,a1+1]   <-  1*mult      
      }
      if( a1 >= node ){
        A[a1,a1]     <-  1*mult
        A[a1,a1+1]   <- -1*mult    
      }    
    }
    
    B <- matrix( 0, nrow=2 , ncol=dim(X)[2] )
    B[1,1]    <- -1*mult
    B[1,node] <-  1*mult    
    
    B[2,node] <-  1*mult        
    B[2,P1]   <- -1*mult
  }
  
  # Return the constraints object
  constraints <- list( A=A , B=B , order=order , decreasing=decreasing , node=node )
  constraints
}
