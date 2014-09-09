

##
## Define the S3 class 'clme' and some methods
##



"summary.clme" <- function( object , alpha=0 , dec.theta=2 , dec.w=3 , 
                            dec.var=4 , dec.p=4 , ...){
  
  clme <- object
  
  if( is.clme(clme)==FALSE ){
    stop("Argument 'object' is not of class clme.")
  }
  
  # Format the variable names
  theta.names <- names(clme$theta)
  th.fix      <- theta.names==""
  if( is.null(theta.names) ){ th.fix <- rep( 1, length(clme$theta) ) }

  idx <- which( th.fix==1 )
  theta.names[idx] <- paste( "theta_" , idx , sep="" )    
  theta.names2 <- theta.names
  
  # Make sure names are equal length
  max.len <- max( nchar(theta.names) )
  for( i in 1:length( theta.names) ){
    cur.len <- nchar(theta.names[i])
    if( cur.len < max.len ){
      for( j in 1:(max.len-cur.len) ){
        theta.names[i] <- paste( theta.names[i] , " " , sep="" )
      }
    }
  }
  
  # Print the global test results
  aa <- round(clme$ts.glb   , dec.w)
  bb <- round(clme$p.value , dec.p)
  
  t.sprint <- paste( "%." , dec.theta , "f" , sep="" )
  ts.sprint <- paste( "%." , dec.w     , "f" , sep="" )
  v.sprint <- paste( "%." , dec.var   , "f" , sep="" )
  p.sprint <- paste( "%." , dec.p     , "f" , sep="" )
  
  if( length(clme$ts.glb)>1 ){
    cat( "Global Tests: ")
    for( i in 1:length(clme$ts.glb) ){
      cat( "\n W" , i , " = " ,  noquote( sprintf(ts.sprint, aa[i] ) ) ,
           "    p = " , noquote( sprintf(p.sprint, bb[i] ) ) , sep="" )
    }
  } else{
    cat( "Global Test: \n W = " ,  noquote( sprintf(ts.sprint, aa ) ) ,
         "    p = " , noquote( sprintf(p.sprint, bb ) )  )
  }
  
  if( is.null(clme$est.order)==FALSE ){
    cat( "\n" , noquote( clme$est.order  ) )
  }
  
  cat( "\n\nIndividual Tests: \n" )
  
  
  # Print the individual test results
  for( i in 1:length(clme$ts.ind) ){
    cat( " Contrast " , i , ": " , sep="" )
    
    aa <- round( clme$ts.ind[i] , dec.w)
    bb <- round( clme$p.value.ind[i] , dec.p)
    
    # Get the actual constraint being tested
    positive <- which( clme$constraints$A[i,] > 0 )
    negative <- which( clme$constraints$A[i,] < 0 )
    
    for( j in 1:length(positive) ){
      if( abs(clme$constraints$A[i,positive[j]]) != 1 ){
        mult <- abs(clme$constraints$A[i,positive[j]])
        cat( mult , "*" , theta.names2[positive[j]] , sep="")
      } else{
        cat( theta.names2[positive[j]] , sep="")
      }    
      
      if( j < length(positive)){ cat(" + " , sep="")  }
    }
    
    cat( " - ")
    
    for( j in 1:length(negative) ){
      if( abs(clme$constraints$A[i,negative[j]]) != 1 ){
        mult <- abs(clme$constraints$A[i,negative[j]])
        cat( mult , "*" , theta.names2[negative[j]] , sep="")
      } else{
        cat( theta.names2[negative[j]] , sep="")
      }
      
      if( j < length(negative)){ cat(" - " , sep="")  }
    }
    
    cat( "\n    W = "  , noquote( sprintf(ts.sprint, aa ) ) ,
         "    p = "    , noquote( sprintf(p.sprint, bb ) ) , "\n" , sep="") 
  }
  
  # Print the thetas:
  cat( "\nTheta Coefficients: \n") 
  for( i in 1:length(clme$theta) ){
        
    cat( " " , theta.names[i] ,  sep="" )    
    cat(" = ")
    
    if( clme$theta[i] > 0 ){ cat(" ") }
    rnd.th  <- round(clme$theta[i], dec.theta )    
    rnd.lcl <- round(clme$theta[i] - qnorm(1-alpha/2)*sqrt(clme$cov.theta[i,i])
                     , dec.theta )
    rnd.ucl <- round(clme$theta[i] + qnorm(1-alpha/2)*sqrt(clme$cov.theta[i,i])
                     , dec.theta )
    
    if( alpha >0 & alpha < 1 ){
      cat( sprintf(t.sprint , rnd.th)  ,  "   ( " , sprintf(t.sprint , rnd.lcl)
           , " , " , sprintf(t.sprint , rnd.ucl) , " )\n" , sep="" )
    } else{
      cat( sprintf(t.sprint , rnd.th)  ,  "\n" , sep="" )
    }
    
    
  }
  
  cat( "\n\n")
  
  cat( "Variances (ssq = sigma^2 , tsq = tau^2): \n") 
  # Print the sigmas and taus
  for( i in 1:length(clme$ssq ) ){
    rnd.var <- round(clme$ssq[i] , dec.var)
    cat( " ssq_" , i , " = " , sprintf(v.sprint , rnd.var) , "\n" , sep="" )
  }
  if( is.null(clme$tsq)==FALSE ){
    for( i in 1:length(clme$tsq ) ){
      rnd.var <- round(clme$tsq[i] , dec.var)
      cat( " tsq_" , i , " = " , sprintf(v.sprint , rnd.var) , "\n" , sep="" )
    }
  }
}

