
clme <- function( x=NULL , theta=numeric(0), ssq=numeric(0), tsq=numeric(0), 
                  cov.theta=matrix(numeric(0)), ts.glb=numeric(0), 
                  ts.ind=numeric(0), p.value=numeric(0), p.value.ind=numeric(0), 
                  constraints=list(), method=character(0), est.order=character(0) ){
  
  if( is.clme(x) ){
    return(x)
  } else{  
        
    result <- list()
    
    result$theta       <- theta
    result$ssq         <- ssq
    result$tsq         <- tsq
    result$cov.theta   <- cov.theta
    result$ts.glb       <- ts.glb
    result$ts.ind       <- ts.ind
    result$p.value     <- p.value
    result$p.value.ind <- p.value.ind
    result$constraints <- constraints
    result$method      <- method
    result$est.order   <- est.order
    
    class(result) <- "clme"
    return(result)
    
  }
}



is.clme <- function(x) inherits(x, "clme")




as.clme <- function( x , ... ){
  
  if( is.clme(x) ){
    return(x)
  } else{
    
    err.flag  <- 0
    flagTheta <- flagSsq <- flagTsq <- flagCov <- flagW1 <- flagW2 <- flagP1 <- flagP2 <- flagConst <- ""
    
    if( !is.numeric(x$theta) ){
      err.flag  <- 1
      flagTheta <- " theta must be numeric\n"
      x$theta   <- numeric(0)
    }
    
    if( !is.numeric(x$ssq) ){
      err.flag <- 1
      flagSsq  <- " ssq must be numeric (scaler or vector)\n"
      x$ssq    <- numeric(0)
    } 
    
    if( !is.null(x$tsq) & !is.numeric(x$tsq) ){
      err.flag <- 1
      flagTsq  <- " if present, tau must be numeric (scaler or vector)\n"
      x$tsq    <- NULL
    }
    
    if( !is.matrix(x$cov.theta) || !is.numeric(x$cov.theta) ||
          nrow(x$cov.theta) != ncol(x$cov.theta) ||
          nrow(x$cov.theta) != length(x$theta)   ||
          sum(sum(abs(x$cov.theta - t(x$cov.theta)))) > sqrt(.Machine$double.eps) ){
      err.flag    <- 1
      flagCov     <- " cov.theta must be square, symmetric, numeric matrix with dimensions equal to length of theta\n"
      x$cov.theta <- matrix( numeric(0) , nrow=length(x$theta) , ncol=length(x$theta) )
    }
    
    if( !is.numeric(x$ts.glb) ){
      err.flag <- 1
      flagW1   <- " ts.glb must be numeric\n"
      x$ts.glb  <- numeric(0)
    } 
    
    if( !is.numeric(x$ts.ind) ){
      err.flag <- 1
      flagW2   <- " ts.ind must be numeric (scaler or vector)\n"
      x$ts.ind  <- numeric(0)
    }
    
    if( !is.numeric(x$p.value) || length(x$p.value) != length(x$ts.glb) ){
      err.flag   <- 1
      flagP1     <- " p.value must be numeric and of same length as ts.glb\n"
      x$p.value  <- numeric(0)
    } 
    
    if( !is.numeric(x$p.value.ind) || length(x$p.value.ind) != length(x$ts.ind) ){
      err.flag       <- 1
      flagP2         <- " p.value.ind must be numeric and of same length as ts.ind\n"
      x$p.value.ind  <- numeric(0)
    } 
    
    if( !is.list(x$constraints) ){
      err.flag        <- 1
      flagConst       <- " constraints must be list\n"
      x$constraints   <- list( A=matrix( numeric(0) ) )
    } else{
      cnames <- names(x$constraints)
      if( sum(cnames=="A") != 1 ){
        err.flag        <- 1
        flagConst       <- " constraints must contain element A\n"
        x$constraints$A <- matrix( numeric(0) , nrow=length(x$ts.ind ))
      }
    }
    
    
    if( err.flag==1 ){
      err.mssg <- paste( "coercing 'x' to class 'clme' produced errors: \n", 
                         flagTheta, flagSsq, flagTsq, flagCov, flagW1,
                         flagW2, flagP1, flagP2, flagConst, "output may not be valid." , sep = "")
      # warning(warn, sys.call(-1))
      warning( err.mssg )
    }
    
    class(x) <- "clme"
    return(x)
    
  }
}

