


pava.simple.order <- function( theta , cov.theta=diag(length(theta)) , 
                               node=NULL , decreasing=FALSE, hp=FALSE ){
  if(hp==TRUE){
    ## HWANG & PEDDADA (1994) METHOD
    p1        <- length(theta)
    iso.theta <- rep( 0 , p1)
    
    if( decreasing==FALSE ){ fun1 <- min ; fun2 <- max }
    if( decreasing==TRUE  ){ fun1 <- max ; fun2 <- min }
    
    for( ii in 1:p1 ){
      tt.vec <- rep( 0 , length(ii:p1) )
      idx1   <- 0
      for( tt in ii:p1 ){
        idx1   <- idx1 + 1
        ss.vec <-  rep( 0 , length(1:tt) )
        idx2   <- 0
        for( ss in 1:tt ){
          idx2   <- idx2 + 1
          sub.th <- theta[ss:tt]
          sub.cv <- solve( cov.theta[ss:tt , ss:tt , drop=F] )
          one    <- rep( 1, length(ss:tt) )
          ss.vec[idx2] <- ( (t(one) %*% sub.cv %*% sub.th) / 
                              (t(one) %*% sub.cv %*% one   ) )
        }
        tt.vec[idx1] <- fun2(ss.vec)
      }
      iso.theta[ii] <- fun1(tt.vec)
    }
    
  } else{
    ## ONLY USE DIAGONAL OF COVARIANCE
    wt <- 1/diag( cov.theta )
    if( is.null(wt) ){ wt <- rep( 1 , length(theta)) }
      iso.theta <- pava( theta , wt , decreasing=decreasing, long.out=FALSE, stepfun=FALSE)
  }
  
  ## Return theta
  return( iso.theta )
  
}




pava.simple.tree <- function( theta , cov.theta=diag(length(theta)) , 
                              node=1 , decreasing=FALSE, hp=FALSE ){
  
  if( hp==TRUE ){
    ## HWANG & PEDDADA (1994) METHOD
    if( is.null(cov.theta) ){ cov.theta <- diag(length(theta)) }
    if( is.null(decreasing)){
      stop("''decreasing'' in pava.simple.tree must be specified")
    }
    
    p1        <- length(theta)
    iso.theta <- theta
    idxs      <- 1:p1
    
    if( decreasing==FALSE ){ fun1 <- min ; fun2 <- max }
    if( decreasing==TRUE  ){ fun1 <- max ; fun2 <- min }
    
    # Reorder points
    theta.ord  <- c( theta[node] , theta[-node] )
    
    cv.ord     <- cov.theta
    cv.ord[,1] <- cov.theta[,node]
    cv.ord[1,] <- cov.theta[node,]
    cv.ord[2:p1,2:p1] <- cov.theta[idxs[-node],idxs[-node]]
    
    # Estimate the node
    node.vec <- vector()
    idx <- 0
    for( kk in 1:p1 ){
      
      combs    <- t( combn( 1:p1 , kk ) )
      row.comb <- apply( combs , 1 , FUN=function(x){ sum(x==node) > 0 } )
      combs    <- combs[ row.comb==T , , drop=FALSE]
      nk       <- nrow(combs)
      
      for( cc in 1:nk ){
        idx  <- idx+1
        inds <- c( combs[cc,])
        
        sub.th <- theta[inds]
        sub.cv <- solve( cov.theta[inds , inds , drop=FALSE] )
        one    <- rep( 1, length(inds) )
        node.vec[idx] <- ( (t(one) %*% sub.cv %*% sub.th) / 
                             (t(one) %*% sub.cv %*% one   ) )
      }
    }
    
    # Update the rest of the thetas
    iso.theta[node] <- fun1(node.vec)  
    for( ii in idxs[-node] ){
      iso.theta[ii] <- fun2(  c(iso.theta[node] , iso.theta[ii])  )    
    }
  } else{
    ## ONLY USE DIAGONAL OF COVARIANCE
    wt <- 1/diag( cov.theta )
    
    if( is.null(wt) ){ wt <- rep( 1 , length(theta)) }
    if( is.null(decreasing)){
      stop("''decreasing'' in pava.simple.tree must be specified")
    }
    
    p1        <- length(theta)
    iso.theta <- theta
    idxs      <- 1:p1
    
    if( decreasing==FALSE ){ fun1 <- min ; fun2 <- max }
    if( decreasing==TRUE  ){ fun1 <- max ; fun2 <- min }
    
    # Reorder points
    theta.ord  <- c( theta[node] , theta[-node] )
    
    wt.ord       <- wt
    wt.ord[1]    <- wt[ node ]
    wt.ord[2:p1] <- wt[ idxs[-node] ]
    
    # Estimate the node
    if( theta[node] == fun1(theta) ){
      iso.theta[node] <- theta[node]
    } else{
      
      node.vec <- vector()
      idx      <- 0
      
      for( kk in 1:p1 ){
        
        combs    <- t( combn( 1:p1 , kk ) )
        row.comb <- apply( combs , 1 , FUN=function(x){ sum(x==node) > 0 } )
        combs    <- combs[ row.comb==T , , drop=FALSE]
        nk       <- nrow(combs)
        
        for( cc in 1:nk ){
          idx  <- idx+1
          inds <- c( combs[cc,])
          
          sub.th <- theta[inds]
          sub.wt <- wt[ inds ]
          node.vec[idx] <- sum( sub.th * sub.wt ) / sum( sub.wt )
        }
      }
      
      iso.theta[node] <- fun1(node.vec)
      
    }
    
    
    # Update the rest of the thetas
    for( ii in idxs[-node] ){
      iso.theta[ii] <- fun2(  c(iso.theta[node] , iso.theta[ii])  )    
    }
  }
  
  ## Return theta
  return(iso.theta)
}



pava.umbrella <- function( theta , cov.theta=diag(length(theta)) ,
                           node=1 , decreasing=FALSE, hp=FALSE ){
  if( hp==TRUE ){
    ## HWANG & PEDDADA (1994) METHOD
    if( is.null(cov.theta) ){ cov.theta <- diag(length(theta)) }
    if( is.null(decreasing)){
      stop("''decreasing'' in pava.umbrella must be specified")
    }
    p1        <- length(theta)
    iso.theta <- theta
    
    if( decreasing==FALSE ){ dec2 <- TRUE  ; fun1 <- max }
    if( decreasing==TRUE  ){ dec2 <- FALSE ; fun1 <- min }
    
    theta.temp1 <- pava.simple.order( iso.theta[1:node]  , 
                                      cov.theta[1:node  , 1:node  , drop=FALSE] ,
                                      decreasing=decreasing , hp=hp)
    
    theta.temp2 <- pava.simple.order( iso.theta[node:p1] , 
                                      cov.theta[node:p1 , node:p1 , drop=FALSE] , 
                                      decreasing=dec2 , hp=hp  )
    
    iso.theta[node] <- fun1( theta.temp1[node] , theta.temp2[1] )
    
    
    iso.theta[1:(node-1)]  <- theta.temp1[1:(node-1)]
    iso.theta[(node+1):p1] <- theta.temp2[2:(p1-node+1)]
    
  } else{
    ## ONLY USE DIAGONAL OF COVARIANCE
    wt <- 1/diag( cov.theta )
    if( is.null(wt) ){ wt <- rep( 1 , length(theta)) }
    if( is.null(decreasing)){
      stop("''decreasing'' in pava.umbrella must be specified")
    }
    p1        <- length(theta)
    iso.theta <- theta
    
    if( decreasing==FALSE ){ dec2 <- TRUE  ; fun1 <- max }
    if( decreasing==TRUE  ){ dec2 <- FALSE ; fun1 <- min }
    
    theta.temp1 <- pava( iso.theta[1:node]  , wt[1:node]  , decreasing=decreasing, long.out=FALSE, stepfun=FALSE)
    theta.temp2 <- pava( iso.theta[node:p1] , wt[node:p1] , decreasing=dec2      , long.out=FALSE, stepfun=FALSE)
    
    
    iso.theta[node] <- fun1( theta.temp1[node] , theta.temp2[1] )
    
    iso.theta[1:(node-1)]  <- theta.temp1[1:(node-1)]
    iso.theta[(node+1):p1] <- theta.temp2[2:(p1-node+1)]
  }

  ## Return theta
  return(iso.theta)
}


