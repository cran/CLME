resid.boot <-
function( Y , X1 , X2=NULL , U=NULL , Nks=dim(X1)[1] , Qs=dim(U)[2] , 
          constraints , nsim=1000 , mq.phi=NULL , seed=NULL ){
  
  if( is.numeric(seed) ){
    set.seed(seed)
  }
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix( cbind( X1,X2 ))
  
  K  <- length(Nks)
  
  if( is.null(constraints$A) ){
    constraints <- create.constraints(X1=X1 , X2=X2 , constraints=constraints)
  }
  
  
  # Initial values
  theta <- ginv( t(X)%*%X )%*%t(X)%*%Y
  ssq   <- vector()
  for( k in 1:K ){
    Yk <- Y[ N1[k]:N2[k] ]
    Xk <- X[ N1[k]:N2[k],]  
    ssq[k] <- sum( (Yk - Xk%*%theta)^2 ) / (Nks[k])
  }
  
  ## Obtain the estimates of epsilon and delta
  ssqvec  <- rep(ssq,Nks)    
  
  XSiX <- t(X) %*% (X/ssqvec)
  XSiY <- t(X) %*% (Y/ssqvec)
  
  if( Q > 0 ){
    if( is.null(mq.phi) ){
      mq.phi <- minque( Y=Y , X1=X1 , X2=X2 , U=U , Nks=Nks , Qs=Qs )    
      tsq    <- mq.phi[1:Q]
    } else{
      tsq <- mq.phi[1:Q]
    }
    tsqvec  <- rep(tsq,Qs)    
    C       <- U * tsqvec
    
    U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
    tusu       <- t(U1) %*% U1
    diag(tusu) <- diag(tusu) + 1/tsqvec
    tusui      <- solve(tusu)   
    XSiU <- t(X) %*% (U/ssqvec)
    USiY <- t(U) %*% (Y/ssqvec)
    XPiX <- XSiX - XSiU%*%tusui%*%t(XSiU)
    XPiY <- XSiY - XSiU%*%(tusui%*%USiY)
  } else{
    XPiX <- XSiX
    XPiY <- XSiY
  }
  
  
  # H <- X%*%ginv( XPiX )%*%(t(X)%*%PsiI)
  #eps  <- c( Y - H%*%Y )  
  Yhat <- X %*% ginv( XPiX )%*%XPiY
  eps  <- c(Y - Yhat )
  nu   <- eps
  
  for( i in 1:K ){
    idx      <- N1[i]:N2[i]
    nu[idx]  <- eps[idx] / sd( eps[idx] )
  }
  
  if( Q > 0 ){
    # xi    <- c( t(C) %*% PsiI %*% eps )  
    USiR <- t(U) %*% (eps/ssqvec)
    USiU <- t(U) %*% (U/ssqvec)
    xi   <- tsqvec * ( USiR - USiU%*%(tusui%*%USiR) )
    delta <- xi
    for( i in 1:Q ){
      idx         <- Q1[i]:Q2[i]
      delta[idx]  <- xi[idx] / sd( xi[idx] )
    }
  }
  
  ## Obtain the estimate under the null
  qp.meq <- dim( constraints$A )[1]
  Dmat   <- XPiX
  Deigen <- eigen(Dmat)
  if( sum(Deigen$values <= sqrt(.Machine$double.eps) ) > 0 ){
    idx.neg.lamda <- which(Deigen$values <= sqrt(.Machine$double.eps) )
    Deigen$values[idx.neg.lamda] <- sqrt(.Machine$double.eps)
    Dmat <- Deigen$vectors %*% diag(Deigen$values) %*% t(Deigen$vectors)
  }
  dvec  <- XPiY
  theta.null <- solve.QP(Dmat, dvec, t(constraints$A), meq=qp.meq )$solution
  
  ## Obtain the bootstrap samples
  Y.boot  <- matrix( NA , nrow=N  , ncol=nsim )
  XT.boot <- X%*%theta.null
  
  if( Q > 0 ){
    Qc <- dim(U)[2]
    for( m in 1:nsim ){
      xi.boot   <- sample( delta , replace=TRUE )
      eps.boot  <- sample( nu    , replace=TRUE )
      
      for( i in 1:Q ){
        idx <- Q1[i]:Q2[i]
        xi.boot[idx]   <- sqrt(tsq[i]) * xi.boot[idx]
      }
      for( i in 1:K ){
        idx <- N1[i]:N2[i]
        eps.boot[idx]  <- sqrt(ssq[i]) * eps.boot[idx]
      }
      Y.boot[,m] <- XT.boot + U%*%xi.boot + eps.boot
    }
    
  } else{
    for( m in 1:nsim ){
      eps.boot  <- sample( nu    , replace=TRUE )
      for( i in 1:K ){
        idx <- N1[i]:N2[i]
        eps.boot[idx]  <- sqrt(ssq[i]) * eps.boot[idx]
      }
      Y.boot[,m] <- XT.boot + eps.boot
    }
  }
  
  ## Return the bootstrap samples
  Y.boot
  
}
