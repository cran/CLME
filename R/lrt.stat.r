

lrt.stat <- function( theta , A , Y , X1 , X2 , U ,
                      tsq , ssq , Nks , Qs , ... ){
  
  # Unpack some values
  N  <- sum(Nks)
  N1 <- 1+cumsum(Nks)-Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1+cumsum(Qs)-Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix( cbind( X1,X2 ))
  
  K  <- length(Nks)
  
  # Create covariance matrices
  ssqvec  <- rep(ssq,Nks)    
  
  XSiX <- t(X) %*% (X/ssqvec)
  XSiY <- t(X) %*% (Y/ssqvec)
  
  if( Q > 0 ){
    tsqvec  <- rep(tsq,Qs)  
    
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
  
  
  # Estimate of theta under the null
  qp.meq <- dim( A )[1]
  Dmat   <- XPiX
  Deigen <- eigen(Dmat)
  if( sum(Deigen$values <= (.Machine$double.eps)^(1/3) ) > 0 ){
    idx.neg.lamda <- which(Deigen$values <= (.Machine$double.eps)^(1/3) )
    Deigen$values[idx.neg.lamda] <- (.Machine$double.eps)^(1/3)
    Dmat <- Deigen$vectors %*% diag(Deigen$values) %*% t(Deigen$vectors)
  }
  dvec  <- XPiY
  theta.null <- solve.QP(Dmat,dvec, t(A) ,meq=qp.meq )$solution
  theta.diff <- theta - theta.null
  
  # Calculate likelihood function
  ts.stat <- t(theta.diff) %*% XPiX %*% theta.diff
  
  # Return test statistic
  c(ts.stat)
  
}


