clme.em.mixed <- function( method, Y, X1, X2 = NULL, U = NULL, Nks = dim(X1)[1],
                     Qs = dim(U)[2], constraints, mq.phi = NULL, tsf = lrt.stat, 
                     tsf.ind = w.stat.ind, pav.alg, hp=FALSE, em.iter = 500, 
                     em.eps = sqrt(.Machine$double.eps), verbose = FALSE ){
  
  if( verbose==TRUE ){
    message("Starting EM algorithm")    
  }
  
  N  <- sum(Nks)
  N1 <- 1 + cumsum(Nks) - Nks
  N2 <- cumsum(Nks)
  
  Q  <- length(Qs)
  Q1 <- 1 + cumsum(Qs) - Qs
  Q2 <- cumsum(Qs)
  
  X  <- as.matrix(cbind(X1, X2))
  theta.names <- NULL
  if( is.null(colnames(X))==FALSE ){
    theta.names <- colnames(X)
  }
  
  K  <- length(Nks)
  
  # Unpack the constraints
  A          <- constraints$A
  B          <- constraints$B
  decreasing <- constraints$decreasing
  node       <- constraints$node
  
  P1 <- dim(X1)[2]
  P2 <- dim(X2)[2]
  
  # Initialize values
  XX    <- t(X) %*% X
  XXi   <- solve( XX )
  theta <- ginv( t(X)%*%X) %*% ( t(X)%*%Y )
  
  R   <- Y-X%*%theta
  ssq <- apply( as.matrix(1:K, nrow=1), 1 , 
                FUN=function(k, N1, N2, R ){ sum( R[N1[k]:N2[k]]^2 ) / Nks[k] } ,
                N1, N2, R)
  
  ssqvec  <- rep(  ssq,Nks)  
    
    UU <- U %*% t(U)
    
    if( is.null(mq.phi) ){
      mq.phi <- minque( Y=Y, X1=X1, X2=X2, U=U, Nks=Nks, Qs=Qs)
    }
    tsq    <- mq.phi[1:Q]
    tsqvec <- rep(tsq,Qs)
    
    XSX <- t(X)%*%(X*ssqvec)    
    XU  <- t( apply( (t(X)%*%(U) ) , 1 , FUN=function(x,tv){ x*sqrt(tv) } , tsqvec) )
    XPX <- (XU%*%t(XU)) + XSX    
  
  cov.theta <- XXi %*% (XPX) %*% XXi
  
  
  theta1 <- theta
  ssq1   <- ssq
  tsq1   <- tsq
  
  
  
  
  # Being the EM Algorithm convergence loop
  CONVERGE  <- 0
  iteration <- 0
  while( CONVERGE==0 ){
    iteration <- iteration+1
    
    R <- Y-X%*%theta
    
    if( verbose==TRUE ){
      message( "EM iteration " , iteration)
    }
      
    # Step 1: Estimate Sigma
      SiR    <- R / ssqvec
      # U' * SigI * U
      U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
      tusu       <- t(U1) %*% U1
      diag(tusu) <- diag(tusu) + 1/tsqvec
      tusui      <- solve(tusu)
      PiR        <- SiR - (U %*% (tusui %*% (t(U)%*%SiR))) * (1/ssqvec)
      # PsiI
      BU    <- tusui %*% t(U)
      UBU   <- apply( as.matrix(1:sum(Nks)) , 1 , FUN=function(kk,uu,bu){ sum( uu[kk,]*bu[,kk] ) } , U, BU )
      SUBUS <- 1/ssqvec - UBU / ssqvec^2
      trace.vec <- PiR^2 - SUBUS      
    
    
    # trace.vec <- diag( PsiI%*%( R%*%t(R) )%*%PsiI - PsiI )
    ssq <- apply( as.matrix(1:K, nrow=1), 1 , 
                  FUN=function(k, ssq, Nks, N1, N2 , trv){
                    idx <- N1[k]:N2[k]
                    ssq[k] + ( (ssq[k]^2)/(Nks[k]) )*sum( trv[idx] )
                    } ,
                  ssq , Nks , N1 , N2 , trace.vec )
    
    ssqvec  <- rep(  ssq,Nks)
        
    # Step 2a: Estimate Thetas
    # Update the blocks
    SiR   <- R / ssqvec
    X1SiR <- t(X1) %*% SiR

      X1SiU <- t(X1) %*% (U/ssqvec)  #  X1SU
      USiR  <- t(U)  %*% SiR
      
      U1         <- apply( U , 2 , FUN=function(x,sq){x*sq} , 1/sqrt(ssqvec) )
      tusu       <- t(U1) %*% U1
      diag(tusu) <- diag(tusu) + 1/tsqvec
      tusui      <- solve(tusu)    
      
      
      #theta[1:P1]  <- theta1[ 1:P1] + ginv(t(X1)%*%SigmaI%*%X1) %*% ((t(X1)%*%PsiI)%*%R )
      theta[1:P1]  <- theta1[1:P1] + ginv( t(X1)%*%(X1/ssqvec) ) %*% (X1SiR - X1SiU%*%(tusui%*%USiR))
      
      
      if( is.null(X2)==FALSE ){
        #theta[(P1+1):(P1+P2)] <- ( theta1[ (P1+1):(P1+P2)] + 
        #                              ginv(t(X2)%*%SigmaI%*%X2) %*% ((t(X2)%*%PsiI)%*%R ) )
        X2SiU <- t(X2) %*% (U/ssqvec)
        X2SiR <- t(X2) %*% SiR
        theta[(P1+1):(P1+P2)] <- ( theta1[(P1+1):(P1+P2)] + 
                                     ginv(t(X2)%*%(X2/ssqvec)) %*% (X2SiR - X2SiU%*%(tusui%*%USiR)) )
      }
    
    
    
    
    # Step 2b: Estimate Tau
      # USiR  <- t(U)  %*% SiR # <-- previously calcualted
      USiU <- t(U)  %*% (U/ssqvec)
      UPiR <- USiR - USiU%*%(tusui%*%USiR)
      trace.vec.tau <- (UPiR^2) - diag(USiU) + diag( USiU%*%tusui%*%USiU )
      
      for( q in 1:Q ){
        # tau.idx <- Q1[q]:Q2[q]
        #Uq      <- as.matrix( U[,tau.idx] )
        #cq      <- dim(Uq)[2]
        #sumdg   <- sum(diag( t(Uq)%*%( PsiI%*%R%*%t(R)%*%PsiI - PsiI )%*%Uq ))      
        #tsq[q]  <- tsq1[q] + ((tsq1[q]^2)/cq)*sumdg
        tau.idx <- Q1[q]:Q2[q]
        cq      <- length(tau.idx)
        tsq[q]  <- tsq1[q] + ((tsq1[q]^2)/cq)*sum(trace.vec.tau[tau.idx])
      }
    
    
    
    XSiX <- t(X) %*% (X/ssqvec)
    XSiY <- t(X) %*% (Y/ssqvec)
    
      tsqvec  <- rep(tsq,Qs)  
      
      XSX <- t(X)%*%(X*ssqvec)    
      XU  <- t( apply( (t(X)%*%(U) ) , 1 , FUN=function(x,tv){ x*sqrt(tv) } , tsqvec) )
      XPX <- (XU%*%t(XU)) + XSX  
      
      XSiU <- t(X) %*% (U/ssqvec)
      USiY <- t(U) %*% (Y/ssqvec)
      XPiX <- XSiX - XSiU%*%tusui%*%t(XSiU)
      XPiY <- XSiY - XSiU%*%(tusui%*%USiY)
      
    
    
    
    ## Apply order constraints / isotonization
      # PAVA constraints   
      if( method=="PAVA" ){
        cov.theta <- XXi %*% (XPX) %*% XXi
        theta[1:P1] <- pav.alg( theta[1:P1] , cov.theta[1:P1,1:P1,drop=FALSE] , node , decreasing , hp) 
      } else{
      # QPE constrained optimization   
      Dmat  <- XPiX
      Deigen <- eigen(Dmat)
      if( sum(Deigen$values <= sqrt(.Machine$double.eps) ) > 0 ){
        idx.neg.lamda <- which(Deigen$values < sqrt(.Machine$double.eps) )
        Deigen$values[idx.neg.lamda] <- sqrt(.Machine$double.eps)
        Dmat <- Deigen$vectors %*% diag(Deigen$values) %*% t(Deigen$vectors)
      }
      # dvec  <- t( t(Y)%*%PsiI%*%X )
      dvec  <- XPiY
      theta <- solve.QP(Dmat, dvec, t(A))$solution    
    }
    
    # Evaluate some convergence criterion
    rel.change <- abs(theta - theta1)/theta1
    if( mean(rel.change) < em.eps || iteration >= em.iter ){
      CONVERGE <- 1
    } else{
      theta1 <- theta
      ssq1   <- ssq
      tsq1   <- tsq
    }
    
  } # End converge loop (while)
  
  if( verbose==TRUE ){
    message("EM Algorithm ran for " , iteration , " iterations." )
  }
  
  theta        <- c(theta)
  names(theta) <- theta.names
  
    XSX <- t(X)%*%(X*ssqvec)    
    XU  <- t( apply( (t(X)%*%(U) ) , 1 , FUN=function(x,tv){ x*sqrt(tv) } , tsqvec) )
    XPX <- (XU%*%t(XU)) + XSX    
  
  cov.theta <- XXi %*% (XPX) %*% XXi
  # cov.theta <- XXi %*% ( t(X) %*% Psi %*% X ) %*% XXi  
  
  # Compute test statistic
  ts.glb <- tsf( theta=theta, cov.theta=cov.theta, B=B, A=A, Y=Y, X1=X1, 
                X2=X2, U=U, tsq=tsq, ssq=ssq, Nks=Nks, Qs=Qs  )
  
  ts.ind <- tsf.ind(theta=theta, cov.theta=cov.theta, B=B, A=A, Y=Y, X1=X1, 
                   X2=X2, U=U, tsq=tsq, ssq=ssq, Nks=Nks, Qs=Qs )
  
  # Return the results
  return.obj <- list(theta=theta, ssq=ssq, tsq=tsq,
                     cov.theta=cov.theta, ts.glb=ts.glb, ts.ind=ts.ind )
  
  return.obj

}
