constrained.lme <-
function( method="QPE" , Y , X1 , X2=NULL , U=NULL , Nks=dim(X1)[1] , 
          Qs=dim(U)[2] , constraints=list() , nsim=1000 ,
          em.eps=sqrt(.Machine$double.eps) , em.iter=500, 
          mq.eps=sqrt(.Machine$double.eps) , mq.iter=500 , 
          verbose=c(FALSE,FALSE,FALSE) , tsf=lrt.stat , tsf.ind=w.stat.ind ,
          pav.alg=NULL, hp=FALSE , seed=NULL
          ){
  
  if( is.numeric(seed) ){
    set.seed(seed)
  }
  
  # Capitalize the method
  method <- toupper(method)
  
  # If only one element for verbose specified, fill the rest with FALSEs
  if( length(verbose)<3 ){
    verbose <- c(verbose, rep(FALSE, 3-length(verbose) ) )
  }  
  
  # Create the full X matrix
  X <- cbind(X1,X2)

  
  ## Assess the constraints
  cust.const <- is.matrix(constraints$A)
  p1 <- dim(X1)[2]
  
  if( cust.const == FALSE ){
    
    # Constraints are non-null, but A and B are not provided
    # Determine which other elements are missing/needed
    
    if( is.null(constraints$order) ){
      warning( "'constraints$order' is NULL, program will run search for 
               ''simple'' and ''umbrella'' orders")
      constraints$order <- c("simple" , "umbrella" )      
    }
    
    if( is.null(constraints$node) ){
      warning( "'constraints$node' is NULL, program will run search for node")
      constraints$node <- 1:p1
    } else{
      search.node <- FALSE
    }
    
    if( is.null(constraints$decreasing) ){
      warning( "'constraints$decreasing' is NULL, program will run search for
               TRUE and FALSE")      
      constraints$decreasing <- c(TRUE,FALSE)      
    }
    
  }
  
  ## Make sure test stat function is okay
  if( is.function(tsf)==FALSE ){
    stop("'tsf' is not a valid function")
  }
  if( is.function(tsf.ind)==FALSE ){
    stop("'tsf.ind' is not a valid function")
  }
  
  ## Revert to LRT if:  w.stat requested, custom constraints given, and no B provided.
  ## Revert to QPE if:  PAVA requested, custom constraints given, and no pav.alg provided.
  if( cust.const==TRUE ){
    if( identical( tsf , w.stat ) & constraints$B==NULL ){
        warning("Williams type statistic selected with custom constraints, but 
              'constraints$B' is NULL. Reverting to LRT ststistic")
      tsf <- lrt.stat
    }
    if( method=="PAVA" & is.function(pav.alg)==FALSE ){
      warning("PAVA requested with custom constraints, but no pava algorithm 
              provided. Reverting to QPE")
      method <- "QPE"
    }
  }
  
  
  
  
  ## Set up search grid if using defaults
  if( cust.const==FALSE ){
    search.grid <- expand.grid( constraints$order , 
                                constraints$decreasing ,
                                constraints$node )  
    search.grid[,1] <- as.character(search.grid[,1])
    
    # Remove duplicates / extraneous
    # "simple" doesn't need node
    idx         <- 1*(search.grid[,1]=="simple"   &  search.grid[,3] > 1)
    search.grid <- search.grid[ idx==0 , , drop=FALSE]
    
    # "umbrella" with node=1 or node=p1 is covered by simple order
    if( sum(constraints$order=="simple") >= 1 ){
      idx <- 1*((search.grid[,1]=="umbrella" & search.grid[,3] == 1) + 
                (search.grid[,1]=="umbrella" & search.grid[,3] == p1))
      search.grid <- search.grid[ idx==0 , , drop=FALSE]
    } else{
      idx <- 1*(search.grid[,1]=="umbrella" & search.grid[,3] == 1)
      search.grid[idx,1] <- rep( "simple" , sum(idx) )
      idx <- 1*(search.grid[,1]=="umbrella" & search.grid[,3] == p1)
      search.grid <- search.grid[ idx==0 , , drop=FALSE]
    }
    
    # Move simple.tree to the bottom
    idx <- search.grid[,1]=="simple.tree"
    if( sum(idx)>0 ){
      search.grid <- rbind( search.grid[idx==0, , drop=FALSE] ,
                            search.grid[idx==1, , drop=FALSE] )
    }
    
    ## A check for duplicate rows here may be wise
    MNK <- dim( search.grid )[1]  
    
  } else{
    MNK <- 1
    loop.const <- est.const <- constraints
    loop.PAV   <- pav.alg
  }
  
  
  ##
  ## End preparation steps, begin the analysis
  ##
  
  ## Obtain tau if needed
  if( is.null(U)==FALSE ){
    if( is.null(Qs) ){ Qs <- dim(U)[2]  }
    mq.phi <- minque( Y=Y , X1=X1 , X2=X2 , U=U , Nks=Nks , Qs=Qs ,
                      mq.eps=mq.eps , mq.iter=mq.iter , verbose=verbose[2] )
  } else{
    mq.phi <- NULL
  }
  
  ## EM for the observed data
  if( verbose[1]==TRUE ){
    print( paste( "Starting EM Algorithm for observed data." , sep=""))
  }
  
  
  if( is.null(U) ){
    clme.em <-   clme.em.fixed
  } else{
    clme.em <-   clme.em.mixed
  }
  
  ## Loop through the search grid
  est.order <- NULL
  ts.max     <- -Inf

  
  for( mnk in 1:MNK ){
    
    if( cust.const==FALSE ){
      
      grid.row <- list( order     = search.grid[mnk,1], 
                        node      = search.grid[mnk,3], 
                        decreasing= search.grid[mnk,2])       

      loop.const <- create.constraints( X1=X1 , X2=X2 , constraints=grid.row  )
      
      if( method=="PAVA" & is.function(pav.alg)==FALSE ){
        loop.PAV <- switch(
            EXPR=search.grid[mnk,1] ,
              "simple"      = pava.simple.order ,
              "simple.tree" = pava.simple.tree ,
              "umbrella"    = pava.umbrella
        )
      }
      
    }
        
    clme.temp <- clme.em( method=method , Y=Y , X1=X1 , X2=X2 , U=U , Nks=Nks ,
                          Qs=Qs , constraints=loop.const , mq.phi=mq.phi ,
                          tsf=tsf , tsf.ind=tsf.ind , pav.alg=loop.PAV , hp=hp , 
                          em.eps=em.eps , em.iter=em.iter , verbose=verbose[3])
    
    # If global test stat is larger, update current estimate of order
    if( cust.const==FALSE ){
      update.max <- (MNK==1) + (clme.temp$ts.glb > ts.max)
    } else{
      update.max <- 1
    }
    
    if( update.max > 0 ){
      ts.max     <- clme.temp$ts.glb
      clme.out  <- clme.temp
      est.order <- mnk
    }
    
  }
  
  ## Get the observed global and individual test stats
  ts.glb <- clme.out$ts.glb
  ts.ind <- clme.out$ts.ind
  
  if( cust.const==FALSE ){
    grid.row <- list( order=search.grid[est.order,1], 
                      node=search.grid[est.order,3],
                      decreasing=search.grid[est.order,2] ) 
    est.const <- create.constraints( X1=X1 , X2=X2 , constraints=grid.row  ) 
  } else{
    est.const <- constraints
  }
  
  if( nsim > 0 ){
    ## Obtain bootstrap samples
    Y.boot <- resid.boot( Y=Y , X1=X1 , constraints=est.const , X2=X2 , U=U ,
                          Nks=Nks , Qs=Qs , nsim=nsim , mq.phi=mq.phi )
    
    ## EM for the bootstrap samples    
    p.value <- rep( 0 , length(ts.glb) )
    pval.ind <- rep( 0 , dim(est.const$A)[1] )
    
    for( m in 1:nsim ){
  
      if( verbose[1]==TRUE ){
        print( paste( "Bootstrap Iteration " , m , " of " , nsim , sep=""))
      }
      
      ## Loop through the search grid
      ts.boot <- -Inf
      
      for( mnk in 1:MNK ){
        
        if( cust.const==FALSE ){
          grid.row <- list( order=search.grid[mnk,1], node=search.grid[mnk,3],
                            decreasing=search.grid[mnk,2] )
          loop.const <- create.constraints( X1=X1 , X2=X2 , 
                                            constraints=grid.row  )
          
          if( method=="PAVA" & is.function(pav.alg)==FALSE ){
            loop.PAV <- switch( 
              EXPR=search.grid[mnk,1] ,
              "simple"      = pava.simple.order ,
              "simple.tree" = pava.simple.tree ,
              "umbrella"    = pava.umbrella
            )
          }
        }
        
        clme.temp <- clme.em( method=method , Y=Y.boot[,m] , X1=X1 , X2=X2 , 
                              U=U , Nks=Nks , Qs=Qs , constraints=loop.const , 
                              mq.phi=mq.phi , tsf=tsf , tsf.ind=tsf.ind , 
                              pav.alg=loop.PAV , hp=hp, em.eps=em.eps , 
                              em.iter=em.iter , verbose=verbose[3])
                
        
        idx <- (clme.temp$ts.glb > ts.boot)
        if( length(idx)>0 ){
          ts.boot[idx] <- clme.temp$ts.glb[idx]
        }
        
        ## Individual test statistics for the estimated order
        #update.max <- (MNK==1) + (mnk==1) + (clme.temp$ts.glb > ts.max)
        #if( update.max > 0 ){
          
        update.ind <- (MNK==1) + (mnk == est.order)
        if( update.ind>0 ){
          ts.ind.boot <- clme.temp$ts.ind 
        }
        
      }
      
      # null.dist[m]      <- ts.boot
      # null.dist.ind[m,] <- ts.ind.boot
      
      p.value  <- p.value  + 1*( ts.boot    >= ts.glb )
      pval.ind <- pval.ind + 1*(ts.ind.boot >= ts.ind )
          
    }
    clme.out$p.value     <- round( p.value/nsim  , 4 )
    clme.out$p.value.ind <- round( pval.ind/nsim , 4 )  
    
  }
  
  # Add some values to the output object
  
  clme.out$constraints  <- list( A=est.const$A , B=est.const$B )
  class(clme.out)       <- "clme"
  names(clme.out$theta) <- colnames(X)
  clme.out$method       <- method
  
  ### Add the distribution to the output object? Seems unnecessary.
  # clme.out$null.dist     <- null.dist
  # clme.out$null.dist.ind <- null.dist.ind 
  
  
  ## Report the estimated order
  if( cust.const == TRUE ){
    clme.out$est.order <- "Custom order restrictions were specified."
  } else{
    
    if(est.const$decreasing){
      inc.dec <- "decreasing"
    } else{
      inc.dec <- "increasing"
    }
    
    if( MNK==1 ){
      est.order <- "Order was "
    } else{
      est.order <- "Estimated order was "
    }
    
    if( est.const$order=="simple" ){
      clme.out$est.order <- paste( est.order , inc.dec , " simple order." , 
                                   sep="" )
    } else{
      clme.out$est.order <- paste( est.order , inc.dec , " " , 
                                   est.const$order , " order with node=" ,
                                   est.const$node , "." , sep="" )
    }
    
  }
  
  ## Return the output object
  clme.out
  
}
