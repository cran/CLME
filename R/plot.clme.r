



"plot.clme" <- 
  function( x , alpha=0.05 , place.leg="below" , inset=0.01,
            ci.wd=0 , plot.names=TRUE , ylim , cex=1.75 , pch=21 , bg="white" , 
            xlab = expression( paste( "Component of " , theta[1] ) ),
            ylab = expression( paste( "Estimated Value of " , theta[1] ) ) , 
           ...){
  
  if( is.clme(x)==FALSE ){
    stop("Argument 'x' is not of class clme.")
  } else{
    
    theta <- c(x$theta)
    A     <- x$constraints$A
    r     <- dim(A)[1]
    p1    <- max( apply( A , 1 , FUN=function(x){ max(which(x != 0)) } ) )
    
    if( is.null(ci.wd) ){ ci.wd <- min( 1/p1 , 1/20 )  }
    
    if( place.leg=="below" ){
      layout( rbind(1,2) , heights=c(7,1) )
    }
    
    # Calculate the confidence intervals    
    theta.lcl <- x$theta - qnorm(1-alpha/2)*sqrt(diag(x$cov.theta) )
    theta.ucl <- x$theta + qnorm(1-alpha/2)*sqrt(diag(x$cov.theta) )
    
    # Pick some reasonable plot limits
    if( missing(ylim) ){
      if( ci.wd>0 ){
        ylim <- c( min(theta.lcl[1:p1]) , max(theta.ucl[1:p1]))
      } else{
        if( min(x$theta[1:p1]) < 0 ){ 
          ymin <- min(x$theta[1:p1])*1.05 
        } else{
          ymin <- min(x$theta[1:p1])/1.05
        }
        if( max(x$theta[1:p1]) < 0 ){ 
          ymax <- max(x$theta[1:p1])/1.05
        } else{
          ymax <- max(x$theta[1:p1])*1.05
        }
        ylim <- c( ymin , ymax )
      }
    }
    
    # The initial plot of the points
    plot( 1:p1 , theta[1:p1] , cex=cex , pch=pch , bg=bg ,
        ylim = ylim , xaxt='n' , xlab = xlab , ylab = ylab , ... )
    
    
    if( plot.names==TRUE & is.null(names(x$theta)[1:p1]) ){
      warning( "Names of theta are null, reverting to indices." )
      plot.names <- FALSE
    }
    if( plot.names==TRUE ){
      theta.names <- names(x$theta)[1:p1]
      axis(side=1, at=1:p1, labels=theta.names , ...)
    } else{
      axis(side=1, at=1:p1, labels=1:p1 , ...)
    }
      
    # Connect the contrasts with solid/dashed lines
    for( i in 1:r){  
      idx <- which( A[i,] !=0 )
      if( x$p.value.ind[i]  > alpha ){ lty <- 1 }
      if( x$p.value.ind[i] <= alpha ){ lty <- 2 }
      points( idx , theta[idx] , lty=lty , lwd=2 , type="l")
    }
    
    ## Add the CIs if necessary
    if( ci.wd>0 ){
      for( i in 1:p1){        
        points( c(i,i) , c(theta.lcl[i],theta.ucl[i]) , type="l"  )
        points( c(i-ci.wd,i+ci.wd) , c(theta.lcl[i],theta.lcl[i]) , type="l"  )
        points( c(i-ci.wd,i+ci.wd) , c(theta.ucl[i],theta.ucl[i]) , type="l"  )        
      }    
    }
    
    # Replot the pointsso the circles are filled
    points( 1:p1 , theta[1:p1] , cex=cex , pch=pch , bg=bg )
    
    ## Put a legend on the plot if requested
    if( place.leg=="below" ){
      safe.mar <- par( no.readonly=TRUE )$mar
      par(mar=c(0, 0, 0, 0) , ...)
      plot.new()
      legend('center','groups', c( paste("p >" , alpha, "  " ) , paste("p <" , alpha, "  " )),
             lty = c(1,2), col=1 , ncol=2 , bty ="o" , ...)
      
      par( mar=safe.mar )
    } else{
      leg.texts <- c("bottom", "bottomleft", "left", "topleft",
                    "top", "topright", "right", "bottomright", "center")
      if( place.leg %in% leg.texts){
        legend( place.leg , legend=c("Non-significant    " , "Significant"), 
                lty = c(1,2), col=1 , ncol=1 , bty ="o" , inset=inset , ...)
      }
    }
    
  }
  
}


