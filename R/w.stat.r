

w.stat <- function(  theta , cov.theta , B , A , ...  ){
  
  cov.contrast <- c( diag( B %*% cov.theta %*% t(B) ) )
  diffs        <- c( B%*%theta )
  W.to.return  <- max( diffs/sqrt(cov.contrast) )
  W.to.return
  
}




w.stat.ind <- function(  theta , cov.theta , B , A , ...  ){
  
  cov.contrast <- c( diag( A %*% cov.theta %*% t(A) ) )
  diffs        <- c( A%*%theta )    
  W.to.return  <- diffs/sqrt(cov.contrast)
  W.to.return
}
