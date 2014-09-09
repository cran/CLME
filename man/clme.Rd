\name{clme}

\alias{clme}
\alias{as.clme}
\alias{is.clme}


\title{Constructor method for objects S3 class clme}
\description{
Creates an object of class clme. Methods \code{as.clme} and \code{is.clme} coerce (with error checking) to class clme and test whether object is of class clme.
}

\usage{
clme( x , theta=numeric(0), ssq=numeric(0), tsq=numeric(0), 
                  cov.theta=matrix(numeric(0)), ts.glb=numeric(0), 
                  ts.ind=numeric(0), p.value=numeric(0), 
                  p.value.ind=numeric(0),  constraints=list(), 
                  method=character(0), est.order=character(0) )
\method{as}{clme}(x, ... )
\method{is}{clme}( x )
}

\arguments{
  \item{x}{ list with the elements corresponding to the output of \code{\link{constrained.lme}}. }
  \item{theta}{vector of coefficients or effects.}
  \item{ssq}{scaler or vector of residual variance estimates.}
  \item{tsq}{scaler or vector of random effect variance estimates.}
  \item{cov.theta}{covariance matrix of theta.}
  \item{ts.glb}{global test statistic(s).}
  \item{ts.ind}{test statistics for individual constraints.}
  \item{p.value}{global p-value.}
  \item{p.value.ind}{p-values for individual constraints.}
  \item{constraints}{list containing information on constrants. Required element is A, the matrix of linear constraints.}
  \item{method}{method of isotonization (not required).}
  \item{est.order}{sentence stating the type of ordering (not required).}  
  \item{...}{space for additional arguments.}  
}


\value{
Returns an object of the class \code{clme}.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\seealso{
\code{\link{CLME-package}}, \code{\link{constrained.lme}}
}

\examples{


clme1 <- clme()

clme2 <- list( theta=c(1,3,4.5,9) )
is.clme( clme2 )
as.clme(clme2)

clme3 <- list( theta=c(1,3,4.5,9) , ssq=1 , tsq=2 , 
               cov.theta=as.matrix( c(1,2,3,4)\%*\%t(c(1,2,3,4)) ) ,
               ts.glb=1 , ts.ind=c(1,2,3) , p.value=0.5 , 
               p.value.ind=c(0.2,0.5,0.8) ,
               constraints=list( A=matrix(0,nrow=3)))
is.clme( clme3 )
as.clme( clme3 )

}

