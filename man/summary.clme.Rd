\name{summary.clme}
\alias{summary.clme}

\title{ S3 method to summarize results for objects of class \code{clme} }
\description{ Summarizes the output of objects of class \code{clme}, such as those produced by \code{\link{constrained.lme}}. Prints a tabulated display of global and individual tests, as well as parameter estimates.}


\usage{
\method{summary}{clme}(object, alpha=0, dec.theta=2, dec.w=3, dec.var=4, dec.p=4, ...)
}

\arguments{
  \item{object}{an object of class \code{clme}. See \code{\link{constrained.lme}} or \code{\link{as.clme}} }
  \item{alpha}{level of significance (produced confidence intervals for parameter estimates).}
  \item{dec.theta}{ number of decimal places to round estimates of \eqn{\theta}{theta} coefficients. }
  \item{dec.w}{ number of decimal places to round test statistic values. }
  \item{dec.var}{ number of decimal places to round estimates of \eqn{\sigma^{2}}{sigma^2} and \eqn{\tau^{2}}{tau^2}. }
  \item{dec.p}{ number of decimal places to round p-values. }
  \item{\dots}{ additional arguments (not currently implemented). }
}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }


\note{ The individual tests are performed on the specified order. If no specific order was specified, then the individual tests are performed on the estimated order. }



\seealso{
\code{\link{CLME-package}}, \code{\link{constrained.lme}}
}


\examples{

set.seed( 42 )

n  <- 5
P1 <- 5

X1 <- diag(P1) \%x\% rep(1,n)
X2 <- as.matrix( rep(1,P1) \%x\% runif(n , 0,2) )
U  <- rep(1,P1) \%x\% diag(n)
X  <- as.matrix( cbind(X1,X2) )

tsq <- 1
ssq <- 0.7

Nks <- dim(X1)[1]
Qs  <- dim(U)[2]

xi <- rnorm( sum(Qs)  , 0 , rep(sqrt(tsq) , Qs)  )
ep <- rnorm( sum(Nks) , 0 , rep(sqrt(ssq) , Nks) )  

thetas <- c(2 , 3 , 3, 3 , 4 , 2 )
Y      <- X\%*\%thetas + U\%*\%xi + ep
const  <- list(order='simple' , decreasing=FALSE)


# Note that 'nsim' has been reduced for illustration

clme.out <- constrained.lme( method='PAVA', Y=Y, X1=X1, X2=X2, U=U,
                     constraints=const , nsim=10 )

summary.clme( clme.out )

}