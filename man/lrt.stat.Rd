\name{lrt.stat}

\alias{lrt.stat}


\title{ Likelihood Ratio Type Test Statistic. }
\description{ Calculates the likeihood ratio type test statistic (under Normality assumption) for a constrained linear mixed effects model. This is the default test statistic for \pkg{CLME}. }

\usage{
lrt.stat( theta , A , Y , X1 , X2 , U , tsq , ssq , Nks , Qs , ... )
}

\arguments{
  \item{theta}{ estimated coefficients.}
  \item{A}{   matrix of linear constraints.}
  \item{Y}{   \eqn{N \times 1}{Nx1} vector of response data. }
  \item{X1}{  \eqn{N \times p_1}{Nxp1} design matrix. }
  \item{X2}{  \eqn{N \times p_2}{Nxp2} matrix of covariates. }
  \item{U}{   \eqn{N \times c}{Nxc} matrix of random effects. }
  \item{tsq}{ estimate of random effect variance (may be vector).}
  \item{ssq}{ estimate of residual variance (may be vector).}
  \item{Nks}{ \eqn{K \times 1}{Kx1} vector of group sizes. }
  \item{Qs}{  \eqn{Q \times 1}{Qx1} vector of group sizes for random effects.}
  \item{\dots}{ additional arguments, to enable custom test statistic functions.}
}


\value{ Output is a numeric scalar. }


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This is an internal function, unlikely to be useful outside of \link{CLME-package}. To define custom functions, the arguments available are:

\code{theta}, \code{cov.theta}, \code{B}, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, , \code{ssq}, \code{Nks}, and \code{Qs}.

These are equivalently named as the arguments to \code{\link{constrained.lme}}. The new arguments are \code{cov.theta}, the covariance matrix of the unrestricted estimators of \eqn{\theta}{theta}, \code{ssq}, the estimate of \eqn{ (\sigma_{1}^{2},\ldots,\sigma_{k}^{2})'}{(sigma1^2,\ldots,sigmak^2)'} and \code{tsq}, the estimate of \eqn{ (\tau_{1}^{2},\ldots,\tau_{q}^{2})'}{(tau1^2,\ldots,tauq^2)'}. 

Custom functions must output a scaler or a vector. Vector output corresponds to testing multiple global hypotheses.

}

\seealso{
\code{\link{w.stat}}, 
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
const  <- create.constraints( X1=X1 , X2=X2 , 
          constraints=list(order='simple' , decreasing=FALSE) )

# Estimate parameters (also computed lrt statistic: ts.glb)
clme.out <- clme.em.all( method='QPE', Y=Y, X1=X1, X2=X2, U=U,
                     constraints=const, tsf=lrt.stat )

# Individually compute lrt statistic
lrt.stat( clme.out$theta , const$A , Y , X1 , X2 , 
          U , clme.out$tsq , clme.out$ssq , Nks , Qs )

}



