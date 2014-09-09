\name{w.stat}

\alias{w.stat}
\alias{w.stat.ind}


\title{ Williams Type Test Statistic. }
\description{ Calculates a Williams type test statistic for a constrained linear mixed effects model. This is the default test statistic for \pkg{CLME}. }

\usage{
w.stat( theta , cov.theta , B , A , ...)
w.stat.ind( theta , cov.theta , B , A , ...)
}

\arguments{
  \item{theta}{ estimated coefficients.}
  \item{cov.theta}{ covariance matrix of the (unconstrained) coefficients. }
  \item{B}{   matrix to obtain the global contrast. See Details. }
  \item{A}{   matrix of linear constraints.}
  \item{\dots}{ additional arguments, to enable custom test statistic functions.}
}


\details{
Not all arguments are used by 

If \eqn{B}{B} is a matrix  with each row \eqn{B_i}{Bi} being a contrast, and \eqn{\Sigma}{Sig} is the covariance matrix of \eqn{\theta}{theta}, \code{w.stat} returns the scaler:

\deqn{W = max \frac{B_i \theta}{ B_{i} \Sigma^{-1} B_{i}^{'} }}{ W = max [Bi * theta] / [ Bi * Sigma^-1 * Bi' }

and \code{w.stat.ind} returns the vector W where:

\deqn{W_{i} = \frac{A_i \theta}{ A_{i} \Sigma^{-1} A_{i}^{'} }}{ Wi = max [Ai * theta] / [ Ai * Sigma^-1 * Ai' }

}



\value{ Output is a numeric scaler or vector. }
\references{

Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}

}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This is an internal function, unlikely to be useful outside of \link{CLME-package}. To define custom functions, the arguments available are:

\code{theta}, \code{cov.theta}, \code{B}, \code{Y}, \code{X1}, \code{X2}, \code{U}, \code{tsq}, , \code{ssq}, \code{Nks}, and \code{Qs}.

These are equivalently named as the arguments to \code{\link{constrained.lme}}. The new arguments are \code{cov.theta}, the covariance matrix of the unrestricted estimators of \eqn{\theta}{theta}, \code{ssq}, the estimate of \eqn{ (\sigma_{1}^{2},\ldots,\sigma_{k}^{2})'}{(sigma1^2,\ldots,sigmak^2)'} and \code{tsq}, the estimate of \eqn{ (\tau_{1}^{2},\ldots,\tau_{q}^{2})'}{(tau1^2,\ldots,tauq^2)'}. 

Custom functions must output a scaler or a vector. Vector output corresponds to testing multiple hypotheses.

}


\examples{

theta  <- exp(1:4/4)
th.cov <- diag(4)
X1     <- matrix( 0 , nrow=1 , ncol=4 )
const  <- create.constraints( X1=X1 , constraints=list(order='simple' ,
          decreasing=FALSE) )

w.stat( theta , th.cov , const$B , const$A )

w.stat.ind( theta , th.cov , const$B , const$A )

}






