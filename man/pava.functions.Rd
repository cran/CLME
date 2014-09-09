\name{pava.functions}
\alias{pava.functions}
\alias{pava.simple.order}
\alias{pava.simple.tree}
\alias{pava.umbrella}



\title{ Apply PAVA Order Constraints }
\description{ Implements a Pool Adjacent Violators Algorithm (PAVA). }
\usage{
pava.simple.order(theta, cov.theta=diag(length(theta)) , node = NULL,
                  decreasing = FALSE , hp = FALSE)
pava.simple.tree(theta, cov.theta=diag(length(theta)), node = 1,
                 decreasing = FALSE , hp = FALSE )
pava.umbrella(theta,cov.theta=diag(length(theta)), node = 1, 
              decreasing = FALSE , hp = FALSE )
}

\arguments{
  \item{theta}{ \eqn{p_1 \times 1}{p1 x 1} vector on which to apply PAVA order constraints.}
  \item{cov.theta}{ covariance matrix of \eqn{\theta}{theta} }
  \item{node}{ node element. Set to \code{NULL} for simple order. }
  \item{decreasing}{ logical indicating whether the initial order constraint is increasing or decreasing.}
  \item{hp}{logical indicating whether to use the full covariance matrix (\code{TRUE}) or just the diagonal elements (\code{FALSE}).}
}

\details{
Refer to \code{\link{create.constraints}} or \code{\link{constrained.lme}} for more information on arguments \code{node} and \code{decreasing}.
}

\value{
The function returns a \eqn{p_1 \times 1}{p1 x 1} vector of order-restricted estimates of \code{theta}.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
For \code{hp=FALSE}, the functions \code{pava.simple.order} and \code{pava.umbrella} act as wrappers for function \code{pava} from package \code{Iso}. For \code{hp=TRUE}, the full covariance matrix of \eqn{\theta}{theta} is used. 

Custom PAVA algorithms (supplied to other functions such as \code{\link{constrained.lme}} or \code{\link{clme.em}} as the argument \code{pava.algorithm}) must take the same four arguments and produce the same output as described here. Often such functions may call the default PAVA functions in a particular way (specific to the application at hand). 
}


\examples{

test.pava <- c(4,2,7,5,9)
wts       <- diag(length(test.pava))

pava.simple.order( test.pava , wts, NULL, FALSE, FALSE)
pava.simple.tree(  test.pava , wts, 3   , TRUE , FALSE)
pava.umbrella(     test.pava , wts, 2   , TRUE , FALSE)

}


