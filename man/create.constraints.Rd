
\name{create.constraints}
\alias{create.constraints}
\title{ Generate common order constraints }
\description{ Automatically generates the constraints in the format used by functions \code{\link{constrained.lme}} or \code{\link{clme.em}}. Allowed orders are simple, simple tree, and umbrella orders.  }

\usage{ create.constraints(X1, X2 = NULL, constraints) }

\arguments{
  \item{X1}{ \eqn{N \times p_1}{Nxp1} design matrix. }
  \item{X2}{ optional. \eqn{N \times p_2}{Nxp2} matrix of covariates. }
  \item{constraints}{ List with the elements \code{order}, \code{node}, and \code{decreasing}. See Details for further information. }
}

\details{

The elements of \code{constraints} are:
\itemize{
\item \code{order}: string. Currently \dQuote{simple}, \dQuote{simple.tree} and \dQuote{umbrella} are supported.
\item \code{node}: numeric scaler, the node of the coefficients.
\item \code{decreasing}: logical. Defines if the initial order relation is decreasing (e.g. \eqn{H_{a}: \theta_{1} \geq \theta_{2} \geq \ldots \geq \theta_{p_{1}} }{Ha: theta_1 >= theta_2 >= ... >= theta_p1}).
}

See \code{\link{constrained.lme}} for more information and a depiction of these three elements.

}


\value{
The function returns a list containing the elements of input argument \code{constraints} as well as
\item{ \code{A} }{matrix of dimension \eqn{r \times p_1}{r x p1} containing the order constraints, where r is the number of linear constraints.}
\item{ \code{B} }{matrix containing the contrasts necessary for computation of the test statistic (may be identical to \code{A}).}

See \code{\link{w.stat}} for more information on \code{B}

}



\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
The functions \code{\link{constrained.lme}} and \code{\link{clme.em}} both utilize the argument \code{constraints}. 

For \code{\link{clme.em}}, this argument is intended to be the output of \code{create.constraints}. It must contain \code{A}, \code{B}, and if \code{method=}\dQuote{PAVA}, values of \code{node} or \code{decreasing} as necessary.

For \code{\link{constrained.lme}} it may either be the output of \code{create.constraints}, or equivalent to the argument \code{constraints} of that function.
}
 

\seealso{
\code{\link{constrained.lme}},
\code{\link{clme.em}}
}

\examples{
n  <- 5
P1 <- 5
X1 <- diag(P1) \%x\% rep(1,n)
X2 <- as.matrix( rep(1,P1) \%x\% runif(n , 0,2) )

# For simple order, the node does not matter
create.constraints( X1=X1 , constraints=list( order='simple' ,
                    decreasing=FALSE ))

# See effect of covariates, compare constraints against decreasing=TRUE
create.constraints( X1=X1 , X2=X2 , constraints=list( order='simple' , 
                    decreasing=TRUE ))

# Umbrella order
create.constraints( X1=X1 , constraints=list( order='umbrella' , node=3
                    , decreasing=FALSE ))

}


