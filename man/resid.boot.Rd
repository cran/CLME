\name{resid.boot}

\alias{resid.boot}

\title{Obtain Residual Bootstrap}

\description{ Generates bootstrap samples of the data vector. }


\usage{ resid.boot(Y, X1, X2 = NULL, U = NULL, 
        Nks = dim(X1)[1], Qs = dim(U)[2], constraints,
        nsim = 1000, mq.phi = NULL , seed=NULL)
        }

\arguments{
  \item{Y}{ \eqn{N \times 1}{Nx1} vector of response data. }
  \item{X1}{ \eqn{N \times p_1}{Nxp1} design matrix. }
  \item{constraints}{ list of constraints. Attempts to call \code{create.constraints} if element \code{A} is \code{NULL}. }
  \item{X2}{ optional. \eqn{N \times p_2}{Nxp2} matrix of covariates. }
  \item{U}{ optional. \eqn{N \times c}{Nxc} matrix of random effects. }
  \item{Nks}{ optional. \eqn{K \times 1}{Kx1} vector of group sizes. }
  \item{Qs}{ optional. \eqn{Q \times 1}{Qx1} vector of group sizes for random effects.}
  \item{nsim}{ optional. Number of bootstrap samples to generate. }
  \item{mq.phi}{ MINQUE estimates of variance parameters. Automatically estimated if \code{mq.phi=NULL}.}
  \item{seed}{ optional. Set the seed for the RNG. }
}


\value{ Output is \eqn{N \ times nsim}{N x nsim} matrix, where each column is a bootstrap sample of the response data \code{Y}. }

\references{
Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}
}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This function is primarily designed to be called by \code{\link{constrained.lme}}. 

By default, homogeneous variances are assumed for the residuals and (if included) random effects. Heterogeneity can be induced using the arguments \code{Nks} and \code{Qs}, which refer to the vectors \eqn{ (n_{1}, n_{2}, \ldots, n_{k}) }{(n1, n2 ,... , nk)} and \eqn{ (c_{1}, c_{2}, \ldots, c_{q}) }{(c1, c2 ,... , cq)}, respectively. See \code{\link{CLME-package}} for further explanation the model and these values.


}

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
const  <- create.constraints( X1=X1 , X2=X2 ,
          constraints=list(order='simple' , decreasing=FALSE) )

boot.sample <- resid.boot( Y=Y, X1=X1, X2=X2, U=U, nsim=10 , constraints=const  )


}


