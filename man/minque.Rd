\name{minque}
\alias{minque}
\title{ MINQUE Algorithm }
\description{ Implements MINQUE algorithm to obtain estimates of variance parameters of a linear mixed effects model.}

\usage{
minque( Y , X1 , X2 = NULL , U = NULL , Nks = dim(X1)[1] ,
        Qs = dim(U)[2], mq.eps = sqrt(.Machine$double.eps) ,
        mq.iter = 500 , verbose = FALSE )
}

\arguments{
  \item{Y}{ \eqn{N \times 1}{N x 1}  vector of response data. }
  \item{X1}{ \eqn{N \times p_1}{N x p1} design matrix. }
  \item{X2}{ optional. \eqn{N \times p_2}{N x p2} matrix of covariates. }
  \item{U}{ optional. \eqn{N \times c}{N x c}  matrix of random effects. }
  \item{Nks}{ optional. \eqn{K \times 1}{K x 1} vector of group sizes. See Details. }
  \item{Qs}{ optional. \eqn{Q \times 1}{Q x 1} vector of group sizes for random effects.}
  \item{mq.eps}{ criterion for convergence for the MINQUE algorithm. }
  \item{mq.iter}{ maximum number of iterations permitted for the MINQUE algorithm. }
  \item{verbose}{ if \code{TRUE}, function prints intermediate messages on progress of the MINQUE algorithm. }
}

\details{
By default, the model assumes homogeneity of variances for both the residuals and the random effects (if included). See the Details in \code{\link{CLME-package}} for more information on how to use the arguments \code{Nks} and \code{Qs} to permit heterogeneous variances. 
}

\value{

The function returns a vector of the form \eqn{(\tau^{2}_{1}, \tau^{2}_{2}, \ldots, \tau^{2}_{q}, \sigma^{2}_{1},\sigma^{2}_{2},\ldots, \sigma^{2}_{k})'}{(tau1^2, tau2^2, \ldots, tauq^2, sigma1^2,sigma2^2,\ldots, sigmak^2)'}. If there are no random effects, then the output is just \eqn{(\sigma^{2}_{1},\sigma^{2}_{2},\ldots, \sigma^{2}_{k})'}{(sigma1^2,sigma2^2,\ldots, sigmak^2)'}.

}


\references{

Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}

}

\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{
This function is called by several other function in \pkg{CLME} to obtain estimates of the random effect variances. If there are no random effects, they will not call \code{minque}.
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

# Assume homogeneity of variances and no covariates or random effects
minque( Y=Y , X1=X1 )

# Assume homogeneity of variances and include covariates and random effects
minque( Y=Y , X1=X1 , X2=X2 , U=U )

# Assume heterogeneity of variances, and include random effects
minque( Y=Y , X1=X1 , X2=X2 , U=U , Nks=Nks , Qs=Qs )

}


