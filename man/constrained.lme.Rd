\name{constrained.lme}
\alias{constrained.lme}

\title{ Constrained Inference for Linear Mixed Effects Models.}
\description{ Implements constrained inference for linear fixed or mixed effects models using distribution-free bootstrap methodology. }

\usage{
constrained.lme( method="QPE" , Y , X1 , X2=NULL , U=NULL ,       
          Nks=dim(X1)[1] , Qs=dim(U)[2] , constraints=list() ,
          nsim=1000 , em.eps=sqrt(.Machine$double.eps) , em.iter=500, 
          mq.eps=sqrt(.Machine$double.eps) , mq.iter=500 , 
          verbose=c(FALSE,FALSE,FALSE) , tsf=lrt.stat , 
          tsf.ind=w.stat.ind , pav.alg=NULL, hp=FALSE , seed=NULL )
}

\arguments{
  \item{method}{ method for applying order constraints. Supported are \dQuote{QPE} and \dQuote{PAVA} }
  \item{Y}{ \eqn{N \times 1}{Nx1} vector of response data. }
  \item{X1}{ \eqn{N \times p_1}{Nxp1} design matrix. }
  \item{X2}{ optional. \eqn{N \times p_2}{Nxp2} matrix of covariates. }
  \item{U}{ optional. \eqn{N \times c}{Nxc} matrix of random effects. }
  \item{Nks}{ optional. \eqn{K \times 1}{Kx1} vector of group sizes. }
  \item{Qs}{ optional. \eqn{Q \times 1}{Qx1} vector of group sizes for random effects.}
  \item{constraints}{optional. List containing the constraints, defaults to empty list. See Details for further information. }
  \item{nsim}{ optional. Number of bootstrap samples to generate. }
  \item{em.eps}{ optional. Criterion for convergence for the EM algorithm. }
  \item{em.iter}{ optional. Maximum number of iterations permitted for the EM algorithm. }
  \item{mq.eps}{ optional. Criterion for convergence for the MINQUE algorithm. }
  \item{mq.iter}{ optional. Maximum number of iterations permitted for the MINQUE algorithm. }
  \item{verbose}{ optional. Vector of 3 logicals. The first causes printing of iteration step, the second two are passed as the \code{verbose} argument
  to the functions \code{\link{minque}} and \code{\link{clme.em}}, respectively. }
  \item{tsf}{ Function to calculate the global test statistic. Defaults to likelihood ratio test statistic (\code{\link{lrt.stat}}). See Details for further information. }
  \item{tsf.ind}{ Function to calculate the test statistic for individual constrats. Defaults to Williams type test statistic (\code{\link{w.stat.ind}}). See Details for further information. }
  \item{pav.alg}{ algorithm to apply PAVA constraints. Only required when custom \code{constraints} are specified and \code{method=}\dQuote{PAVA}. Ignored when \code{method=}\dQuote{QPE}. See Details for further information. }
  \item{hp}{logical to indicate whether weights for PAVA should be use the full covariance matrix or the diagonal elements. Ignored when \code{method=}\dQuote{QPE}.}
\item{seed}{ set the seed for the RNG.}
}

\details{
If any random effects are included, the function then computes MINQUE estimates of variance parameters. After, \code{\link{clme.em}} is run to obtain the observed values, followed by \code{\link{resid.boot}}. Finally, \code{\link{clme.em}} is run on each of the bootstrap samples and a \eqn{p}{p}-value is computed.
}

\value{
The output of \code{constrained.lme} is an object of the class \code{clme}, which is list with elements:

\item{\code{method}}{ the method that was used to apply order constraints (QPE or PAVA)}
\item{\code{theta}}{ estimates of \eqn{\theta}{theta} coefficients}
\item{\code{ssq}}{ vector of estimate(s) of residual variance, \eqn{\sigma^{2}_{i}}{sigma.i^2}.}
\item{\code{tsq}}{ vector of estimate(s) of random effects variance, \eqn{\tau^{2}_{i}}{tau.i^2}.}
\item{\code{cov.theta}}{ the covariance matrix of \eqn{\theta}{theta}}
\item{\code{ts.glb}}{ test statistic for the global hypothesis.}
\item{\code{ts.ind}}{ vector of test statistics for each of the constraints.}
\item{\code{p.value}}{ p-value for the global hypothesis}
\item{\code{p.value.ind}}{ vector of p-values for each of the constraints}
\item{\code{constraints}}{ list containing the constraints (\code{A}) and the contrast for the global test (\code{B}).}
\item{\code{est.order}}{ sentence describing the estimated order (or whether custom constraints were specified).}
}


\references{

Farnan, L., Ivanova, A., and Peddada, S. D. (2014).
  Linear Mixed Efects Models under Inequality Constraints with Applications.
  \emph{PLOS ONE}, 9(1). e84778. doi: 10.1371/journal.pone.0084778
  \url{http://www.plosone.org/article/info:doi/10.1371/journal.pone.0084778}

}
\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }

\note{

The argument \code{constraints} is a list containing the order restrictions. The elements are \code{order}, \code{node}, \code{decreasing}, \code{A}, and \code{B}, though not all are necessary. The function can calculate the last two for default orders (simple, simple tree, or umbrella). For default orders, \code{constraints} should be a list containing any subset of \code{order}, \code{node}, and \code{descending}. See the figure below for a depiction of these values; the pictured \code{node} of the simple tree orders (middle column) is 1, and the \code{node} for the umbrella orders (right column) is 3. These may be vectors (e.g. order=('simple','umbrella') ). If any of these three are missing, the function will test for all possible values of the missing element(s), excluding simple tree.

For non-default orders, the elements \code{A} and \code{B} should be provided. \code{A} is an \eqn{r \times p_1}{r x p1} matrix (where r is the number of linear constraints, \eqn{0 < r}{0 < r} and \eqn{p_1}{p1} is the dimension of \eqn{\theta_1}{theta1}, i.e. \code{p1=ncol(X1)}). Additionally, \code{B} should hold the contrasts needed for calculating the test statistic. Element \code{B} is only needed if \code{tsf=w.stat}.

Two arguments are functions and may require more attention:
\itemize{
\item \code{tsf} is a function to calculate the desired test statistic. The default function calculates likelihood ratio type test statistic. A Williams type test statistic, which is the maximum of the test statistic over the constraints in \code{constraints\$B}, is also available, and custom functions may be defined. See \code{\link{w.stat}} for details.

\item \code{pav.alg} is a function to implement the Pool Adjacent Violators Algorithm. This is unnecessary if \code{method='QPE'} Default functions are included for simple, simple tree, and umbrella orders (see \code{\link{pava.functions}}). If a custom \code{A} matrix is provided but this argument is not a valid function, the function will revert to using QPE estimation. See \code{\link{pava.functions}} for details on PAVA for defining custom functions.

}

\figure{OrderPlot.jpg}{Plot of Orders.}

A user-friendly display of the output from \code{constrained.lme} can be accessed either by \code{summary} or \code{\link{summary.clme}}.
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

args(constrained.lme)

# Note that 'nsim' has been reduced for illustration

clme.out1 <- constrained.lme( method='QPE', Y=Y, X1=X1, X2=X2, U=U,
                     constraints=const, nsim=10 , tsf=lrt.stat )

clme.out2 <- constrained.lme( method='PAVA', Y=Y, X1=X1, X2=X2, U=U,
                     constraints=const , nsim=10 )


\dontrun{
# Let program estimate the order ('constraints' not specified)
# Takes longer to run
clme.out3 <- constrained.lme( method='PAVA', Y=Y, X1=X1, X2=X2, U=U,
                     nsim=10 , em.eps=0.0001)
}


}


