\name{plot.clme}

\alias{plot.clme}

\title{ S3 method to plot objects of class \code{clme}}

\description{ Generates a basic plot of estimated \eqn{\theta}{theta} coefficients which are subject to constraints. Lines indicate individual constraints (not global tests) and significance.}

\usage{
\method{plot}{clme}( x , alpha=0.05 , place.leg="below" , inset=0.01,
            ci.wd=0, plot.names=TRUE, ylim, cex=1.75, pch=21, bg="white" ,  
            xlab = expression(paste("Component of ", theta[1])),
            ylab = expression(paste("Estimated Value of " , theta[1])), 
           ...)
}

\arguments{
  \item{x}{ object of class 'clme'.}
  \item{alpha}{ significance level of the test.}
  \item{place.leg}{ character indicating placement of legend. See Details.}
  \item{inset}{ inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword. }
  \item{ci.wd}{ width of bars for the confidence limits. }
  \item{plot.names}{ logical, whether to plot the group names on the X-axis.}
  \item{ylim}{ limits of the y axis.}
  \item{cex}{ size of plotting symbols.}
  \item{pch}{ plotting symbols.}
  \item{bg}{ background (fill) color of the plotting symbols.}
  \item{xlab}{ label of the x axis.}
  \item{ylab}{ label of the y axis.}
  \item{\dots}{ additional plotting arguments. }
}

\details{
All of the individual contrasts in the \code{constraints\$A} matrix are tested and plotted. The global contrasts from \code{constraints\$B} are not shown (unless they happen to align with the individual contrasts). Only the elements of \eqn{\theta}{theta} which appear in any constraints (e.g. the elements of \eqn{\theta_{1}}{theta1}) are plotted. Coefficients for the covariates are not plotted.

Solid lines denote no significant difference, while dashed lines denote statistical significance. Significance is determined by the individual p-value being less than or equal to the supplied \eqn{\alpha}{alpha} threshold. By default a legend denoting the meaning of solid and dashed lines will be placed below the graph. Argument \code{place.leg} may be set to a legend keyword (e.g. \code{legend=''bottomright''}) to place it inside the graph at the specified location. Setting \code{place.leg} to \code{FALSE} or to a non-supported keyword suppresses the legend.

Confidence intervals for the coefficients may be plotted, but the user is cautioned about these. They are individual confidence intervals, and are computed using the covariance matrix of the unconstrained estimates of \eqn{\theta_{1}}{theta1}. As a result, the confidence intervals will be wider and sometimes appear to conflict with the significance tests.
}


\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }



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

plot.clme( clme.out )



}