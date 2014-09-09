\name{rat.blood}

\docType{data}

\alias{rat.blood}

\title{ Experiment on mice}

\description{
  This data set contains the data from an experiment on 24 Sprague-Dawley rats from Cora et al (2012).
}

\usage{rat.blood}
\format{A matrix containing 242 observations on 39 variables.

 \tabular{rl}{
    [,1]   \tab ID for rat (factor). \cr
    [,2]   \tab time period (in order, 0 , 6, 24, 48, 72, 96 hours). \cr
    [,3]   \tab storage temperature reference (\code{''Ref''}) vs. room temperature (\code{''RT''}) \cr
    [,4]   \tab sex, male (\code{''Male''}) vs. female (\code{''Female''}). Coded as \code{''Female''=1}. \cr
    [,5]   \tab response variable \code{wbc}. \cr
    [,6]   \tab response variable \code{rbc}. \cr
    [,7]   \tab response variable \code{hgb}. \cr
    [,8]   \tab response variable \code{hct}. \cr
    [,9]   \tab response variable \code{spun}. \cr
    [,10]   \tab response variable \code{mcv}. \cr
    [,11]   \tab response variable \code{mch}. \cr
    [,12]   \tab response variable \code{mchc}. \cr
    [,13]   \tab response variable \code{plts}. \cr
  }
}

\references{
Cora M, King D, Betz L, Wilson R, and Travlos G (2012). Artifactual changes in Sprauge-Dawley rat hematologic parameters after storage of samples at 3 C and 21 C. Journal of the American Association for Laboratory Animal Science, 51(5), 616-621. URL \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3447451/}.
 
}

\seealso{
\code{\link{CLME-package}}, \code{\link{constrained.lme}}
}

\keyword{datasets}