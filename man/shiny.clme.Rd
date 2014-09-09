\name{shiny.clme}

\alias{shiny.clme}
\alias{shiny.server}
\alias{shiny.ui}

\title{ Shiny GUI for CLME }
\description{ A graphical user interface to run \pkg{CLME}, built from the \pkg{shiny} package. }

\usage{
shiny.clme()
shiny.server(input,output)
shiny.ui()
}

\arguments{
  \item{input}{ input from GUI.}
  \item{output}{ output to GUI.}
}


\details{
Currently the GUI does not allow specification of custom orders for the alternative hypothesis. Future versions will enable this capability.

The data file should be a CSV file with the first row being a header. The format should be: 

\deqn{\left[ Y : X_{1} : X_{2} : U \right]}{ [ Y : X1 : X2 : U ] },

where ':' denotes column-binding (i.e., the first column is Y, the next p1 columns are the columns of X1, etc). These must be formatted appropriately. For example, \eqn{X_{1}}{X1} should not be a factor, but rather a matrix where each column is an indicator for membership in a factor level.

Currently the user interface requires manual input for the number of columns of \eqn{X_{1}}{X1}, \eqn{X_{2}}{X2}, and \eqn{U}{U}. If there are no covariates or random effects, set the number of columns for those terms to zero (default).

}



\author{ Casey M. Jelsema <casey.jelsema@nih.gov> }



\examples{

\dontrun{ shiny.clme() }

}





