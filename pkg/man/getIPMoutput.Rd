\name{getIPMoutput}
\alias{getIPMoutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets IPM output from a list of P matrices (only survival and size information) or full IPMs (P matrices + F matrices; the latter include sexual reproduction information).
}
\description{
Gets synthetic values including life expectancy, passage time, and, if a fecundity matrix is available, population growth rate (lambda), stable stage distribution, reproductive output, etc. 
}
\usage{
getIPMoutput(PmatrixList, targetSize = c(), FmatrixList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PmatrixList}{
    List of survival-size (T) matrices for which summary statistics desired.
}
  \item{targetSize}{
    Target size for passage time. If none is provided defaults to the median of the IPM meshpoints.
}
  \item{FmatrixList}{
    List of fecundity (F) matrices corresponding to the list of P matrices, facultative. When this information is included, population growth rate (lambda), and stable stage distribution will also be provided.
}
}
\value{
  \item{LE}{matrix of life expectancies, columns correspond to
  meshpoints, rows corresponding to each element of the list of P matrices}
  \item{pTime}{matrix of passage times to the targetSize from each of
    the meshpoints (columns) and for each element in the P matrix list
    (columns).}
  \item{lambda}{vector of population growth rates corresponding to value obtained
    combining each element of the list of P matrices with the
    corresponding element in the list of F matrices; if no F matrix list
    is provided, it returns a vector of NAs.}
  \item{stableStage}{matrix of stable stage distributions rows corresponding to values obtained.
    combining each element of the list of Pmatrices with the
    corresponding element in the list of Fmatrices; if no Fmatrix list
    is provided, this is a matrix of NAs.}
}
\references{
Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{getIPMoutputDirect}} 
}
\examples{
dff <- generateData()

## TODO: Cory example

}
