\name{createMPMTmatrix}
\alias{createMPMTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds a discrete T matrix from the types of data structures used in IPMpack.
}
\description{
A method that uses bins to construct a discrete survival matrix population model of class IPM using the data structures used by IPMpack.
}
\usage{
createMPMTmatrix(dataf, bins, nEnv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
	a dataframe with columns `size', `sizeNext' (`size' is size at t, `sizeNext' is size at t+1); facultatively, dataf may include `covariate' and `covariatel' for discrete covariates, indicating values at t, and at t+1, respectively; these must take values of sequential integers, starting at `1'.
	}
\item{bins}{
  numeric, limits of the desired discrete (size) classes.
}
\item{nEnv}{
  currently not used.
}
}
\details{
  The number of individuals moving from one bin to the next in one time step is identified and probabilities entered into a matrix cell. To obtain a classic matrix population model this can be added to values defined by createMPMFmatrix. Note that the assumption that all
  individuals for whom is.na(sizeNext) died.   
}
\value{
  an object of class IPMmatrix.
}
\references{
For information on T matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

For information on collapsing large matrices: Salguero-Gomez & Plotkin. 2010. Matrix dimensions bias demographic inferences: implications for comparative plant demography. The American Naturalist 176, p710-720.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createIPMTmatrix}}, \code{\link{createMPMFmatrix}}
}
\examples{
dff <- generateData()
matrixCuts <- c(1, 4, 10, 15)
Tmatrix <- createMPMTmatrix(dff, bins = matrixCuts)
Tmatrix
}
