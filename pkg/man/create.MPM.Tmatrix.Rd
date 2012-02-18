\name{create.MPM.Tmatrix}
\alias{create.MPM.Tmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build discrete T matrix
}
\description{
A method that uses bins to construct a discrete survival transition matrix of class IPM using the
data structures used by IPMpack
}
\usage{
create.MPM.Tmatrix(dataf,bins,n.env=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
	a dataframe with columns `size',
    `sizenext' (`size' is size at t, `sizenext' is size at t+1); facultatively, dataf may include `covariate'
    and `covariatel' for discrete covariates, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. }
\item{bins}{
  limits of the desired size classes
}
\item{n.env}{
  currently not used
}
}
\details{
  The number of individuals moving from one bin to the next in one
  time-step is identified and probabilities entered into a matrix. To
  get a classic Matrix Population model this can be added to values
  defined by create.MPM.Fmatrix. Note that the assumption that all
  individuals for whom is.na(sizenext) died.   
}
\value{
  an object of class IPM.matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.IPM.Tmatrix}}, \code{\link{create.MPM.Fmatrix}}
}
\examples{
dff<-generateData()
Tmatrix <-create.MPM.Tmatrix(dff,bins=c(1,4,10,15))
Tmatrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
