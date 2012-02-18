\name{MeanLifeExpect}
\alias{MeanLifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to calculate mean life expectancy
}
\description{
  Provided with a Tmatrix, defining survival transitions across size,
  provides a vector defining life-expectancy in units of the time-step used, for each of the size bins. 
}
\usage{
MeanLifeExpect(IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object defining survival transitions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a vector of life-expectancies each corresponding to a value in Tmatrix@meshpoints
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. Sinauer. 

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LifeExpect}}, \code{\link{create.IPM.Tmatrix}} 
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObjGeneral(dff),
survObj=makeSurvObjGeneral(dff))
PassageTime(5,Tmatrix)

Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObjGeneral(dff),
survObj=makeSurvObjGeneral(dff))
MeanLifeExpect(Tmatrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
