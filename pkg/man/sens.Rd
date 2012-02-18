\name{sens}
\alias{sens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Function to estimate element by element sensitivity
 }
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sens(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{elas}}
}
\examples{
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff))
ss <- sens(Tmatrix)

image(Tmatrix@meshpoints, Tmatrix@meshpoints,t(ss), 
	xlab="size", ylab="sizeNext")
contour(Tmatrix@meshpoints, Tmatrix@meshpoints,t(ss),add=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
