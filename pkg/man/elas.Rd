\name{elas}
\alias{elas}
\alias{sens}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Functions to estimate elasticity and sensitivity of a matrix
}
\description{
  Functions to estimate elasticity and sensitivity of each element of the IPM
}
\usage{
elas(A)
sens(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    a matrix defining all transitions (growth, survival, fecundity) between sizes
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a matrix
}
\references{
Caswell, 2001, Matrix Population Models: construction, analysis,
interpretation. Sinauer. %add page!
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
Modified following code developed by Mark Rees, Dylan Childs. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}
}
\examples{
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))

Fmatrix <-create.IPM.Fmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),fecObj=makeFecObj(dff))

IPM <- Tmatrix + Fmatrix

elast <- elas(IPM)
image(Tmatrix@meshpoints,Tmatrix@meshpoints,t(elast), 
	main="elasticity", xlab="size", ylab="sizeNext")

senst <- sens(IPM)
image(Tmatrix@meshpoints,Tmatrix@meshpoints,t(senst), 
	main="sensitivity", xlab="size", ylab="sizeNext")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
