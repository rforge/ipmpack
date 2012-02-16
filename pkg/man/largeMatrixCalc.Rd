\name{largeMatrixCalc}
\alias{largeMatrixCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate lambda and stable size distribution efficiently
}
\description{
Method to calculate lambda and stable size distrbution where a large
number of bins are used
}
\usage{
largeMatrixCalc(Tmatrix, Fmatrix, tol = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
object of class IPM.matrix describing survival transitions
}
  \item{Fmatrix}{
object of class IPM.matrix describing fertility transitions
}
  \item{tol}{
    tolerance for convergence, defaults to 1e-08
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{lambda}{Population rate of increase}
  \item{stable.dist}{Stable size distribution}
  \item{h1}{size bin width}
}
\references{
From Rees & Ellner, 2009, American Naturalist. 
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObjGeneral(dff),
survObj=makeSurvObjGeneral(dff))
Fmatrix <-create.IPM.Fmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),fecObj=makeFecObjGeneral(dff))

largeMatrixCalc(Tmatrix,Fmatrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
