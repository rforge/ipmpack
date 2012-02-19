\name{largeMatrixCalc}
\alias{largeMatrixCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate population growth rate (lambda) and stable size distribution in a computationally efficient way when the number of bins in the IPM is large
}
\description{
Method to calculate lambda and stable size distribution where a large
number of bins are used in the IPM, as it may be the case with species that vary largely in size, or models that include size x age interactions.
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
	object of class IPM.matrix describing fecundity transitions
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
Rees & Ellner. 2009. Integral projection models for populations in temporally varying environments. The American Naturalist 79 p575-594

Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p377-502

Garcia, Dahlgren & Ehrlen. 2011. No evidence of senescence in a 300-year-old mountain herb. Journal of Ecology 99 p1424Ð1430

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
Modified from Appendix A in Rees & Ellner 2009 (see references)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
Fmatrix <-create.IPM.Fmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),fecObj=makeFecObj(dff))

largeMatrixCalc(Tmatrix,Fmatrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
