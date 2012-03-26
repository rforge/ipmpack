\name{largeMatrixCalc}
\alias{largeMatrixCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates population growth rate (lambda) and stable stage distribution in a computationally efficient way when the number of bins in the IPM is large.
}
\description{
Method to calculate population growth rate (lambda) and stable stage distribution where a large number of bins are used in the IPM, as it may be the case with species that vary largely in size, or models that include size x age interactions.
}
\usage{
largeMatrixCalc(Tmatrix, Fmatrix, tol = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
	object of class IPMmatrix describing survival transitions.
}
  \item{Fmatrix}{
	object of class IPMmatrix describing fecundity transitions.
}
  \item{tol}{
    tolerance for convergence, defaults to 1e-08.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{lambda}{Population rate of increase.}
  \item{stableDist}{Stable stage distribution.}
  \item{h1}{size bin width.}
}
\references{
Rees and Ellner. 2009. Integral projection models for populations in temporally varying environments. The American Naturalist 79, p575-594.

Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p377-502.

Garcia, Dahlgren, Ehrlen. 2011. No evidence of senescence in a 300-year-old mountain herb. Journal of Ecology 99, p1424-1430.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
Modified from Appendix A in Rees and Ellner 2009 (see references).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dff <- generateData()
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff))
Fmatrix <- createIPMFmatrix(minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE),fecObj = makeFecObj(dff))

largeMatrixCalc(Tmatrix, Fmatrix)

plot(largeMatrixCalc(Tmatrix, Fmatrix)$stableDist, ylab = "Stable stage distribution", 
	xlab = "Continuous (e.g. size) stage", type = "l", col = "blue", lty = 1, ylim = c(0:1))

%ROB TO DO: I will change the appropriate section in base so reproductive output is also attainable and plottable here for IPMpack v.2.0
%#Add line for reproductive effort with dash line and red color
%#points(largeMatrixCalc(Tmatrix,Fmatrix)$reprodValue-0.1,type="l",col="red",lty=2)
%legend(1, 0.9,c("Stable distribution","Reproductive value"),col=c("blue","red"),lty=c(1,2),merge=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
