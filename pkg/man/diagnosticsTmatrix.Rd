\name{diagnosticsTmatrix}
\alias{diagnosticsTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create a series of diagnostic graphs for a Tmatrix
}
\description{
  Displays the effects of increasing number of bins and size range on a
  number of predictions from the Tmatrix to verify that sufficient
  resolution and size range are being used. 
}
\usage{
diagnosticsTmatrix(Tmatrix, growObj, survObj, dff, integrate.type, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
    an IPM.matrix object
}
  \item{growObj}{
    the growth object used to construct the IPM.matrix object
}
  \item{survObj}{
    the survival object used to construct the IPM.matrix object    
}
  \item{dff}{
    the dataframe from which the survival and growth objects were constructed, with the usual headings
}
  \item{integrate.type}{
    integration type, defaults to "midpoint" (which uses prob. density
    function); other option is "cumul" (which uses the cumul. density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
}
\details{
  The first panel shows the size range in the
  data (histogram) and the size range used in the IPM (red dashed
  lines); the second panel shows survival (y axis) plotted against the sums of columns of the current
  Tmatrix (black) or an extended Tmatrix (red) where minSize set to 0.5 its current
  value (or 1.5 times its current value, if it is negative), the max size set to 1.5 its current value, and the number of bins set to 1.5 the current value; the dashed line is the
  1:1 line; the third panel shows predicted life expectancy for the current
  (black) and extended IPM as above (red); the last three plots
  show numerical predicted density function (histogram) and the expected
  value from the normal density function (red line) for small, medium,
  and large sizes. 
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
\code{\link{create.IPM.Tmatrix}}
}
\examples{
dff<-generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff), dff=dff)

#compare with...
Tmatrix <-create.IPM.Tmatrix(n.big.matrix=8,minsize=min(dff$size,na.rm=TRUE),
maxsize=0.5*max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff), dff=dff)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
