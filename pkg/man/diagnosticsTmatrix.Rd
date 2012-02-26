\name{diagnosticsTmatrix}
\alias{diagnosticsTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to create a series of diagnostic graphs for a T matrix
}
\description{
  Displays the effects of increasing number of bins and size range on a
  number of predictions from the T matrix to verify that sufficient
  resolution and size range are being used. 
}
\usage{
diagnosticsTmatrix(Tmatrix, growObj, survObj, dff, integrateType, correction)
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
  \item{integrateType}{
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
  Top left panel: size range in the data (histogram) and the size range used in the IPM (red dashed lines)
  Top center panel: survival (y axis) plotted against the sums of columns of the current
  T matrix (black) or an extended T matrix (red) where minSize is set to x 0.5 its current
  value (or x 1.5 its current value, if it is negative); the max size set to x 1.5 its current value, and the number of bins set to 1.5 the current value; the dashed line is the
  1:1 ratio
  Top right panel: predicted life expectancy for the current (black) and extended IPM as above (red)
  Bottom panels: numerical predicted density function (histogram) and the expected
  value from the normal density function (red line) for (x 0.5) small, (x 1) medium,
  and (x 1.5) large sizes.
  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81 p694-708

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167 p410-428

For effects of mesh size on IPM output: Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98 p345-355

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
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff), dff=dff)

#compare with the following example where mesh size has an important effect on output:
Tmatrix <-create.IPM.Tmatrix(nBigMatrix=8,minSize=min(dff$size,na.rm=TRUE),
maxSize=0.5*max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff), survObj=makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff), dff=dff)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
