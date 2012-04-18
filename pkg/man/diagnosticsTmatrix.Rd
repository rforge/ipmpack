\name{diagnosticsTmatrix}
\alias{diagnosticsTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a series of diagnostic graphs for a T matrix.
}
\description{
  Displays the effects of increasing number of bins and continuous (size) stage range on a number of predictions from the T matrix to verify that sufficient resolution and continuous stage range are being used. 
}
\usage{
diagnosticsTmatrix(Tmatrix, growObj, survObj, dff, integrateType, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tmatrix}{
    an IPMmatrix object.
}
  \item{growObj}{
    the growth object used to construct the IPMmatrix object.
}
  \item{survObj}{
    the survival object used to construct the IPMmatrix object.
}
  \item{dff}{
    the dataframe from which the survival and growth objects were constructed.
}
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density function); other option is "cumul" (which uses the cumulative density function).
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for survival at that size. 
  }
}
\details{
  Top left panel: continuous (size) stage range in the data (histogram) and the size range used in the IPM (red dashed lines).
  Top center panel: survival (y axis) plotted against the sums of columns of the current.
  T matrix (black) or an extended T matrix (red) where minSize is set to x 0.5 its current value value (or x 1.5 its current value, if it is negative); the maxSize set to x 1.5 its current value, and the number of bins set to 1.5 the current value; the dashed line is the 1:1 ratio.
  Top right panel: predicted life expectancy for the current (black) and extended IPM as above (red).
  Bottom panels: numerical predicted density function (histogram) and the expected value from the normal density function (red line) for (x 0.5) small, (x 1) medium, and (x 1.5) large sizes.
  
  Not defined for compound T matrices, or matrices with discrete classes.
  
  Note that if the pdfs corresponding to large sizes are lower than the histogram, this indicates
  that the maximum size is such that some fraction of the pdf is being lost.  
  
  
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81, p694-708.

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167, p410-428.

For effects of mesh size on IPM output: Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createIPMTmatrix}}
}
\examples{
# Example where mesh size does not have a major effect on model output:
dff <- generateData()
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm= TRUE),
	maxSize = max(dff$size, na.rm=TRUE), growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff)

# Compare with the following example where mesh size has an important effect on output:
Tmatrix <- createIPMTmatrix(nBigMatrix = 8, minSize = min(dff$size, na.rm= TRUE), 
	maxSize = 0.5*max(dff$size, na.rm=TRUE), growObj = makeGrowthObj(dff), 
		survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff)

#with cumul
Tmatrix <- createIPMTmatrix(nBigMatrix = 10, minSize = min(dff$size, na.rm= TRUE), 
	maxSize = max(dff$size, na.rm=TRUE), growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), 
	integrateType="cumul")
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff, integrateType="cumul")

#with truncated increment (fit is quite bad, lots of loss from integration)
gr <- makeGrowthObj(dff)
gr1 <- new("growthObjTruncIncr")
gr1@fit$coefficients <- gr@fit$coefficients
gr1@fit$sigmax2 <- (summary(gr@fit)$sigma)^2
Tmatrix <- createIPMTmatrix(nBigMatrix = 30, minSize = min(dff$size, na.rm= TRUE), 
	maxSize = max(dff$size, na.rm=TRUE), growObj = gr1, 
		survObj = makeSurvObj(dff), integrateType="midpoint")
diagnosticsTmatrix(Tmatrix, growObj = gr1, survObj = makeSurvObj(dff), dff = dff)

#with log increment
Tmatrix <- createIPMTmatrix(nBigMatrix = 50, minSize = min(dff$size, na.rm= TRUE), 
	maxSize = max(dff$size, na.rm=TRUE), growObj = makeGrowthObj(dff,responseType = "logincr"), 
		survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff,responseType = "logincr"), 
	survObj = makeSurvObj(dff), dff = dff)


}
