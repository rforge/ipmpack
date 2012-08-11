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
diagnosticsTmatrix(Tmatrix, growObj, survObj, dff, integrateType, correction, cov = data.frame(covariate = 1))
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
  integration type, defaults to ``midpoint" (which uses probability density function); other option is ``cumul" (which uses the cumulative density function).
 }
 \item{correction}{
  correction type, defaults to ``none"; option is ``constant" which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for survival at that size. 
 }
 \item{cov}{
  a data-frame with one row containing all covariates; defaults to 1, and will be ignored if covariates do not exist in growth and survival objects
}
 \item{sizesToPlot}{
 a vector containing desired sizes to plot growth resolution for (second panel); if not supplied, the function will use the quantiles
}
}
\details{
IPMpack includes a useful function {\tt diagnosticsTmatrix} that provides a
series of plots indicative of whether bin choice and size range is adequate.
Applying this function as a preliminary step before obtaining demographic and evolutionary output from IPMs can help identify basic problems in the fitting of vital rate functions or in the creation of the IPM matrices before proceeding. Several common problems can be diagnosed with the panelled figure produced by {\tt diagnosticsTmatrix} (Figure
\ref{fig:diag}).

The top-left panel shows the range of the data and the range of the state
variable fitted in the matrix.  If these are mis-matched, the limits of the
data used in building the Tmatrix can be adjusted with the {\tt minSize} and
{\tt maxSize} arguments in {\tt createIPMTmatrix}.  One common problem
can result from the loss of parts of the continuous distribution
when binning. This can mean that the sum of the columns of the matrix do not
match the fitted survival.  The top-middle panel of Figure \ref{fig:diag} will
show this discrepency if the red and black lines do not overlay.  To correct
this, {\tt createIPMTmatrix} has an argument {\tt correction = "constant"}
which ensures that the columns sum to the fitted survival by
multiplying every column in the Integral
Projection Model by the value that will ensure this. The top-right panel checks whether extending the size range included in the
matrix and increasing the number of bins (by increasing {\tt
  nBigMatrix} and thereby having narrower bins) does not alter basic predictions from the IPM.
The three lower panels show the discretized IPM (histograms) and the
theoretical density function.  If the theoretical density function
curve is very distant from the histograms, increasing the {\tt nBigMatrix}
argument may correct this discrepency. (Note that {\tt
  diagnosticsTmatrix} plots a range of three standard deviations around the starting value for the
final three plots, so in some cases, where extreme changes in size are
possible via growth, the histogram may be off the
$x$ axis scale). 
 
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
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm = TRUE),
	maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff)

# Compare with the following example where mesh size has an important effect on output:
Tmatrix <- createIPMTmatrix(nBigMatrix = 8, minSize = min(dff$size, na.rm = TRUE), 
	maxSize = 0.5*max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), 
		survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff)

#with cumul
Tmatrix <- createIPMTmatrix(nBigMatrix = 10, minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), 
	integrateType = "cumul")
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), dff = dff, integrateType = "cumul")

#with truncated increment (fit is quite bad, lots of loss from integration)
gr <- makeGrowthObj(dff)
gr1 <- new("growthObjTruncIncr")
gr1@fit$coefficients <- gr@fit$coefficients
gr1@fit$sigmax2 <- (summary(gr@fit)$sigma)^2
Tmatrix <- createIPMTmatrix(nBigMatrix = 30, minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE), growObj = gr1, 
		survObj = makeSurvObj(dff), integrateType = "midpoint")
diagnosticsTmatrix(Tmatrix, growObj = gr1, survObj = makeSurvObj(dff), dff = dff)

#with log increment
Tmatrix <- createIPMTmatrix(nBigMatrix = 50, minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff,Formula = logincr~size), 
		survObj = makeSurvObj(dff))
diagnosticsTmatrix(Tmatrix, growObj = makeGrowthObj(dff,Formula = logincr~size), 
	survObj = makeSurvObj(dff), dff = dff)


}
