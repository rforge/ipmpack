\name{contourPlot}
\alias{contourPlot}
\alias{sens}
\alias{elas}
%Comment out when I make the function LTRE available
%\alias{LTRE}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produces a graph with the smooth kernel and associated value label.
}
\description{
  Outputs a graph with a smoothed surface of the IPM kernel or other chosen matrix. This graph also contains a label for values of changes on stage conditional on survival and per-capita sexual/clonal contributions. This function is applicable to the different components of the IPM separately (i.e. Pmatrix, Fmatrix, C matrix) or to the kernel of the sensitivity and elasticity IPM. This function can also be applied with lifecycles with discrete stages (e.g. seedbank in "hyperDataCovSubset") but the user must be careful with the interpretation of perturbation analyses between discrete and continuous stages since they are logically on a different scale.
}
\usage{
contourPlot(M, meshpts, upper, lower, color)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
    a matrix of the type P matrix, F matrix, C matrix, the full IPM (= P + F + C), its sensitivity or elasticity. The function transposes the matrix internally so the user does not have to worry about having size in time t as the x axis and size in time t+1 as the y axis.
  \item{meshpts}{
    meshpoints of the kernel to be plotted.
  \item{upper}{
    the highest value to be displayed in the "z" axis. If plotting the full IPM, it is recommended to set this value carefully so that the P matrix (which is constrained between 0 and 1) is not visually "lost" due to potentially high value of the F/C matrices (e.g., individuals may contribute to the population sexually/clonally with potentially very high values of recruits). Values of the kernel above the upper threshold color scheme default to white color.
  \item{lower}{
    the lowest value to be displayed in the "z" axis. If plotting P, F, IPM, mesh-level sensitivity/elasticity values, it is recommended to set this value to 0 since those kernels are always non-negative. If using this function to display mesh-level LTRE differences, negative mesh-cells might be produced (e.g. depicting a negative effect of a treatment IPM relative to the control IPM in a specific region of the IPM, like in recruitment for instance) and thus the user is recommended to set the lowest negative value. Values of the kernel above the lower threshold color scheme default to white color
  \item{color}{
    color scheme with which to display the z values of the kernel. Options are extracted directly from the palette options available via the library "graphics", in function "filled.contour" and library "grDevices", in "palette". Options are: rainbow, heat.colors, terrain.colors, topo.colors and cm.colors.
}

}
\value{
  a graph.
}
\references{
For more details on plotting contour surfaces: the Graphics R package. http://stat.ethz.ch/R-manual/R-patched/library/graphics/html/00Index.html

On futher literature for LTREs: Caswell, 2001, Matrix Population Models: construction, analysis,
interpretation. 2nd ed. Sinauer. p206-256.

de Kroon, Plaisier, van Groenendael & Caswell. 1986. Elasticity: the relative contribution of demographic parameters to population growth rate. Ecology 67, p1427-1431.

de Kroon, van Groenendael & Ehrlen. 2000. Elasticities: a review of methods and model limitations. Ecology 81, p607-618.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}, \code{\link{elas}}
}
\examples{
dff <- generateData()
Pmatrix <- createIPMPmatrix(minSize = min(dff$size, na.rm=TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff),
survObj = makeSurvObj(dff))

Fmatrix <- createIPMFmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), fecObj = makeFecObj(dff))

IPM <- Pmatrix + Fmatrix

#Plot P matrix contour plot
minSize <- 0
maxSize <- max(Pmatrix)
contourPlot(Pmatrix,Pmatrix@meshpoints,minSize,maxSize,heat.colors)

#Plot F matrix contour plot
minSize <- 0
maxSize <- max(Fmatrix)
contourPlot(Fmatrix,Fmatrix@meshpoints,minSize,maxSize,heat.colors)

#Plot IPM contour plot
contourPlot(IPM,Pmatrix@meshpoints,minSize,maxSize,heat.colors)

#Plot IPM sensitivity contour plot
senIPM <- sens(IPM)
minSize <- 0
maxSize <- max(senIPM)
contourPlot(senIPM,Pmatrix@meshpoints,minSize,maxSize,heat.colors)

#Plot IPM elasticity contour plot
elaIPM <- elas(IPM)
minSize <- 0
maxSize <- max(elaIPM)
contourPlot(senIPM,Pmatrix@meshpoints,minSize,maxSize,heat.colors)

}
