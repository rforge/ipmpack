\name{plotGrowthModelComp}
\alias{plotGrowthModelComp}
\alias{plotSurvModelComp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots compared models built with \code{growthModelComp} and \code{survModelComp}.
}
\description{
Function plots compared models built with \code{growthModelComp} and \code{survModelComp}. This can be invoked directly from \code{growthModelComp} and \code{survModelComp} with the argument \code{makePlot = TRUE}.
}
\usage{
plotGrowthModelComp(dataf, treatN, expVars, grObj, testType, summaryTable, plotLegend = TRUE, mainTitle = "")
		
plotSurvModelComp(dataf, treatN, svObj, expVars, testType, summaryTable, plotLegend = TRUE, mainTitle = "") 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grObj/survObj}{
    a list with the objects of the class growth object or survival object equal to \code{treatN}. 
  }
  \item{summaryTable}{
    dataframe output from \code{growthModelComp} and \code{survModelComp} that contains linear predictor and testType scores (see \code{growthModelComp} and \code{survModelComp}).
  }
  \item{dataf}{
    dataframe containing \code{size} and \code{sizeNext} 
  }
  \item{expVars}{
    vector, list of covariates.  Defaults to \code{c("1",  "size", "size + size2")}.
  }
   \item{respType}{
    character string identifying the response variable for the growthModelComp.  wither \code{"sizeNext"}, \code{"logincr"} or \code{"incr"}. Defaults to \code{sizeNext}. 
  }
  \item{plotLegend}{
    logical indicated whether a legend is created. If \code{TRUE}, positions the legend in \code{"topleft"} for growth models and \code{"bottomleft"} for survival models. 
  }
  \item{mainTitle}{
    string to place as the \code{main} attribute in plots (if \code{makePlot = TRUE}.  defaults to \code{NULL}.
  }
 }
\details{
Plots multiple growth and survival objects returned from \code{growthModelComp} and \code{survModelComp}.  See \code{plotGrowthModelComp} and \code{plotSurvModelComp} for more details.  
\value{
a plot object.
}
\references{
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{growthModelComp}},\code{\link{growthModelComp}}
}
\examples{
# Data with size and sizeNext
dff <- generateData()

grModels <- growthModelComp(dff, makePlot = FALSE)
plotGrowthModelComp(grObj = grModels$growthObjects, summaryTable = grModels$summaryTable, dataf = dff, expVars = as.character(grModels$summaryTable[,1]))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
