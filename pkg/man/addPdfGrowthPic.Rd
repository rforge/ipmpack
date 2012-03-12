\name{addPdfGrowthPic}
\alias{addPdfGrowthPic}
\alias{survModelComp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds probability density functions of density function of size or increment given current size and growth to plots.
}
\description{
Function generates pdfs (probability density functions) corresponding to chosen sizes and adds them to a figure using growth methods.}
\usage{
addPdfGrowthPic(respType = "sizeNext", sizesPlotAt = c(20, 50, 60),
                            sizeRange = c(20, 400), incrRange = c(-10, 50), 
                            scalar = 100, growthObjList,
                            cols = 1:5,
                            cov = 1,
                            minShow = 1e-2,
                            jitt = 2,  
                            ...)
		
survModelComp <- function(dataf, expVars = c("1", "size", "size + size2"), testType = "AIC", makePlot = FALSE, mainTitle = "",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{respType}{
    character string identifying the response variable for the growthModelComp.  wither \code{"sizeNext"}, \code{"logincr"} or \code{"incr"}. Defaults to \code{sizeNext}. 
  }
  \item{sizesPlotAt}{
    vector, list of sizes at which pdfs should be plotted.
  }
  \item{sizeRange}{
    sizeRange for which pdf should be estimated
  }
  \item{incrRange}{
    increment range for which pdf should be estimated
  }
  \item{scalar}{
    value by which pdf may be multiplied to improve visibility
  }
  \item{growthObjList}{
    list of growth objects for which pdfs are desired to be plotted
  }
  \item{cols}{
    colours corresponding to the list of growth objects for plotting
  }
  \item{cov}{
    covariate value (if growth objects contain a discrete covariate, or
    other); defaults to 1; will be ignored if no covariates are fitted
  }
  \item{minShow}{
    minimum value below which pdf lines will not be shown (to avoid ugly
    vertical lines)
  }
  \item{jitt}{
    amount by which sequential pdfs should be separated on the x axis, for
    visibility 
  }
  \item{...}{
    extra arguments to pass to the plot function.
  }
}
\details{
}
\value{
  none
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
\code{\link{makeGrowthObj}},\code{\link{makeSurvObj}},\code{\link{plotGrowthModelComp}}, \code{\link{plotSurvModelComp}}
}
\examples{
# Data with size and sizeNext
dff <- generateData()

a1 <- growthModelComp(dff, makePlot = TRUE)
addPdfGrowthPic(respType = "sizeNext", sizesPlotAt = c(2, 6, 10), scalar = 1, jitt = 0.1, sizeRange = c(-5, 25), growthObjList = a1$growthObjects, cols = 2:5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
