\name{growthModelComp}
\alias{growthModelComp}
\alias{survModelComp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compares growth and survival objects built from different covariate sets.
}
\description{
Function compares model fits for growth and survival objects built with different linear combinations of covariates (plotting currently restricted to transforms of size; but comparison can include any chosen covariates).  Growth can have multiple response forms. 
Returns a list containing a summary table of covariates and scores, and another list containing all of the growth (or survival) objects used in the comparison.}
\usage{
growthModelComp(dataf, expVars = c("1", "size", "size + size2"), 
	regressionType = "constantVar", respType = "sizeNext", testType = "AIC", 
		makePlot = FALSE, mainTitle = "",...)
		
survModelComp(dataf, expVars = c("1", "size", "size + size2"), testType = "AIC", 
	makePlot = FALSE, mainTitle = "",ncuts=20,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
 	dataframe with columns size, surv, and the growth response variable of choice
   }
  \item{expVars}{
    vector, list of covariates.  Defaults to \code{c("1",  "size", "size + size2")}.
  }
  \item{regressionType}{
    character string identifying whether the type of regression run will have constant or changing variance (for \code{growthModelComp}.  Defaults to \code{constantVar}.
  }
  \item{respType}{
    character string identifying the response variable for the growthModelComp.  wither \code{"sizeNext"}, \code{"logincr"} or \code{"incr"}. Defaults to \code{sizeNext}. 
  }
  \item{testType}{
    character string identifying the metric used to compare models.  Can be any string that uses \code{loglike} from the lm or glm object.  For example \code{"AIC"} or \code{"BIC"}.  Defaults to \code{"AIC"}.
  }
  \item{makePlot}{
    logical whether to make plots with the comparison building. Defaults to \code{FALSE}.
  }
  \item{mainTitle}{
    string to place as the \code{main} attribute in plots (if \code{makePlot = TRUE}.  defaults to \code{NULL}.
  }
  \item{ncuts}{
 for survModelComp, number of cuts in the data-set to be used in plotting
  }
  \item{...}{
 	additional arguments for plotting (ylim, col, etc)
   }
 }
\details{
Both \code{growthModelComp} and \code{survModelComp} use a dataframe that has variables \code{size} and \code{sizeNext} to build a series of nested models.  The default will build growth or survival objects with an intercept, an intercept and size, an an intercept with size and size^2 terms. 
%Covariates can also be added, and require that these variables have the term "cov" in the variable name. This is not TRUE? right? 

The models build use only \code{lm} or \code{glm} (and not \code{mcmcGLMM} for example) to estimate maximum likelihood estimates of functions.  The testType (default \code{"AIC"} uses the \code{loglike} output from the lm or glm objects to score the model.  

Plotting calls the functions \code{plotGrowthModelComp} or \code{plotSurvModelComp} to plot the objects.  These functions can also be called after building the model comparison lists that are returned.  If called outside of the initial building functions, they need to receive the \code{GrowthObjects} or \code{SurvObjects} list in the outputList from the build function.  See \code{plotGrowthModelComp} and \code{plotSurvModelComp} for more details.  
}
\value{
  a list with a summary table of class \code{dataframe} with models and scores and list of containing the objects of class grObj and survObj for each model.
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

growthModelComp(dff, makePlot = TRUE)
survModelComp(dff, makePlot = TRUE)

growthModelComp(dff, makePlot = TRUE, regressionType = "declineVar")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
