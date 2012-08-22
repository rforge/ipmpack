\name{makeDiscreteTrans}
\alias{makeDiscreteTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds matrix for transitions between discrete (e.g. seedbank) and continuous (e.g. size) stages.
}
\description{
  Function that takes a data frame reflecting a population that includes continuous (e.g., size) and discrete (e.g., diapause)
  classes of individuals, and builds a transition matrix from this for movement between discrete and continuous stages (providing just a single value for the continuous stages; the detail of movement between continuous stages (e.g., sizes) is generated elsewhere). This object can then be used as an argument in the function to "createIPMPmatrix" to build a P matrix that contains both discrete and continuous stages. 

  The dataframe must contain columns  "stage" and "stageNext" with values of the names of the discrete classes and the name "continuous" where appropriate, in the current and subsequent time step. Continuous categories are also defined by their measurements, contained in "size" and "sizeNext". This is necessary for defining the mean and variance of the continuous stage of individuals emerging from discrete stages.
}
\usage{
makeDiscreteTrans(dataf, continuousToDiscreteExplanatoryVariables = "size")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns "stage", "stageNext" (both FACTORS (use \code{as.factor}) containing either names of discrete stages or "continuous" or "dead"), "size", "sizeNext", and "surv" (continuous and boolean variables respectively).}
  \item{continuousToDiscreteExplanatoryVariables}{a character defining the relationship defining size influences the probability of individuals in the continuous class moving to one of the discrete classes.}
}
\value{
  an object of class "discreteTrans" with columns corresponding to all the discrete and the one continuous stage.
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81, p694-708.

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167, p410-428.

%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createIPMPmatrix}}

}
\examples{
dff <- generateDataDiscrete()
discTrans <- makeDiscreteTrans(dff)
createIPMPmatrix(nBigMatrix = 10, growObj = makeGrowthObj(dff), 
	survObj = makeSurvObj(dff), discreteTrans = discTrans)
}
