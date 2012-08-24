\name{makeDiscreteTrans}
\alias{makeDiscreteTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds matrix for transitions between discrete (e.g. seedbank, dormant) and continuous (e.g. size) stages.
}
\description{
  Function that takes a data frame reflecting a population that includes continuous (e.g., size) and discrete (e.g., diapause)
  classes of individuals, and builds a transition matrix from this for movement between discrete and continuous stages (providing just a single value for the continuous stages; the detail of movement between continuous stages (e.g., sizes) is generated elsewhere). This object can then be used as an argument in the function to "createIPMPmatrix" to build a P matrix that contains both discrete and continuous stages. 

  The dataframe must contain columns  "stage" and "stageNext" with values of the names of the discrete classes and the name "continuous" where appropriate, in the current and subsequent time step. Continuous categories are also defined by their measurements, contained in "size" and "sizeNext". This is necessary for defining the mean and variance of the continuous stage of individuals emerging from discrete stages. Alternatively, you can enter the transition matrix for the discrete stages in the 'discreteTrans' argument.
}
\usage{
makeDiscreteTrans(dataf, stages = NA, discreteTrans = NA, meanToCont = NA, sdToCont = NA, continuousToDiscreteExplanatoryVariables = "size")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns "stage", "stageNext" (both FACTORS (use \code{as.factor}) containing either names of discrete stages or "continuous" or "dead"), "size", "sizeNext", and "surv" (continuous and boolean variables respectively).}
  \item{stages}{a character vector with the names of the discrete classes. Normally this argument does not have to be used as the names of discrete classes are extracted from the data or entered discreteTrans matrix.}
  \item{discreteTrans}{a matrix with transition probabilities between the discrete and continuous stages. The column names should match the stage/stageNext names in the data file. The names discrete stages should be in alphabetical order, followed by 'continuous'. The row names should be the same, but with 'dead' added at the bottom. Thus, this matrix represents all the fates of individuals from the various classes.}
  \item{meanToCont}{a matrix containing the mean sizeNext values for individuals moving from discrete classes to the continuous classes (should contain NA when no individuals move from a particular discrete class to a continuous sizeNext). The column names should be equal (in alphabetical order) to the names of the discrete classes (so no continuous).}
  \item{sdToCont}{a matrix containing the sd sizeNext values for individuals moving from discrete classes to the continuous classes (should contain NA when no individuals move from a particular discrete class to a continuous sizeNext). The column names should be equal (in alphabetical order) to the names of the discrete classes (so no continuous).}
  \item{continuousToDiscreteExplanatoryVariables}{a character defining the relationship defining size influences the probability of individuals in the continuous class moving to any of the discrete classes.}
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
	
# Often some data needed for makeDiscreteTrans are not available as individual records in your data file. For instance the fate of seeds in the seedbank of the generateDataDiscrete() example:
dffNew <- subset(dff,!dff$stage%in%c("seedAge1","seedOld"))
# Now the transition matrix needs to be entered as an argument, as well as the mean and sd of sizeNext values for discrete to continuous transitions. 
discTransNew <- makeDiscreteTrans(dffNew,
					discreteTrans = matrix(c(.26,0,0,.74,0,0,0,.44,.07,.48,0,0,.43,.11,.45,.04,0,0,.43,.54),
					nrow=5, ncol=4, dimnames = list(c("dormant","seedAge1","seedOld","continuous","dead"),
					c("dormant","seedAge1","seedOld","continuous"))),
					meanToCont = matrix(c(mean(dffNew$sizeNext[dffNew$stage=="dormant"&dffNew$stageNext=="continuous"],na.rm=TRUE),3.4,2.4),
					nrow=1, ncol=3, dimnames = list(1,c("dormant","seedAge1","seedOld"))), 
					sdToCont = matrix(c(sd(dffNew$sizeNext[dffNew$stage=="dormant"&dffNew$stageNext=="continuous"],na.rm=TRUE),2.2,1.8),
					nrow=1, ncol=3, dimnames = list(1,c("dormant","seedAge1","seedOld"))))
					 


	
}
