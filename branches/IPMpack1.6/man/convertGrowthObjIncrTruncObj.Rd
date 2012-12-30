\name{convertGrowthObjIncrTruncObj}
\alias{convertGrowthObjIncrTruncObj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert an increment based growth object to a truncated increment based growth object.
}
\description{
Function that takes a growth object fitted with a linear regression with increment as the response variable, and returns the corresponding truncated increment growth object (useful if, for example, shrinkage is not thought to occur) 
}
\usage{
convertGrowthObjIncrTruncObj(growthObj)
		
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growthObj}{
	the growth object for which the conversion is required
  }
}
\value{
 a growth object of class truncated increment}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\details{
Note that fitting a truncated normal cannot be done by fitting a normal distribution and cutting off the tails; truncation redistributes the probability density between the size limits which modifies size transition probabilities.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeGrowthObj}},\code{\link{makeSurvObj}},\code{\link{plotGrowthModelComp}}, \code{\link{plotSurvModelComp}}
}
\examples{
# Data with size and sizeNext
dff <- generateData()

gr1 <- makeGrowthObj(dff,incr~size)
gr2 <- convertGrowthObjIncrTruncObj(gr1)


}
