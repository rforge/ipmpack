\name{varPassageTime}
\alias{varPassageTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates variation in passage time.
}
\description{
Function to take a T matrix (either compound or not) and estimate variance in passage time to a chosen continuous stage value.
}
\usage{
varPassageTime(chosen.size, IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chosen.size}{
The continuous stage value of interest.
}
  \item{IPM.matrix}{
The Tmatrix (compound or not).
}
}
\details{
Note how variation in passage time for values exactly equal to the chosen size (targetSize) are low, because of way the conditionals are framed. Passage time from values > than targetSize should be ignored (space to the right of the red vertical line in example below), unless dealing with an organism that is able to display retrogression.
}
}
\value{
Numeric vector corresponding to variance in passage time from each of the meshpoints in the IPM (so both size, and if a compound matrix, size from different environments). 
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p119.

Metcalf, Horvitz, Tuljapurkar & Clark. 2009. A time to grow and a time to die: a new way to analyze the dynamics of size, light, age and death of tropical trees. Ecology 90, p2766-2778.

For bias in this estimation where variance in growth is small relative to the size range: Zuidema, Jongejans, Chien, During & Schieving. 2010. Integral Projection Models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.

For species with shrinkage: Salguero-Gomez & Casper. 2010. Keeping shrinkage in the demographic loop. Journal of Ecology 98, p313-323.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PassageTime}}, \code{\link{create.IPM.Tmatrix}} 
}
\examples{
# With continuous (e.g. size) stage
dff <- generateData()
Tmatrix <- create.IPM.Tmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize <- max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff))
targetSize <- 8
vP <- varPassageTime(targetSize, Tmatrix)

plot(Tmatrix@meshpoints, vP, type = "l", xlab="Continuous (e.g. Size) stage", ylab = "Variance in passage time", col = "dark gray")
abline(v = targetSize, col = "red")

# With a compound matrix (discrete stages not yet defined)
Tmatrix <- create.compound.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), envMatrix = makeEnvObj(dff),
growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"),
survObj = makeSurvObj(dff, explanatoryVariables = "size+size2+covariate"), discreteTrans = 1)
targetSize <- 8
vP <- varPassageTime(targetSize, Tmatrix)

plot(vP, ylab = "Passage time", xlab = "Continuous (e.g. size) stage", type = "l", col = "dark gray", ylim = c(0, max(vP)))
abline(v = targetSize, col = "red")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
