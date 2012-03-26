\name{getListRegObjects}
\alias{getListRegObjects}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates growth objects reflecting distribution of parameters from lm or glm.
}
\description{
Function generates pdfs corresponding to estimated parameters and the
variance covarianace matrix that defines them using a multivariate
normal distribution and then defines the corresponding growth or
survival objects.}
\usage{
getListRegObjects(Obj,nsamp=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Obj}{
    a growth or survival or fertility object with a slot named "fit"
    containing an lm or glm, etc
  }
  \item{nsamp}{
    desired number of samples from the posterior multivariate normal
  }

}
\details{
}
\value{
  list of objects containing an lm or glm
}
\references{
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeGrowthObj}},\code{\link{makeSurvObj}},\code{\link{getIPMoutputDirect}}
}
\examples{
# Data with size and sizeNext
dff <- generateData()
gr1 <- makeGrowthObj(dff)
sv1 <- makeSurvObj(dff)

grList <- getListRegObjects(gr1, nsamp=10)
svList <- getListRegObjects(sv1, nsamp=10)


res <- getIPMoutputDirect(survObjList = svList,
                               growObjList = grList,
                               targetSize = 6,
                               nBigMatrix = 100,
                               minSize = 1.1*min(dff$size,na.rm=TRUE),
                               maxSize = 1.1*max(dff$size,na.rm=TRUE),
                               integrateType="midpoint",correction="constant", 
                               storePar=FALSE)


  
     # Plot out with different colours for different rows:
     par(mfrow = c(2,2), bty = "l", pty = "s")
     plot(res$meshpoints, res$LE[1, ], xlab = "Continuous (e.g. size) stage", 
     	ylab = "Life expectancy", type= "l", ylim = range(res$LE, na.rm=TRUE))
     for (j in 1:nrow(res$LE)) points(res$meshpoints,res$LE[j, ],col=j, type = "l")

     plot(res$meshpoints, res$pTime[1, ], xlab="Continuous (e.g. Size) stage", 
     	ylab = "Passage time", type = "l", ylim = range(res$pTime, na.rm=TRUE))
     for (j in 1:nrow(res$pTime)) points(res$meshpoints,res$pTime[j, ],col=j,  type = "l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
