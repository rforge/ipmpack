\name{getListRegObjects}
\alias{getListRegObjects}
\alias{getListRegObjectsFec}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates growth objects reflecting distribution of parameters from lm or glm.
}
\description{
Function generates pdfs corresponding to estimated parameters and the
variance covariance matrix that defines them using a multivariate
normal distribution and then defines the corresponding growth or
survival objects.}
\usage{
getListRegObjects(Obj,nsamp=1000)
getListRegObjectsFec(Obj,nsamp=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Obj}{
    a growth or survival object with a slot named "fit" for \code{getListRegObjects} containing an lm or glm, etc or a fertility object with a slot named "fitFec" for \code{getListRegObjectsFec}, likewise. 
   
  }
  \item{nsamp}{
    desired number of samples from the posterior multivariate normal
  }

}
\value{
  list of growth or survival objects containing an lm or glm; or fertility objects likewise
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
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
fv1 <- makeFecObj(dff,Transform="log")

#higher nsamp are advisable for a detailed analysis
grList <- getListRegObjects(gr1, nsamp=10)
svList <- getListRegObjects(sv1, nsamp=10)
fvList <- getListRegObjectsFec(fv1, nsamp=10)


res <- getIPMoutputDirect(survObjList = svList,
                               growObjList = grList,
                               fecObjList = fvList,
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
     
     plot(res$meshpoints, Re(res$stableStage[1, ]), xlab="Continuous (e.g. Size) stage", 
     	ylab = "Stable distribution", type = "l", ylim = range(Re(res$stableStage), na.rm=TRUE))
     for (j in 1:nrow(res$stableStage)) points(res$meshpoints,Re(res$stableStage[j, ]),
     	col=j,  type = "l")
     	
     hist(res$lambda,col="grey", xlab=expression(lambda), main="")
     	
     
     
     
}
