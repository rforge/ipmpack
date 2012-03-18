\name{SilwoodDataClean}
\alias{SilwoodDataClean}
\docType{data}
\title{
Silwood Monocarp Dataset
}
\description{
Demographic measurements from several monocarpic plants from Silwood
}
\usage{data(SilwoodDataClean)}
\format{
  The format is:
 chr "SilwoodDataClean"
}
\details{
data-frame with headings exactDate (date measurement) exactDateNext
(date next measurement) id (unique plant id) Species (species name) Site
(location within Silwood park where plant was measured) rtcr (root crown
diameter, measured with callipers) rtcrNext (root crown diameter at next
census time) ll (length of longest leaf) llNext (length of longest leaf
at next census) rosetteDiam (rosette Diameter) rosetteDiamNext (rosette
diameter at next census) flowered (0 for not flowered, 1 if did) surv (0
if died, one if survived, NAs if not know). 
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(SilwoodDataClean)
print(head(SilwoodDataClean))
plot(SilwoodDataClean$rtcr,SilwoodDataClean$rtcrNext,xlab="size now", ylab="size next", pch=19,log="xy")
## maybe str(SilwoodDataClean) ; plot(SilwoodDataClean) ...
}
\keyword{datasets}

