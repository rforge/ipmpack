\name{SilwoodDataClean}
\alias{SilwoodDataClean}
\docType{data}
\title{
Silwood Monocarp Dataset
}
\description{
Demographic data of several monocarpic plants from Silwood (UK)
}
\usage{data(SilwoodDataClean)}
\format{
  The format is:
 chr "SilwoodDataClean"
}
\details{
data-frame with headings:
	
	- exactDate: date measurement
	
	- exactDateNext: date next measurement
	
	- id: unique plant id
	
	- Species: species name
	
	- Site: location within Silwood park where plant was measured
	
	- rtcr: root crown diameter, measured with caliper (UNITS)
	
	- rtcrNext: root crown diameter at next census time (UNITS)
	
	- ll: length of longest leaf (UNITS)
	
	- llNext: length of longest leaf at next census (UNITS)
	
	- rosetteDiam: rosette Diameter (UNITS)
	
	- rosetteDiamNext: rosette diameter at next census (UNITS)
	
	- flowered: probability of reproduction (0: individual did not flowered, 1: individual flowered)
	
	- surv: survival (0: dead, 1: survival, NA: not known)
}
\references{
Data used in: Metcalf, C.J.E., Rees, M., Alexander, J.M., Rose, K.E. 2006. Growth-survival trade-offs and allometries in rosette-forming perennials. Funct. Ecol. 20, 217-225.
}
\examples{
data(SilwoodDataClean)
print(head(SilwoodDataClean))
plot(SilwoodDataClean$rtcr,SilwoodDataClean$rtcrNext,
xlab="size now", ylab="size next", pch=19,log="xy")
## maybe str(SilwoodDataClean) ; plot(SilwoodDataClean) ...
}
\keyword{datasets}

